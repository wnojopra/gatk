package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.samples.Trio;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

public final class TransmittedSingleton extends PedigreeAnnotation implements InfoFieldAnnotation {
    protected final Logger warning = LogManager.getLogger(this.getClass());
    private Set<Trio> trios;
    private final int HI_GQ_THRESHOLD = 20;
    private final int HI_DP_THRESHOLD = 20;
    private final double CALL_RATE_THRESHOLD = 0.90;

    public TransmittedSingleton(final Set<Trio> trios) {
        super((Set<String>) null);
        this.trios = Collections.unmodifiableSet(new LinkedHashSet<>(trios));
    }

    public TransmittedSingleton(final GATKPath pedigreeFile){
        super(pedigreeFile);
    }

    public TransmittedSingleton(){
        super((Set<String>) null);
    }

    private Set<Trio> initializeAndGetTrios() {
        if (trios == null) {
            trios = getTrios();
        }
        return trios;
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.TRANSMITTED_SINGLETON, GATKVCFConstants.NON_TRANSMITTED_SINGLETON);
    }
    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(vc);
        Set<Trio> trioSet = initializeAndGetTrios();
        if (!vc.isBiallelic() || trioSet.isEmpty()) {
            return Collections.emptyMap();
        }
        long highQualCalls = vc.getGenotypes().stream().filter(gt -> gt.getGQ() > HI_GQ_THRESHOLD).count();
        if ((double) highQualCalls / vc.getNSamples() < CALL_RATE_THRESHOLD) {
            return Collections.emptyMap();
        }
        final List<String> transmittedSingletonParent = new ArrayList<>();
        final List<String> nonTransmittedSingletonParent = new ArrayList<>();
        for (final Trio trio : trioSet) {
            if (vc.isBiallelic() &&
                    PossibleDeNovo.contextHasTrioGQs(vc, trio)) {
                final boolean childIsHighGQHet = vc.getGenotype(trio.getChildID()).isHet() && vc.getGenotype(trio.getChildID()).getGQ() >= HI_GQ_THRESHOLD;
                final boolean momIsHighGQHet = vc.getGenotype(trio.getMaternalID()).isHet() && vc.getGenotype(trio.getMaternalID()).getGQ() >= HI_GQ_THRESHOLD;
                final boolean dadIsHighGQHet = vc.getGenotype(trio.getPaternalID()).isHet() && vc.getGenotype(trio.getPaternalID()).getGQ() >= HI_GQ_THRESHOLD;

                final boolean momIsHighGQHomRef = vc.getGenotype(trio.getMaternalID()).isHomRef() && vc.getGenotype(trio.getMaternalID()).getGQ() >= HI_GQ_THRESHOLD;
                final boolean dadIsHighGQHomRef = vc.getGenotype(trio.getPaternalID()).isHomRef() && vc.getGenotype(trio.getPaternalID()).getGQ() >= HI_GQ_THRESHOLD;
                final boolean childIsHighGQHomRef = vc.getGenotype(trio.getChildID()).isHomRef() && vc.getGenotype(trio.getChildID()).getGQ() >= HI_GQ_THRESHOLD;

                final boolean childIsHighDepth = vc.getGenotype(trio.getChildID()).getDP() >= HI_DP_THRESHOLD;
                final boolean momIsHighDepth = vc.getGenotype(trio.getChildID()).getDP() >= HI_DP_THRESHOLD;
                final boolean dadIsHighDepth = vc.getGenotype(trio.getChildID()).getDP() >= HI_DP_THRESHOLD;

                if (childIsHighDepth && momIsHighDepth && dadIsHighDepth &&
                        vc.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0) == 2) {
                    if (childIsHighGQHet && momIsHighGQHet && dadIsHighGQHomRef) {
                        transmittedSingletonParent.add(trio.getMaternalID());
                    } else if (childIsHighGQHet && dadIsHighGQHet && momIsHighGQHomRef) {
                        transmittedSingletonParent.add(trio.getPaternalID());
                    }
                }
                //TODO: This only works for trios (not quads or other more complicated family structures that would effect number of singletons for parents or transmission to multiple kids)
                if (childIsHighDepth && momIsHighDepth && dadIsHighDepth &&
                        vc.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0) == 1) {
                    if (childIsHighGQHomRef && momIsHighGQHet && dadIsHighGQHomRef) {
                        nonTransmittedSingletonParent.add(trio.getMaternalID());
                    } else if (childIsHighGQHomRef && dadIsHighGQHet && momIsHighGQHomRef) {
                        nonTransmittedSingletonParent.add(trio.getPaternalID());
                    }
                }
            }
        }
        final Map<String, Object> attributeMap = new LinkedHashMap<>(1);
        if (!transmittedSingletonParent.isEmpty()) {
            attributeMap.put(GATKVCFConstants.TRANSMITTED_SINGLETON, transmittedSingletonParent);
        }
        if (!nonTransmittedSingletonParent.isEmpty()) {
            attributeMap.put(GATKVCFConstants.NON_TRANSMITTED_SINGLETON, nonTransmittedSingletonParent);
        }
        return attributeMap;
    }
}
