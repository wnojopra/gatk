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
    private final int hi_GQ_threshold = 20;

    public TransmittedSingleton(final Set<Trio> trios, final double minGenotypeQualityP) {
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
        return Arrays.asList(GATKVCFConstants.TRANSMITTED_SINGLETON);
    }
    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(vc);
        Set<Trio> trioSet = initializeAndGetTrios();
        if (!vc.isBiallelic() || trioSet.isEmpty() || (vc.isBiallelic() && vc.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0) != 2)) {
            return Collections.emptyMap();
        }
        final List<String> transmittedSingletonChildren = new ArrayList<>();
        for (final Trio trio : trioSet) {
            if (vc.isBiallelic() &&
                    PossibleDeNovo.contextHasTrioGQs(vc, trio)) {
                final boolean childIsHighGQHet = vc.getGenotype(trio.getChildID()).isHet() && vc.getGenotype(trio.getChildID()).getGQ() > hi_GQ_threshold;
                final boolean momIsHighGQHet = vc.getGenotype(trio.getMaternalID()).isHet() && vc.getGenotype(trio.getMaternalID()).getGQ() > hi_GQ_threshold;
                final boolean dadIsHighGQHet = vc.getGenotype(trio.getPaternalID()).isHet() && vc.getGenotype(trio.getPaternalID()).getGQ() > hi_GQ_threshold;

                if ((childIsHighGQHet && momIsHighGQHet) || (childIsHighGQHet && dadIsHighGQHet)) {
                    transmittedSingletonChildren.add(trio.getChildID());
                }
            }
        }
        final Map<String, Object> attributeMap = new LinkedHashMap<>(1);
        if (!transmittedSingletonChildren.isEmpty()) {
            attributeMap.put(GATKVCFConstants.TRANSMITTED_SINGLETON, transmittedSingletonChildren);
        }
        return attributeMap;
    }
}
