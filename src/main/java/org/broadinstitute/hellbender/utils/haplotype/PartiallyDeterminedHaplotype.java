package org.broadinstitute.hellbender.utils.haplotype;

import com.google.common.annotations.VisibleForTesting;
import com.netflix.servo.util.Objects;
import htsjdk.samtools.Cigar;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

public final class PartiallyDeterminedHaplotype extends Haplotype {
    private static final long serialVersionUID = 1L;

    //COMMENTS DIRECTLY RIPPED OUT OF DRAGEN CODE:
    // This is the array which gives information if the base in the haplotype is
    // resolved or not.
    // 8-bit field
    //  0   0   0   0   0   0   0   0
    // SNP  DELETE  A   C   G   T   N
    // If SNP is true, then bits 3-7 are active indicating which bases are also possible
    // Bits 1-2 indicate delete states
    // 01 - Delete region start
    // 10 - Delete region stops
    // 11 - Delete has length 1, start and stop position is the same
    // All 8 bits 0 indicate either a resolved position or if w
    final public static byte SNP = 1;
    final static byte DEL_START = 2;
    final static byte DEL_END = 4;
    final static byte A = 8;
    final static byte C = 16;
    final static byte G = 32;
    final static byte T = 64;
    final static byte N = (byte) 128;

    public byte[] getAlternateBases() {
        return alternateBases;
    }

    private byte[] alternateBases;
    private List<VariantContext> constituentEvents;
    private VariantContext importantVariant;

    /**
     * Compares two haplotypes first by their lengths and then by lexicographic order of their bases.
     */
    public static final Comparator<Haplotype> SIZE_AND_BASE_ORDER =
            Comparator.comparingInt((Haplotype hap) -> hap.getBases().length)
                    .thenComparing(Allele::getBaseString);

    /**
     *
     * @param base
     * @param isRef
     * @param pdBytes
     * @param constituentEvents
     */
    public PartiallyDeterminedHaplotype(final Haplotype base, boolean isRef, byte[] pdBytes, List<VariantContext> constituentEvents, VariantContext eventWithVariant, Cigar cigar) {
        super(base.getBases(), isRef, base.getAlignmentStartHapwrtRef(), cigar);
        this.setGenomeLocation(base.getGenomeLocation());
        this.alternateBases = pdBytes;
        this.constituentEvents = constituentEvents;
        this.importantVariant = eventWithVariant;
    }


    @VisibleForTesting
    public static byte[] getPDBytesForHaplotypes(Allele refAllele, Allele altAllele) {
        // Asserting we are either indel OR SNP
        byte[] output = new byte[altAllele.length() == refAllele.length() ? refAllele.length() : refAllele.length()-1 ];
        // SNP case
        if (altAllele.length() == refAllele.length()){
            output[0] += SNP;
            switch (altAllele.getBases()[0]) {
                case 'A':
                    output[0] += A;
                    break;
                case 'C':
                    output[0] += C;
                    break;
                case 'T':
                    output[0] += T;
                    break;
                case 'G':
                    output[0] += G;
                    break;
                default:
                    throw new RuntimeException("Found unexpected base in alt alleles");
            }

        // DEL case
        } else {
            output[0] += DEL_START;
            output[output.length - 1] += DEL_END; //TODO guardrail this
        }
        return output;
    }


    @Override
    public byte[] getDisplayBases() {
        byte[] bytes = getBases().clone();
        for (int i = 0; i < bytes.length; i++) {
            if ((alternateBases[i] & ~3) != 0) { // Don't change the byte if
                bytes[i] = 'P';
            }
        }
        return bytes;
    }

    @Override
    public String toString() {
        String output = "HapLen:"+length() +", "+new String(getDisplayBases());
        output = output + "\nUnresolved Bases["+alternateBases.length+"] "+Arrays.toString(alternateBases);
        return output + "\n"+getCigar().toString()+" "+ constituentEvents.stream()
                .map(v ->(v==this.importantVariant?"*":"")+getDRAGENDebugVariantContextString((int)getStartPosition()).apply(v) )
                .collect(Collectors.joining("->"));
    }

    public static Function<VariantContext, String> getDRAGENDebugVariantContextString(final int offset) {
        return v -> "(" + Integer.toString(v.getStart() - offset + (v.isSimpleDeletion()? 1 : 0))+ (v.isSimpleInsertion()?".5":"") + ",Rlen=" + v.getLengthOnReference()+"," + v.getAlternateAlleles() + ")";
    }

    // NOTE, we don't have to trim this because we form these after trimming

    @Override
    public boolean equals( final Object h ) {
        return h instanceof PartiallyDeterminedHaplotype
                && getUniquenessValue() == ((Haplotype) h).getUniquenessValue()
                && isReference() == ((Haplotype) h).isReference()
                && Arrays.equals(getBases(), ((Haplotype) h).getBases())
                && Arrays.equals(alternateBases, ((PartiallyDeterminedHaplotype) h).alternateBases);
    }

    @Override
    public int hashCode() {
        // TODO is this good enough for a hash function?
        return Objects.hash(Arrays.hashCode(getBases()),Arrays.hashCode(alternateBases));
    }

}
