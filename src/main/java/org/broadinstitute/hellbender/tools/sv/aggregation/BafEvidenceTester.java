package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class BafEvidenceTester {

    private static final Median MEDIAN = new Median();
    private static final StandardDeviation STDEV = new StandardDeviation();
    private final KolmogorovSmirnovTest KS_TEST =  new KolmogorovSmirnovTest();
    private static final Random RAND = new Random(42);

    private final int minSnpCarriers;
    private final int minBafCount;
    private final int maxBafCount;
    private final double pSnp;
    private final double pMaxHomozygous;

    public BafEvidenceTester(final int minSnpCarriers, final int minBafCount, final int maxBafCount, final double pSnp, final double pMaxHomozygous) {
        Utils.validateArg(minSnpCarriers > 0, "minSnpCarriers must be positive");
        Utils.validateArg(minBafCount > 0, "minBafCount must be positive");
        Utils.validateArg(pSnp > 0 && pSnp < 1, "pSnp must be a probability on (0, 1)");
        Utils.validateArg(pMaxHomozygous > 0 && pMaxHomozygous < 1, "pMaxHomozygous must be a probability on (0, 1)");
        this.minSnpCarriers = minSnpCarriers;
        this.minBafCount = minBafCount;
        this.maxBafCount = maxBafCount;
        this.pSnp = pSnp;
        this.pMaxHomozygous = pMaxHomozygous;
    }

    public BafResult calculateLogLikelihood(final SVCallRecord record, final List<BafEvidence> evidence, final Set<String> excludedSamples, final int flankSize) {
        if (!record.isSimpleCNV()) {
            return null;
        }
        if (evidence == null || evidence.isEmpty()) {
            return null;
        }
        final Set<String> carrierSamples = Sets.difference(record.getCarrierSampleSet(), excludedSamples);
        final Set<String> allSamples = Sets.difference(record.getAllSamples(), excludedSamples);
        final List<BafEvidence> filteredEvidence = evidence.stream().filter(baf -> allSamples.contains(baf.getSample())).collect(Collectors.toList());

        final List<BafEvidence> innerBaf = new ArrayList<>(filteredEvidence.size());
        final List<BafEvidence> beforeBaf = new ArrayList<>(filteredEvidence.size());
        final List<BafEvidence> afterBaf = new ArrayList<>(filteredEvidence.size());
        final Map<String, SampleStats> sampleStats = allSamples.stream().collect(Collectors.toMap(s -> s, s -> new SampleStats()));
        for (final BafEvidence baf : filteredEvidence) {
            if (baf.getStart() < record.getPositionA()) {
                sampleStats.get(baf.getSample()).beforeHetCount++;
                beforeBaf.add(baf);
            } else if (baf.getStart() >= record.getPositionB()) {
                sampleStats.get(baf.getSample()).afterHetCount++;
                afterBaf.add(baf);
            } else {
                sampleStats.get(baf.getSample()).innerHetCount++;
                innerBaf.add(baf);
            }
        }

        if (record.getType() == StructuralVariantType.DEL) {
            return calculateDeletionTestStatistic(record.getLength(), sampleStats, carrierSamples, flankSize);
        } else {
            if (record.getId().equals("ref_panel_1kg_depth_DUP_chr11_78")) {
                int x = 0;
            }
            if (record.getId().equals("ref_panel_1kg_depth_DUP_chr11_231")) {
                int x = 0;
            }
            return calculateDuplicationTestStatistic(innerBaf, carrierSamples);
        }
    }

    public SVCallRecord applyToRecord(final SVCallRecord record, final BafResult result) {
        Utils.nonNull(record);
        if (result == null) {
            return record;
        }
        final Map<String, Object> attributes = new HashMap<>();
        if (record.getType() == StructuralVariantType.DEL) {
            attributes.put(GATKSVVCFConstants.BAF_HET_RATIO_ATTRIBUTE, result.getDelHetRatio());
        } else {
            final Integer q = EvidenceStatUtils.probToQual(result.getDupMannWhitneyP(), (byte) 99);
            attributes.put(GATKSVVCFConstants.BAF_MWU_QUALITY_ATTRIBUTE, q);
        }
        return SVCallRecordUtils.copyCallWithNewAttributes(record, attributes);
    }

    private BafResult calculateDuplicationTestStatistic(final List<BafEvidence> evidence, final Set<String> carrierSamples) {
        final List<BafEvidence> frequencyFilteredEvidence = new ArrayList<>();
        final Iterator<BafEvidence> iter = evidence.iterator();
        final List<BafEvidence> buffer = new ArrayList<>();
        int pos = -1;
        while (iter.hasNext()) {
            final BafEvidence baf = iter.next();
            if (baf.getStart() != pos) {
                if (buffer.size() >= minSnpCarriers) {
                    frequencyFilteredEvidence.addAll(buffer);
                }
                buffer.clear();
                pos = baf.getStart();
            }
            buffer.add(baf);
        }
        if (buffer.size() >= minSnpCarriers) {
            frequencyFilteredEvidence.addAll(buffer);
        }

        double[] carrierBaf = frequencyFilteredEvidence.stream().filter(baf -> carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).toArray();
        if (carrierBaf.length < minBafCount) {
            return null;
        } else if (carrierBaf.length > maxBafCount) {
            //carrierBaf = shrinkArray(carrierBaf, maxBafCount);
            //carrierBaf = downsampleArray(carrierBaf, maxBafCount);
        }
        double[] nullBaf = frequencyFilteredEvidence.stream().filter(baf -> !carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).toArray();
        if (nullBaf.length < minBafCount) {
            return null;
        } else if (nullBaf.length > maxBafCount) {
            //nullBaf = shrinkArray(nullBaf, maxBafCount);
            //nullBaf = downsampleArray(nullBaf, maxBafCount);
        }
        //final double stat = KS_TEST.kolmogorovSmirnovStatistic(carrierBaf, nullBaf);
        //final double p = KS_TEST.monteCarloP(stat, carrierBaf.length, nullBaf.length, true, 10000);
        //final double p = kolmogorovSmirnovTest(stat, carrierBaf, nullBaf);
        //return BafResult.createDuplicationResult(p, stat);
        final MannWhitneyU test = new MannWhitneyU();
        final MannWhitneyU.Result stat = test.test(carrierBaf, nullBaf, MannWhitneyU.TestType.TWO_SIDED);
        return BafResult.createDuplicationResult(stat.getP());
        //final MannWhitneyUTest test = new MannWhitneyUTest();
        //final double u = (long) carrierBaf.length * nullBaf.length - test.mannWhitneyU(nullBaf, carrierBaf);

        //final double p = calculateAsymptoticPValue(u, carrierBaf.length, nullBaf.length);
        //return BafResult.createDuplicationResult(p, u);
        /*
        final double nullMedian = MEDIAN.evaluate(nullBaf);
        final double carrierMedian = MEDIAN.evaluate(carrierBaf);
        final double nullStd = STDEV.evaluate(nullBaf);
        final double carrierStd = STDEV.evaluate(carrierBaf);
        return (nullMedian - carrierMedian) / Math.sqrt(nullStd * nullStd + carrierStd * carrierStd);
        */
    }

    public double kolmogorovSmirnovTest(final double stat, double[] x, double[] y) {
        final long lengthProduct = (long) x.length * y.length;
        if (lengthProduct < 200) {
            return KS_TEST.exactP(stat, x.length, y.length, true);
        }
        if (lengthProduct < 10000) {
            return KS_TEST.monteCarloP(stat, x.length, y.length, true, 1000000);
        }
        return KS_TEST.approximateP(stat, x.length, y.length);
    }

    private static double[] shrinkArray(final double[] arr, final int size) {
        if (arr.length <= size) {
            return arr;
        }
        final double[] sorted = Arrays.copyOf(arr, arr.length);
        Arrays.sort(sorted);
        final int k = arr.length / size;
        int j = k / 2;
        final double[] out = new double[size];
        for (int i = 0; i < size; i++) {
            out[i] = sorted[j];
            j += k;
        }
        return out;
    }

    private static double[] downsampleArray(final double[] arr, final int size) {
        for (int i = arr.length - 1; i >= arr.length - size; i--) {
            final double tmp = arr[i];
            final int j = RAND.nextInt(arr.length);
            arr[i] = arr[j];
            arr[j] = tmp;
        }
        return Arrays.copyOf(arr, size);
    }


    private BafResult calculateDeletionTestStatistic(final int length,
                                                     final Map<String, SampleStats> sampleStats,
                                                     final Set<String> carrierSamples,
                                                     final int flankSize) {
        //final double pSnp = 0.001; //Math.min(50. / length, 0.0005);
        //int totalInnerCount = 0;
        final List<Double> nullRatios = new ArrayList<>();
        final List<Double> carrierRatios = new ArrayList<>();

        final BinomialDistribution binomialDistributionFlank = new BinomialDistribution(flankSize, pSnp);
        final BinomialDistribution binomialDistributionInner = new BinomialDistribution(length, pSnp);

        for (final Map.Entry<String, SampleStats> entry : sampleStats.entrySet()) {
            final String sample = entry.getKey();
            final SampleStats stats = entry.getValue();
            final double pFlank = binomialDistributionFlank.cumulativeProbability(Math.min(stats.beforeHetCount, stats.afterHetCount));
            final double pInner = binomialDistributionInner.cumulativeProbability(stats.innerHetCount);
            if (!(pInner < pMaxHomozygous && pFlank < pMaxHomozygous)) {
                stats.deletionRatio = Math.log(stats.innerHetCount + 1.);
                if (carrierSamples.contains(sample)) {
                    carrierRatios.add(stats.deletionRatio);
                } else {
                    nullRatios.add(stats.deletionRatio);
                }
            }
            //totalInnerCount += stats.innerHetCount;
        }

        if (carrierRatios.isEmpty() || nullRatios.isEmpty()) {
            return null;
        }
        final double medianCarrier = MEDIAN.evaluate(carrierRatios.stream().mapToDouble(Double::doubleValue).toArray());
        final double medianNull = MEDIAN.evaluate(nullRatios.stream().mapToDouble(Double::doubleValue).toArray());
        return BafResult.createDeletionResult(medianNull - medianCarrier);

        /*
        if (nullRatios.size() <= 10 || totalInnerCount < 10) {
            return null;
        }
        final double minNullRatio = nullRatios.stream().mapToDouble(d -> d).min().getAsDouble();
        final double maxNullRatio = nullRatios.stream().mapToDouble(d -> d).max().getAsDouble();
        if (maxNullRatio - minNullRatio < 0.0001) {
            return null;
        }
        final List<SampleStats> carrierStats = carrierSamples.stream().map(sampleStats::get).filter(s -> s.deletionRatio != null).collect(Collectors.toList());
        if (carrierStats.isEmpty() || carrierStats.size() > nullRatios.size()) {
            return null;
        }
        final double[] carrierRatiosArr = carrierSamples.stream().map(sampleStats::get).filter(s -> s.deletionRatio != null).mapToDouble(s -> s.deletionRatio).toArray();
        return median.evaluate(carrierRatiosArr);
         */
    }

    private static final class SampleStats {
        public int beforeHetCount = 0;
        public int innerHetCount = 0;
        public int afterHetCount = 0;
        public Double deletionRatio = null;
    }

    public static final class BafResult {
        private final Double delHetRatio;
        private final Double dupMannWhitneyP;

        private BafResult(final Double delHetRatio, final Double dupMannWhitneyP) {
            this.delHetRatio = delHetRatio;
            this.dupMannWhitneyP = dupMannWhitneyP;
        }

        private static BafResult createDeletionResult(final double delP) {
            return new BafResult(delP, null);
        }

        private static BafResult createDuplicationResult(final double dupP) {
            return new BafResult(null, dupP);
        }

        public Double getDelHetRatio() {
            return delHetRatio;
        }

        public Double getDupMannWhitneyP() {
            return dupMannWhitneyP;
        }
    }

    protected boolean isRegionOfHomozygosity(final int beforeHetCount, final int innerHetCount, final int afterHetCount, final double length, final double threshold) {
        return innerHetCount < threshold * length &&
                (beforeHetCount  < threshold * length || afterHetCount < threshold * length);
    }

    protected double calculateDeletionRatio(final int beforeHetCount, final int innerHetCount, final int afterHetCount, final double length, final double threshold) {
        final int flankHetCount = Math.min(beforeHetCount, afterHetCount);
        return Math.log10(innerHetCount + length * threshold) - Math.log10(flankHetCount + length * threshold);
    }
}
