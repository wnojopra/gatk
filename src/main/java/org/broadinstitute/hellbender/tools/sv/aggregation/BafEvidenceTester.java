package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class BafEvidenceTester {

    private static final Median MEDIAN = new Median();

    private final int minSnpCarriers;
    private final double pSnp;
    private final double pMaxHomozygous;

    public BafEvidenceTester(final int minSnpCarriers, final double pSnp, final double pMaxHomozygous) {
        this.minSnpCarriers = minSnpCarriers;
        this.pSnp = pSnp;
        this.pMaxHomozygous = pMaxHomozygous;
    }

    public Double calculateLogLikelihood(final SVCallRecord record, final List<BafEvidence> evidence, final Set<String> excludedSamples, final int flankSize) {
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
            return calculateDuplicationTestStatistic(innerBaf, carrierSamples);
        }
    }

    public SVCallRecord applyToRecord(final SVCallRecord record, final Double result) {
        Utils.nonNull(record);
        if (result == null) {
            return record;
        }
        final Map<String, Object> attributes = new HashMap<>();
        if (record.getType() == StructuralVariantType.DEL) {
            attributes.put(GATKSVVCFConstants.BAF_STAT_DEL_ATTRIBUTE, result);
        } else {
            attributes.put(GATKSVVCFConstants.BAF_STAT_DUP_ATTRIBUTE, result);
        }
        return SVCallRecordUtils.copyCallWithNewAttributes(record, attributes);
    }

    private Double calculateDuplicationTestStatistic(final List<BafEvidence> evidence, final Set<String> carrierSamples) {

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

        final double[] nullBaf = frequencyFilteredEvidence.stream().filter(baf -> !carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).map(d -> Math.min(d, 1.0 - d)).toArray();
        if (nullBaf.length == 0) {
            return null;
        }
        final double[] carrierBaf = frequencyFilteredEvidence.stream().filter(baf -> carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).map(d -> Math.min(d, 1.0 - d)).toArray();
        if (carrierBaf.length == 0) {
            return null;
        }
        final double meanNull = MEDIAN.evaluate(nullBaf);
        final double meanCarrier = MEDIAN.evaluate(carrierBaf);
        return meanNull - meanCarrier;
    }

    private Double calculateDeletionTestStatistic(final int length,
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
        return medianNull - medianCarrier;

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

    protected boolean isRegionOfHomozygosity(final int beforeHetCount, final int innerHetCount, final int afterHetCount, final double length, final double threshold) {
        return innerHetCount < threshold * length &&
                (beforeHetCount  < threshold * length || afterHetCount < threshold * length);
    }

    protected double calculateDeletionRatio(final int beforeHetCount, final int innerHetCount, final int afterHetCount, final double length, final double threshold) {
        final int flankHetCount = Math.min(beforeHetCount, afterHetCount);
        return Math.log10(innerHetCount + length * threshold) - Math.log10(flankHetCount + length * threshold);
    }
}
