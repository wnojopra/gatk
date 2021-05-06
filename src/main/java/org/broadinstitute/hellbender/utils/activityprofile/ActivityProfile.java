package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Class holding information about per-base activity scores for
 * assembly region traversal
 */
public class ActivityProfile {
    protected final List<ActivityProfileState> stateList;
    protected final double activeProbThreshold;
    protected SAMFileHeader samHeader;

    /**
     * Create a empty ActivityProfile, restricting output to profiles overlapping intervals, if not null
     * @param activeProbThreshold threshold for the probability of a profile state being active
     */
    public ActivityProfile(final double activeProbThreshold, final SAMFileHeader header) {
        stateList = new ArrayList<>();
        this.activeProbThreshold = activeProbThreshold;
        samHeader = header;
    }

    @Override
    public String toString() {
        return "ActivityProfile{start=" + startLoc() + ", stop=" + stopLoc() + '}';
    }

    /**
     * How many profile results are in this profile?
     * @return the number of profile results
     */
    public int size() {
        return stateList.size();
    }

    /**
     * Is this profile empty? (ie., does it contain no ActivityProfileStates?)
     * @return true if the profile is empty (ie., contains no ActivityProfileStates)
     */
    public boolean isEmpty() {
        return stateList.isEmpty();
    }

    /**
     * Get the span of this activity profile, which is from the start of the first state to the stop of the last
     * @return a potentially null SimpleInterval.  Will be null if this profile is empty
     */
    public SimpleInterval getSpan() {
        return isEmpty() ? null : startLoc().spanWith(stopLoc());
    }

    // the start is just the location of the first state in stateList
    private String contig() {
        return stateList.isEmpty() ? null : stateList.get(0).getLoc().getContig();
    }

    // the start is just the location of the first state in stateList
    private SimpleInterval startLoc() {
        return stateList.isEmpty() ? null : stateList.get(0).getLoc();
    }

    // the stop (inclusive) is just the location of the last state in stateList
    private SimpleInterval stopLoc() {
        return stateList.isEmpty() ? null : stateList.get(stateList.size() - 1).getLoc();
    }

    public int getEnd() {
        return stopLoc().getEnd();
    }

    /**
     * Add the next ActivityProfileState to this profile.
     *
     * Must be contiguous with the previously added result, or an IllegalArgumentException will be thrown
     *
     * @param state a well-formed ActivityProfileState result to incorporate into this profile
     */
    public void add(final ActivityProfileState state) {
        Utils.nonNull(state);
        final SimpleInterval loc = state.getLoc();
        Utils.validateArg(stateList.isEmpty() || loc.getStart() == stopLoc().getStart() + 1,
                () -> "Adding non-contiguous state: " + loc + " after " + stopLoc());
        stateList.add(state);
    }

    /**
     * Get the next completed assembly regions from this profile, and remove all states supporting them from this profile
     *
     * Takes the current profile and finds all of the active / inactive from the start of the profile that are
     * ready.  By ready we mean unable to have their probability modified any longer by future additions to the
     * profile.  The regions that are popped off the profile take their states with them, so the start of this
     * profile will always be after the end of the last region returned here.
     *
     * The regions are returned sorted by genomic position.
     *
     * This function may not return anything in the list, if no regions are ready
     *
     * No returned region will be larger than maxRegionSize.
     *
     * @param assemblyRegionExtension the extension value to provide to the constructed regions
     * @param maxRegionSize the maximize size of the returned region
     * @param atEndOfInterval if true, we are at the end of a contig or called interval and may return a region whose end isn't
     *                        sufficiently far from the end of the stateList.
     * @param phasingBuffer
     * @return a non-null list of active regions
     */
    public List<AssemblyRegion> popReadyAssemblyRegions(final int assemblyRegionExtension, final int maxRegionSize, final boolean atEndOfInterval, int phasingBuffer) {
        Utils.validateArg(assemblyRegionExtension >= 0, () -> "assemblyRegionExtension must be >= 0 but got " + assemblyRegionExtension);
        Utils.validateArg( maxRegionSize > 0, () -> "maxRegionSize must be >= 1 but got " + maxRegionSize);

        final List<AssemblyRegion> regions = new ArrayList<>();


        // keep popping regions until we don't have enough states left to proceed and need to wait for more
        while ( !stateList.isEmpty()  && (atEndOfInterval || stateList.size() >= maxRegionSize + 2 * phasingBuffer)) {
            final Pair<Integer, Boolean> endAndActivity = getEndAndActivityOfNextRegion(maxRegionSize, phasingBuffer);
            final int endInclusive = endAndActivity.getLeft();
            final SimpleInterval nextRegionInterval = new SimpleInterval(contig(), startLoc().getStart(), startLoc().getStart() + endInclusive);// we need to create the active region, and clip out the states we're extracting from this profile
            stateList.subList(0, endInclusive + 1).clear();
            regions.add(new AssemblyRegion(nextRegionInterval, endAndActivity.getRight(), assemblyRegionExtension, samHeader));
        }
        return regions;
    }

    // at what index (inclusive) does the next region to pop end and is it active
    private Pair<Integer, Boolean> getEndAndActivityOfNextRegion(final int maxRegionSize, final int phasingBuffer) {
        final int maxLookAhead = maxRegionSize + 2 * phasingBuffer;
        final List<Integer> variants = IntStream.range(0, size()).limit(maxLookAhead).filter(n -> getProb(n) > activeProbThreshold).boxed().collect(Collectors.toList());
        if (variants.isEmpty()) {
            return ImmutablePair.of(Math.min(size(), maxRegionSize), false);
        } else if (variants.get(0) > phasingBuffer) {
            return ImmutablePair.of(Math.min(variants.get(0) - phasingBuffer, maxRegionSize), false);
        } else {
            int maxGapIndex = 0;
            int maxGap = 0;
            for (int n = 0; n < variants.size(); n++) {
                final int variant = variants.get(n);

                if (variant > maxRegionSize) {  // ignore variants in the phasingBuffer padding region
                    break;
                }
                final int nextVariant = (n < variants.size() - 1) ? variants.get(n+1) : size();
                final int gap = nextVariant - variant;

                if (gap > maxGap) {
                    maxGap = gap;
                    maxGapIndex = n;
                }

                if (gap < 2 * phasingBuffer) {
                    return ImmutablePair.of(variant + phasingBuffer, true);
                }
            }

            // if we haven't found anything with a large enough gap, we use the biggest gap
            return ImmutablePair.of(variants.get(maxGapIndex) + maxGap/2, true);
        }
    }

    /**
     * Helper function to get the probability of the state at offset index
     * @param index a valid offset into the state list
     * @return the isActiveProb of the state at index
     */
    private double getProb(final int index) {
        Utils.validIndex(index, stateList.size());

        return stateList.get(index).isActiveProb();
    }
}
