package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

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

    // --------------------------------------------------------------------------------
    //
    // routines to get active regions from the profile
    //
    // --------------------------------------------------------------------------------

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
     * @param minRegionSize the minimum region size, in the case where we have to cut up regions that are too large
     * @param maxRegionSize the maximize size of the returned region
     * @param atEndOfInterval if true, we are at the end of a contig or called interval and may return a region whose end isn't
     *                        sufficiently far from the end of the stateList.
     * @return a non-null list of active regions
     */
    public List<AssemblyRegion> popReadyAssemblyRegions( final int assemblyRegionExtension, final int minRegionSize, final int maxRegionSize, final boolean atEndOfInterval ) {
        Utils.validateArg(assemblyRegionExtension >= 0, () -> "assemblyRegionExtension must be >= 0 but got " + assemblyRegionExtension);
        Utils.validateArg( minRegionSize > 0, () -> "minRegionSize must be >= 1 but got " + minRegionSize);
        Utils.validateArg( maxRegionSize > 0, () -> "maxRegionSize must be >= 1 but got " + maxRegionSize);

        final List<AssemblyRegion> regions = new ArrayList<>();

        // keep popping regions until we don't have enough states left to proceed and need to wait for more
        while ( !(stateList.isEmpty() || (atEndOfInterval && stateList.size() < maxRegionSize) )) {
            final ActivityProfileState first = stateList.get(0);
            final boolean isActiveRegion = first.isActiveProb() > activeProbThreshold;
            final int sizeOfNextRegion = findEndOfRegion(isActiveRegion, minRegionSize, maxRegionSize, atEndOfInterval);
            final SimpleInterval nextRegionInterval = new SimpleInterval(first.getLoc().getContig(), first.getLoc().getStart(), first.getLoc().getStart() + sizeOfNextRegion);// we need to create the active region, and clip out the states we're extracting from this profile

            stateList.subList(0, sizeOfNextRegion + 1).clear();
            regions.add(new AssemblyRegion(nextRegionInterval, isActiveRegion, assemblyRegionExtension, samHeader));
        }
        return regions;
    }

    /**
     * Find the end of the current region, returning the index into the element isActive element, or -1 if the region isn't done
     *
     * The current region is defined from the start of the stateList, looking for elements that have the same isActiveRegion
     * flag (i.e., if isActiveRegion is true we are looking for states with isActiveProb > threshold, or alternatively
     * for states < threshold).  The maximize size of the returned region is maxRegionSize.  If forceConversion is
     * true, then we'll return the region end even if this isn't safely beyond the max prob propagation distance.
     *
     * Note that if isActiveRegion is true, and we can construct an assembly region > maxRegionSize in bp, we
     * find the further local minimum within that max region, and cut the region there, under the constraint
     * that the resulting region must be at least minRegionSize in bp.
     *
     * @param isActiveRegion is the region we're looking for an active region or inactive region?
     * @param minRegionSize the minimum region size, in the case where we have to cut up regions that are too large
     * @param maxRegionSize the maximize size of the returned region
     * @param forceConversion if true, we'll return a region whose end isn't sufficiently far from the end of the
     *                        stateList.  Used to close out the assembly region when we've hit some kind of end (such
     *                        as the end of the contig)
     * @return the index into stateList of the last element of this region, or -1 if it cannot be found
     */
    private int findEndOfRegion(final boolean isActiveRegion, final int minRegionSize, final int maxRegionSize, final boolean forceConversion) {
        int endOfActiveRegion = findFirstActivityBoundary(isActiveRegion, maxRegionSize);

        if ( isActiveRegion && endOfActiveRegion == maxRegionSize ) {
            // we've run to the end of the region, let's find a good place to cut
            endOfActiveRegion = findBestCutSite(endOfActiveRegion, minRegionSize);
        }

        // we're one past the end, so i must be decremented
        return endOfActiveRegion - 1;
    }

    /**
     * Find the the local minimum within 0 - endOfActiveRegion where we should divide region
     *
     * This algorithm finds the global minimum probability state within the region [minRegionSize, endOfActiveRegion)
     * (exclusive of endOfActiveRegion), and returns the state index of that state.
     * that it
     *
     * @param endOfActiveRegion the last state of the current active region (exclusive)
     * @param minRegionSize the minimum of the left-most region, after cutting
     * @return the index of state after the cut site (just like endOfActiveRegion)
     */
    private int findBestCutSite(final int endOfActiveRegion, final int minRegionSize) {
        Utils.validateArg(endOfActiveRegion >= minRegionSize, "endOfActiveRegion must be >= minRegionSize");
        Utils.validateArg(minRegionSize >= 0, "minRegionSize must be >= 0");

        int minI = endOfActiveRegion - 1;
        double minP = Double.MAX_VALUE;

        for ( int i = minI; i >= minRegionSize - 1; i-- ) {
            double cur = getProb(i);
            if ( cur < minP && isMinimum(i) ) {
                minP = cur;
                minI = i;
            }
        }

        return minI + 1;
    }

    /**
     * Find the first index into the state list where the state is considered ! isActiveRegion
     *
     * Note that each state has a probability of being active, and this function thresholds that
     * value on activeProbThreshold, coloring each state as active or inactive.  Finds the
     * largest contiguous stretch of states starting at the first state (index 0) with the same isActive
     * state as isActiveRegion.  If the entire state list has the same isActive value, then returns
     * maxRegionSize
     *
     * @param isActiveRegion are we looking for a stretch of active states, or inactive ones?
     * @param maxRegionSize don't look for a boundary that would yield a region of size > maxRegionSize
     * @return the index of the first state in the state list with isActive value != isActiveRegion, or maxRegionSize
     *         if no such element exists
     */
    private int findFirstActivityBoundary(final boolean isActiveRegion, final int maxRegionSize) {
        Utils.validateArg(maxRegionSize > 0, "maxRegionSize must be > 0");

        final int nStates = stateList.size();
        int endOfActiveRegion = 0;

        while ( endOfActiveRegion < nStates && endOfActiveRegion < maxRegionSize ) {
            if ( getProb(endOfActiveRegion) > activeProbThreshold != isActiveRegion ) {
                break;
            }
            endOfActiveRegion++;
        }

        return endOfActiveRegion;
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

    /**
     * Is the probability at index in a local minimum?
     *
     * Checks that the probability at index is <= both the probabilities to either side.
     * Returns false if index is at the end or the start of the state list.
     *
     * @param index the index of the state we want to test
     * @return true if prob at state is a minimum, false otherwise
     */
    private boolean isMinimum(final int index) {
        Utils.validIndex(index, stateList.size());

        if ( index == stateList.size() - 1 ) {
            // we cannot be at a minimum if the current position is the last in the state list
            return false;
        }
        else if ( index < 1 ) {
            // we cannot be at a minimum if the current position is the first or second
            return false;
        }
        else {
            final double indexP = getProb(index);
            return indexP <= getProb(index+1) && indexP < getProb(index-1);
        }
    }
}
