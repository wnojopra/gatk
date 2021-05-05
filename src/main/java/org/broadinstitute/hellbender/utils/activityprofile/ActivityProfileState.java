package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

/**
 * Captures the probability that a specific locus in the genome represents an "active" site containing
 * real variation.
 */
public final class ActivityProfileState {
    private final SimpleInterval loc;
    private double activeProb;
    private final Type resultState;
    private final int size;

    public double isActiveProb() {
        return activeProb;
    }

    /**
     * @return The type of the value returned by {@link #getSize}
     */
    public Type getResultState() {
        return resultState;
    }

    /**
     * @return Numeric value associated with {@link #getResultState}. If {@link #getResultState} is HIGH_QUALITY_SOFT_CLIPS,
     *         this is the number of bp affected by the soft clips
     */
    public int getSize() {
        return size;
    }

    /**
     * The type of the value returned by {@link #getSize}
     */
    public enum Type {
        NONE
    }

    /**
     * Create a new ActivityProfileState at loc with probability of being active of activeProb
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param activeProb the probability of being active (between 0 and 1)
     */
    public ActivityProfileState(final SimpleInterval loc, final double activeProb) {
        this(loc, activeProb, Type.NONE, 0);
    }

    /**
     * Create a new ActivityProfileState at loc with probability of being active of activeProb that maintains some
     * information about the result state
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param activeProb the probability of being active (between 0 and 1)
     */
    public ActivityProfileState(final SimpleInterval loc, final double activeProb, final Type resultState, final int size) {
        Utils.validateArg(loc.size() == 1, "Location for an ActivityProfileState must have to size 1 bp but saw " + loc);
        this.size = ParamUtils.isPositiveOrZero(size, "Size may not be negative");
        this.loc = loc;
        this.activeProb = activeProb;
        this.resultState = resultState;
    }

    /**
     * The offset of state w.r.t. our current region's start location
     * @param regionStartLoc the start of the region, as a Locatable
     * @return the position of this profile relative to the start of this region
     */
    public int getOffset(final Locatable regionStartLoc) {
        Utils.nonNull(regionStartLoc);
        return getLoc().getStart() - regionStartLoc.getStart();
    }

    /**
     * Get the locus associated with the ActivityProfileState
     * @return the locus associated with the ActivityProfileState as a SimpleInterval
     */
    public SimpleInterval getLoc() {
        return loc;
    }

    @Override
    public String toString() {
        return "ActivityProfileState{" +
                "loc=" + loc +
                ", activeProb=" + activeProb +
                ", resultState=" + resultState +
                ", size=" + size +
                '}';
    }
}
