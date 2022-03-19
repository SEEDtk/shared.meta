/**
 *
 */
package org.theseed.metabolism;

import java.util.Set;

/**
 * A protein rating estimates the effectiveness of a protein in affecting a pathway.  The protein
 * will be marked as either triggering (should be overexpressed) or branching (should be knocked out).
 *
 * The protein is associated with a weight, based on the sum of its weights for all the reactions of
 * interest.  The rating itself starts with a weight of 0, and then then user adds reaction weights
 * to it, modified by considering for the protein's importance to the reaction rule.
 *
 * The rating is ordered by weight, then type (insert before delete), then protein ID.
 *
 * @author Bruce Parrello
 *
 */
public class ProteinRating implements Comparable<ProteinRating>, IReactionSource {

    // FIELDS
    /** ID of the protein */
    private String proteinId;
    /** TRUE for triggering, FALSE for branching */
    private boolean insert;
    /** accumulated weight */
    private double weight;
    /** relevant reaction */
    private Reaction reaction;
    /** TRUE if the reaction is reversed, else FALSE */
    private boolean reversed;
    /** target compound ID */
    private String target;

    /**
     * Create a protein rating.
     *
     * @param protein		ID of the protein
     * @param triggering	TRUE for a triggering weight, FALSE for a branching weight
     */
    public ProteinRating(String protein, boolean triggering) {
        this.proteinId = protein;
        this.insert = triggering;
        this.weight = 0.0;
        this.reaction = null;
        this.reversed = false;
        this.target = null;
    }

    /**
     * Process a weight.
     *
     * @param other		weight to process
     * @param react		relevant reaction
     * @param reversed	TRUE if the reaction is reversed
     * @param targetId	target compound ID
     */
    public void add(double other, Reaction react, boolean reversed, String targetId) {
        if (other > this.weight) {
            this.weight = other;
            this.reaction = react;
            this.reversed = reversed;
            this.target = targetId;
        }
    }

    /**
     * @return the string representation of the protein
     */
    public String getProteinSpec() {
        String retVal = (this.insert ? this.proteinId : "D" + this.proteinId);
        return retVal;
    }

    /**
     * @return the protein ID
     */
    public String getProteinId() {
        return this.proteinId;
    }

    /**
     * @return the weight
     */
    public double getWeight() {
        return this.weight;
    }

    @Override
    public int compareTo(ProteinRating o) {
        int retVal = Double.compare(o.weight, this.weight);
        if (retVal == 0) {
            retVal = Boolean.compare(o.insert, this.insert);
            if (retVal == 0)
                retVal = this.proteinId.compareTo(o.proteinId);
        }
        return retVal;
    }

    /**
     * @return the relevant reaction
     */
    @Override
    public Reaction getReaction() {
        return this.reaction;
    }

    /**
     * @return the target compound ID
     */
    public String getTarget() {
        return this.target;
    }

    /**
     * @return TRUE if the relevant reaction is reversed
     */
    @Override
    public boolean isReversed() {
        return this.reversed;
    }

    @Override
    public Set<String> getSpecial() {
        return Set.of(this.target);
    }

}
