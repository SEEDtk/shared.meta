/**
 *
 */
package org.theseed.metabolism;

import java.util.HashMap;
import java.util.Map;

/**
 * This object rates a compound according to its importance to a pathway.  The rating
 * includes the total stoichiometry of the compound's use as an input, and the number
 * of reactions preceding the compound's use in the pathway.
 *
 * @author Bruce Parrello
 *
 */
public class CompoundRating implements Comparable<CompoundRating> {

    // FIELDS
    /** compound ID */
    private String compoundId;
    /** total stoichiometry */
    private int usage;
    /** number of preceding pathway reactions */
    private int effort;
    /** TRUE if the compound is common, else FALSE */
    private boolean common;
    /** TRUE if the compound is the goal of the pathway, else FALSE */
    private boolean goal;
    /** number of times compound is produced as an output */
    private int output;
    /** weighted sum of efforts for each stoichiometric occurrence */
    private double weight;
    /** weight for being the goal */
    private static double GOAL_WEIGHT = 1.1;
    /** weight for compound efforts */
    private static double EFFORT_WEIGHT = 1.0;
    /** null rating for a weightless compound */
    public static final CompoundRating IRRELEVANT = new CompoundRating("*", true, false);

    /**
     * Construct a blank rating for a compound in a pathway.
     *
     *  @param compound		compound of interest
     *  @param common		TRUE if the compound is common, else FALSE
     *  @param goal			TRUE if the compound is the goal, else FALSE
     */
    protected CompoundRating(String compound, boolean common, boolean goal) {
        // Initialize the fields.
        this.compoundId = compound;
        this.common = common;
        this.usage = 0;
        this.effort = 0;
        this.output = 0;
        this.weight = 0.0;
        this.goal = goal;
    }

    /**
     * @return a sorted list of the compound ratings for all compounds in the specified path
     *
     * @param path		pathway whose compounds are desired
     * @param model		underlying model for the path
     */
    public static Map<String, CompoundRating> getRatingMap(Pathway path, MetaModel model) {
        var retVal = new HashMap<String, CompoundRating>(path.size() * 4);
        // Get the common-compound set.
        var commons = model.getCommons();
        // First, we put in the goal compound.
        String goal = path.getOutput();
        CompoundRating goalRating = new CompoundRating(goal, commons.contains(goal), true);
        goalRating.weight = path.size() * EFFORT_WEIGHT * GOAL_WEIGHT;
        retVal.put(goal, goalRating);
        // Now, loop through the path elements, processing each compound found.
        int position = 0;
        for (Pathway.Element element : path) {
            // The position indicates the effort to reach the compound.  It is used to give
            // higher weight to compounds used later in the chain.
            position++;
            // Remember whether or not the reaction is reversed.
            Reaction react = element.getReaction();
            boolean reversed = element.isReversed();
            // Loop through the compounds in this reaction.
            for (Reaction.Stoich stoich : react.getMetabolites()) {
                // Get the compound rating for this compound.  Note that if we have to
                // create the rating, it cannot be the goal compound, which was added first.
                CompoundRating rating = retVal.computeIfAbsent(stoich.getMetabolite(),
                        x -> new CompoundRating(x, false, false));
                int coeff = stoich.getCoeff();
                // Now we add this reaction element.  If it is an output, we only increment the
                // output count.
                if (stoich.isProduct() != reversed)
                    rating.output += coeff;
                else {
                    // Here we have an input count, and we need to count lots more stuff.
                    rating.weight += EFFORT_WEIGHT * position * coeff;
                    if (rating.effort < position) rating.effort = position;
                    rating.usage += coeff;
                }
            }
        }
        // Sort and return the compounds.
        return retVal;
    }

    @Override
    public int compareTo(CompoundRating o) {
        // The goal always wins.
        int retVal = Boolean.compare(o.goal, this.goal);
        if (retVal == 0) {
            // Uncommon is more important than common.
            retVal = Boolean.compare(this.common, o.common);
            if (retVal == 0) {
                // Compare by weight.  A higher weight is more important.
                retVal = Double.compare(o.weight, this.weight);
                if (retVal == 0) {
                    // Compare by usage.  A higher usage count is more important.
                    retVal = o.usage - this.usage;
                    if (retVal == 0) {
                        // Compare by output.  A lower output count is more important.
                        retVal = this.output - o.output;
                        if (retVal == 0) {
                            // Compare by effort.  A higher effort is more important.
                            retVal = o.effort - this.effort;
                            if (retVal == 0) {
                                // Finally, sort by ID.
                                retVal = this.compoundId.compareTo(o.compoundId);
                            }
                        }
                    }
                }
            }
        }
        return retVal;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.compoundId == null) ? 0 : this.compoundId.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof CompoundRating)) {
            return false;
        }
        CompoundRating other = (CompoundRating) obj;
        if (this.compoundId == null) {
            if (other.compoundId != null) {
                return false;
            }
        } else if (!this.compoundId.equals(other.compoundId)) {
            return false;
        }
        return true;
    }

    /**
     * @return the rated compound ID
     */
    public String getCompoundId() {
        return this.compoundId;
    }

    /**
     * @return the number of units of the compound used by the pathway
     */
    public int getUsage() {
        return this.usage;
    }

    /**
     * @return the effort count (minimum number of pathway steps before the compound is required)
     */
    public int getEffort() {
        return this.effort;
    }

    /**
     * @return TRUE if this compound is common
     */
    public boolean isCommon() {
        return this.common;
    }

    /**
     * @return the number of units of the compound output by the pathway
     */
    public int getOutput() {
        return this.output;
    }

    /**
     * @return the weight of the compound
     */
    public double getWeight() {
        return this.weight;
    }

    /**
     * @return TRUE if this compound is the goal
     */
    public Object isGoal() {
        return this.goal;
    }

}
