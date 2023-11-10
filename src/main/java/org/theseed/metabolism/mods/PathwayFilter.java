/**
 *
 */
package org.theseed.metabolism.mods;

import org.theseed.basic.ParseFailureException;
import org.theseed.metabolism.MetaModel;
import org.theseed.metabolism.Pathway;

/**
 * A pathway filter contains various criteria for a pathway.  The basic pathway
 * search is to find a path from compound A to compound B.  The filter will accept
 * or reject paths according to various ancillary criteria.
 *
 * @author Bruce Parrello
 */
public abstract class PathwayFilter extends Modifier {

    /**
     * This enumeration indicates the types of filters.
     */
    public static enum Type {
        /** include certain reactions in the pathway */
        REACTIONS {
            @Override
            public PathwayFilter create(String line) {
                return new IncludePathwayFilter(line);
            }

        },
        /** avoid certain compounds in the pathway */
        AVOID {
            @Override
            public PathwayFilter create(String line) {
                return new AvoidPathwayFilter(line);
            }

        };

        /**
         * @return a pathway filter of the specified type
         *
         * @param line		parameter line
         *
         * @throws ParseFailureException
         */
        public abstract PathwayFilter create(String line) throws ParseFailureException;

    }

    public PathwayFilter(String line) {
        super(line);
    }

    /**
     * This method is called at each stage of the search.  For avoid-type
     * filters, it allows the filter to cut off a pathway as soon as it
     * crosses into a region to be avoided.  (Return FALSE to do this.)
     *
     * @param path		pathway to check
     *
     * @return TRUE if this pathway should be kept
     */
    public abstract boolean isPossible(Pathway path);

    /**
     * This method is called when a pathway has reached its goal.  For
     * include-type filters, it allows the filter to reject a pathway if
     * it does not include everything needed.
     *
     * @param path		pathway to check
     *
     * @return TRUE if this pathway is good
     */
    public abstract boolean isGood(Pathway path);

    @Override
    public void updateModel(MetaModel model) {
        model.addFilter(this);
    }


}
