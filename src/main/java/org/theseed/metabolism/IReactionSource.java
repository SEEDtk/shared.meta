/**
 *
 */
package org.theseed.metabolism;

import java.util.Set;

/**
 * This is an interface used to facilitate the display of reaction formulas in various media.
 *
 * @author Bruce Parrello
 *
 */
public interface IReactionSource {

    /**
     * @return the reaction whose formula is to be displayed
     */
    public Reaction getReaction();

    /**
     * @return TRUE if the reaction is reversed, else FALSE
     */
    public boolean isReversed();

    /**
     * @return the set of IDs for compounds that should be boldfaced
     */
    public Set<String> getSpecial();

}
