/**
 *
 */
package org.theseed.metabolism.mods;

import org.theseed.metabolism.MetaModel;
import org.theseed.metabolism.Reaction;

/**
 * This is the base class for all flow modifiers.  A flow modifier determines whether or not the active-
 * direction indicator on a reaction in a metabolic model needs to be updated.   Flow modifiers are
 * generally read from a text file consisting of a command that identifies the type followed by
 * a parameter string.  The required constructor expects the parameter string as input.
 *
 * @author Bruce Parrello
 *
 */
public abstract class FlowModifier extends Modifier {

    /**
     * Construct a flow modifier from an input line.
     *
     * @param line		parameter string to parse
     */
    public FlowModifier(String line) {
        super(line);
    }

    /**
     * Update the active-direction status of a reaction.
     *
     * @param reaction		reaction to check
     */
    public abstract void setActiveDirections(Reaction reaction);

    @Override
    public void updateModel(MetaModel model) {
        for (Reaction react : model.getAllReactions())
            this.setActiveDirections(react);
    }


}
