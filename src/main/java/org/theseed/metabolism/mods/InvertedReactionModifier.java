/**
 *
 */
package org.theseed.metabolism.mods;

import org.theseed.metabolism.Reaction;

/**
 * @author Bruce Parrello
 *
 */
public class InvertedReactionModifier extends FlowModifier {

    public InvertedReactionModifier(String line) {
        super(line);
        // TODO Auto-generated constructor stub
    }

    @Override
    public void setActiveDirections(Reaction reaction) {
        // TODO code for setActiveDirections

    }

    @Override
    public String getParms() {
        // TODO code for getParms
        return null;
    }

    @Override
    protected boolean checkEqual(Object other) {
        // TODO code for checkEqual
        return false;
    }
    // FIELDS
    // TODO data members for InvertedReactionModifier

    // TODO constructors and methods for InvertedReactionModifier
}
