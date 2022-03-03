/**
 *
 */
package org.theseed.metabolism.mods;

import org.theseed.metabolism.MetaModel;

/**
 * A modifier is loaded from a modification file and provides a temporary change to the way the
 * model works.  Some modifiers affect the active direction.  Others change the way compounds are
 * treated or place additional requirements on a pathway.
 *
 * @author Bruce Parrello
 *
 */
public abstract class Modifier {

    /**
     * Construct a new modifier from a parameter line and a metabolic model.
     *
     * @param line		parameter line to parse (optional)
     */
    public Modifier(String line) { }

    /**
     * @return the parameter line for this modifier
     */
    protected abstract String getParms();

    /**
     * Update the model with this modifier.
     *
     * @param model		model to update
     */
    public abstract void updateModel(MetaModel model);

    /**
     * @return TRUE if this modifier is equivalent to the other
     *
     * @param other		other modifier to check (which must be the same type as this one)
     */
    protected abstract boolean checkEqual(Object other);


    @Override
    public final boolean equals(Object other) {
        boolean retVal = false;
        if (other != null && other.getClass() == this.getClass()) {
            // Here we are the same type.  Use the subclass's comparison tool.
            retVal = this.checkEqual(other);
        }
        return retVal;
    }



}
