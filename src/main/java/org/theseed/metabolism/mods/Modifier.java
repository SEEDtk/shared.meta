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

    // FIELDS
    /** TRUE if this modifier is active, else FALSE */
    private boolean active;

    /**
     * Construct a new modifier from a parameter line and a metabolic model.
     *
     * @param line		parameter line to parse (optional)
     */
    public Modifier(String line) {
        // The modifier stays on until we turn it off.
        this.active = true;
    }

    /**
     * @return the parameter line for this modifier
     */
    public abstract String getParms();

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

    /**
     * @return TRUE if this modifier should be used, FALSE if it is turned off
     */
    public boolean isActive() {
        return this.active;
    }

    /**
     * Specify whether or not this modifier should be used.
     *
     * @param active 	TRUE to turn on the modifier, FALSE to suppress it
     */
    public void setActive(boolean active) {
        this.active = active;
    }

    /**
     * @return the command code for this modifier
     */
    public abstract String getCommand();

}
