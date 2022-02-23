/**
 *
 */
package org.theseed.metabolism.mods;

import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.theseed.metabolism.Reaction;
import org.theseed.metabolism.Reaction.ActiveDirections;

/**
 * This modifier takes a list of reaction IDs and prevents them from being reversed.
 *
 * @author Bruce Parrello
 *
 */
public class OneWayModifier extends FlowModifier {

    // FIELDS
    /** list of reactions that cannot be reversed */
    private Set<String> reactions;

    /**
     * Construct a one-way modifier from the specified input line.
     *
     * @param line		input line containing reactions to make one-way (space-delimited)
     */
    public OneWayModifier(String line) {
        super(line);
        String[] reacts = StringUtils.split(line);
        this.reactions = Set.of(reacts);
    }

    @Override
    public void setActiveDirections(Reaction reaction) {
        if (this.reactions.contains(reaction.getBiggId()))
            reaction.setActive(ActiveDirections.FORWARD);
    }

    @Override
    public String getParms() {
        return StringUtils.join(" ", reactions);
    }

    @Override
    protected boolean checkEqual(Object other) {
        var otherActual = (OneWayModifier) other;
        return this.reactions.equals(otherActual.reactions);
    }

    @Override
    public int hashCode() {
        final int prime = 29;
        int result = 1;
        result = prime * result + ((this.reactions == null) ? 0 : this.reactions.hashCode());
        return result;
    }


}
