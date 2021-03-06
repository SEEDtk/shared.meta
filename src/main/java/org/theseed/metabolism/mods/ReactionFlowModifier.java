/**
 *
 */
package org.theseed.metabolism.mods;

import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.theseed.metabolism.Reaction;
import org.theseed.metabolism.Reaction.ActiveDirections;

/**
 * This modifier takes a list of reaction IDs and sets their active directions to a predetermined value.
 *
 * @author Bruce Parrello
 *
 */
public class ReactionFlowModifier extends FlowModifier {

    // FIELDS
    /** list of reactions that cannot be reversed */
    private Set<String> reactions;
    /** direction to specify for all reactions */
    private ActiveDirections dir;

    /**
     * Construct a one-way modifier from the specified input line.
     *
     * @param line			input line containing reactions to make one-way (space-delimited)
     * @param model			target model
     * @param direction		active direction to specify for the reactions
     */
    public ReactionFlowModifier(String line, ActiveDirections direction) {
        super(line);
        String[] reacts = StringUtils.split(line);
        this.reactions = Set.of(reacts);
        this.dir = direction;
    }

    @Override
    public void setActiveDirections(Reaction reaction) {
        if (this.reactions.contains(reaction.getBiggId()))
            reaction.setActive(this.dir);
    }

    @Override
    public String getParms() {
        return StringUtils.join(reactions, ' ');
    }

    @Override
    protected boolean checkEqual(Object other) {
        var otherActual = (ReactionFlowModifier) other;
        return this.dir == otherActual.dir && this.reactions.equals(otherActual.reactions);
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = this.getClass().hashCode();
        result = prime * result + this.dir.ordinal();
        result = prime * result + ((this.reactions == null) ? 0 : this.reactions.hashCode());
        return result;
    }

    @Override
    public String getCommand() {
        return this.dir.getCommand();
    }


}
