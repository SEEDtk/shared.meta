/**
 *
 */
package org.theseed.metabolism.mods;

import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.theseed.metabolism.Pathway;

/**
 * This pathway filter will only accept a path if it contains certain reactions.
 *
 * @author Bruce Parrello
 *
 */
public class IncludePathwayFilter extends PathwayFilter {

    // FIELDS
    /** set of required reactions */
    private Set<String> required;

    public IncludePathwayFilter(String line) {
        super(line);
        this.required = Set.of(StringUtils.split(line));
    }

    /**
     * Create a pathway filter that requires the specified reactions.
     *
     * @param includes	array of bigg IDs for the reactions to include
     *
     */
    public IncludePathwayFilter(String... includes) {
        super(null);
        this.required = new TreeSet<String>();
        for (String reaction : includes)
            this.required.add(reaction);
    }

    @Override
    public boolean isPossible(Pathway path) {
        // We only make the inclusion check at the end.
        return true;
    }

    @Override
    public boolean isGood(Pathway path) {
        // We take advantage of the fact that a pathway is not allowed to have
        // the same reaction twice.  If we find each reaction once, then we pass
        // the pathway.
        int found = 0;
        for (Pathway.Element part : path) {
            if (required.contains(part.getReaction().getBiggId()))
                found++;
        }
        return (found >= required.size());
    }

    @Override
    public String getParms() {
        return StringUtils.join(required, " ");
    }

    @Override
    protected boolean checkEqual(Object other) {
        IncludePathwayFilter o = (IncludePathwayFilter) other;
        return this.required.equals(o.required);
    }

    @Override
    public String getCommand() {
        return "INCLUDE";
    }

}
