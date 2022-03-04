/**
 *
 */
package org.theseed.metabolism.mods;

import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.theseed.metabolism.Pathway;


/**
 * This pathway filter will reject a pathway immediately if it tries to go
 * through one of a set of prohibited compounds.
 *
 * @author Bruce Parrello
 *
 */
public class AvoidPathwayFilter extends PathwayFilter {

    // FIELDS
    /** IDs for the set of compounds to avoid */
    private Set<String> prohibited;

    public AvoidPathwayFilter(String line) {
        super(line);
        this.prohibited = Set.of(StringUtils.split(line));
    }

    /**
     * Create a pathway filter for avoiding compounds in a specific set
     *
     * @param prohibit	array of prohibited compounds
     */
    public AvoidPathwayFilter(String... prohibit) {
        super(null);
        this.prohibited = new TreeSet<String>();
        for (String compound : prohibit)
            this.prohibited.add(compound);
    }

    @Override
    public boolean isPossible(Pathway path) {
        // We only check the terminus, since this method is called after each
        // reaction is added.
        Pathway.Element terminus = path.getLast();
        return ! (this.prohibited.contains(terminus.getOutput()));
    }

    @Override
    public boolean isGood(Pathway path) {
        // Here we make sure the path is good.
        boolean retVal = true;
        final int n = path.size();
        for (int i = 0; i < n && retVal; i++) {
            Pathway.Element segment = path.getElement(i);
            retVal = ! (this.prohibited.contains(segment.getOutput()));
        }
        return retVal;
    }

    @Override
    public String getParms() {
        return StringUtils.join(this.prohibited, " ");
    }

    @Override
    protected boolean checkEqual(Object other) {
        AvoidPathwayFilter o = (AvoidPathwayFilter) other;
        return (this.prohibited.equals(o.prohibited));
    }

    @Override
    public String getCommand() {
        return "AVOID";
    }

}
