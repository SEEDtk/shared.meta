/**
 *
 */
package org.theseed.metabolism.mods;

import java.util.HashSet;
import java.util.Set;
import org.apache.commons.lang3.StringUtils;
import org.theseed.metabolism.MetaModel;
import org.theseed.metabolism.Pathway;

/**
 * This pathway filter lists the only uncommon compounds that can appear on the input side of a new reaction.  The
 * parameter line is a space-delimited list of exceptions.  If it is empty, then NO uncommon compounds are
 * allowed as inputs other than the reaction's mainline input.
 *
 * @author Bruce Parrello
 *
 */
public class CofactorFilter extends PathwayFilter {

    // FIELDS
    /** set of uncommons to use */
    private Set<String> uncommons;
    /** full set of permissible inputs */
    private Set<String> commons;

    /**
     * Create a new co-factor filter.
     *
     * @param line		space-delimited list of acceptable uncommon inputs
     */
    public CofactorFilter(String line) {
        super(line);
        this.uncommons = Set.of(StringUtils.split(line));
    }

    /**
     * We override updateModel so that we have access to its common-compound list.
     *
     * @param model		model to update
     */
    @Override
    public void updateModel(MetaModel model) {
        // Get all the common compounds.
        this.commons = model.getCommons();
        // Add the permissible uncommons.
        this.commons.addAll(this.uncommons);
        // Add this filter to the model.
        super.updateModel(model);
    }

    @Override
    public boolean isPossible(Pathway path) {
        // Loop through the pathway, accumulating legal compounds.
        Set<String> legal = new HashSet<String>(this.commons);
        var iter = path.iterator();
        boolean retVal = true;
        // Add the global input to the legal input set.
        legal.add(path.getInput());
        // Verify all the path elements.
        while (retVal && iter.hasNext()) {
            var element = iter.next();
            retVal = element.getInputs().stream().allMatch(x -> legal.contains(x));
            // Add the outputs at this stage to the legal list.
            if (retVal)
                legal.addAll(element.getOutputs());
        }
        return retVal;
    }

    @Override
    public boolean isGood(Pathway path) {
        // The isPossible does all our checking for us.
        return true;
    }

    @Override
    public String getParms() {
        return StringUtils.join(this.uncommons, ' ');
    }

    @Override
    protected boolean checkEqual(Object other) {
        CofactorFilter o = (CofactorFilter) other;
        return this.uncommons.equals(o.commons);
    }

    @Override
    public String getCommand() {
        return "COFACTORS";
    }

}
