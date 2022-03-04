/**
 *
 */
package org.theseed.metabolism.mods;

import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.theseed.metabolism.MetaModel;

/**
 * This modifier specifies a set of compounds as commonly-occurring.  Such compounds can never be in
 * the main line of a pathway, as they are considered ubiquitous.  Certain reports also omit common
 * compounds from the output.
 *
 * @author Bruce Parrello
 *
 */
public class CommonCompoundModifier extends Modifier {

    // FIELDS
    /** set of compounds to be designated as common */
    private Set<String> commons;

    public CommonCompoundModifier(String line) {
        super(line);
        this.commons = Set.of(StringUtils.split(line));
    }

    @Override
    public String getParms() {
        return StringUtils.join(this.commons, " ");
    }

    @Override
    public void updateModel(MetaModel model) {
        model.addCommons(this.commons);
    }

    @Override
    protected boolean checkEqual(Object other) {
        CommonCompoundModifier o = (CommonCompoundModifier) other;
        return (this.commons.equals(o.commons));
    }

    @Override
    public String getCommand() {
        return "COMMONS";
    }

}
