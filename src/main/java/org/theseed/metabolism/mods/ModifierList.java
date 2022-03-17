/**
 *
 */
package org.theseed.metabolism.mods;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.metabolism.MetaModel;
import org.theseed.metabolism.Reaction.ActiveDirections;

/**
 * This object represents a list of modifiers.  The modifiers are applied to a metabolic model to
 * to suppress directional flow of certain reactions and to otherwise modify construction of pathways.
 * The modifier list is stored in a tab-delimited flat file with headers. The command is in the
 * first column and the parameter string in the second.  A leading pound sign (#) means the modifier
 * should start deactivated.
 *
 * @author Bruce Parrello
 *
 */
public class ModifierList implements Iterable<Modifier> {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ModifierList.class);
    /** list of (command, flow modifier) pairs */
    private List<Modifier> modifiers;

    /**
     * This enum describes the various flow modifier types.
     */
    public static enum Command {
        /** insure specified compounds only appear on the input side of a reaction */
        FORWARD {
            @Override
            public Modifier create(String line) {
                return new ForwardOnlyModifier(line);
            }
        },
        /** suppress one or more reactions */
        SUPPRESS {
            @Override
            public Modifier create(String line) {
                return new ReactionFlowModifier(line, ActiveDirections.NEITHER);
            }
        },
        /** prevent reversing of one or more reactions */
        ONEWAY {
            @Override
            public Modifier create(String line) {
                return new ReactionFlowModifier(line, ActiveDirections.FORWARD);
            }
        },
        /** require reversing of one or more reactions */
        INVERTED {
            @Override
            public Modifier create(String line) {
                return new ReactionFlowModifier(line, ActiveDirections.REVERSE);
            }
        },
        /** avoid one or more compounds in the main pathway */
        AVOID {
            @Override
            public Modifier create(String line) {
                return new AvoidPathwayFilter(line);
            }
        },
        /** allow only certain uncommon cofactors */
        COFACTORS {
            @Override
            public Modifier create(String line) {
                return new CofactorFilter(line);
            }
        },
        /** specify one or more compounds as common */
        COMMONS {
            @Override
            public Modifier create(String line) {
                return new CommonCompoundModifier(line);
            }
        };

        /**
         * @return a flow modifier of this type
         *
         * @param line		parameter string
         */
        public abstract Modifier create(String line);

    }

    /**
     * Create a modifier list from an input file stream.
     *
     * @param tabReader		tabbed line reader stream containing the modifier specs
     *
     * @throws IOException
     */
    public ModifierList(TabbedLineReader tabReader) throws IOException {
        this.readModifiers(tabReader);
    }

    /**
     * Create a modifier list from an input file.
     *
     * @param inFile		tab-delimited file containing modifier commands and parameters
     */
    public ModifierList(File inFile) throws IOException {
        try (TabbedLineReader inStream = new TabbedLineReader(inFile)) {
            this.readModifiers(inStream);
        }
    }

    /**
     * Create an immutable, empty modifier list.
     */
    public ModifierList() {
        this.modifiers = Collections.emptyList();
    }

    /**
     * Create a modifier list backed by a list of Modifier objects.  Changes to the list
     * will be reflected automatically in this object.
     *
     * @param list	underlying list of modifiers to use
     */
    public ModifierList(List<Modifier> list) {
        this.modifiers = list;
    }

    /**
     * Read a modifier list from a tab-delimited input stream.
     *
     * @param tabReader		tabbed line reader stream containing the modifier specs
     *
     * @throws IOException
     */
    private void readModifiers(TabbedLineReader tabReader) throws IOException {
        this.modifiers = new ArrayList<>();
        for (TabbedLineReader.Line line : tabReader) {
            // We convert the command to upper case so it matches the enum name.
            String command = line.get(0).toUpperCase();
            String parms = line.get(1);
            boolean active = true;
            if (command.startsWith("#")) {
                active = false;
                command = command.substring(1);
            }
            this.addModifier(command, parms, active);
        }
    }

    /**
     * Add a flow modifier build from a command and a parameter string.
     *
     * @param command	command (all upper case)
     * @param parms		parameter string
     * @param active	TRUE if the modifier should be active, else FALSE
     *
     * @throws IOException
     */
    private void addModifier(String command, String parms, boolean active) throws IOException {
        Command commandCode;
        // An invalid command code is rethrown as an IO exception.
        try {
            commandCode = Command.valueOf(command);
        } catch (IllegalArgumentException e) {
            throw new IOException("Invalid flow modifier command code \"" + command + "\".");
        }
        // Build the modifier from the command and parameter string.
        Modifier modifier = commandCode.create(parms);
        modifier.setActive(active);
        this.modifiers.add(modifier);
    }

    /**
     * Process a metabolic model and update the active directions for all of its reactions.
     *
     * @param model		metabolic model to update
     *
     * @return the number of modifications made
     */
    public void apply(MetaModel model) {
        // Clean the model.
        model.clearMods();
        // Apply the modifiers.
        for (Modifier mod : this.modifiers) {
            if (mod.isActive())
                mod.updateModel(model);
        }
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.modifiers == null) ? 0 : this.modifiers.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof ModifierList)) {
            return false;
        }
        ModifierList other = (ModifierList) obj;
        if (this.modifiers == null) {
            if (other.modifiers != null) {
                return false;
            }
        } else if (!this.modifiers.equals(other.modifiers)) {
            return false;
        }
        return true;
    }

    /**
     * @return the number of flow modifiers
     */
    public int size() {
        return this.modifiers.size();
    }

    /**
     * Save the modifiers to a file.
     *
     * @param saveFile	file to contain the modifiers
     *
     * @throws IOException
     */
    public void save(File saveFile) throws IOException {
        try (PrintWriter writer = new PrintWriter(saveFile)) {
            writer.println("command\tparms");
            StringBuilder line = new StringBuilder(80);
            for (Modifier mod : this.modifiers) {
                line.setLength(0);
                // Inactive modifiers are commented out.
                if (! mod.isActive())
                    line.append("#");
                // Format the modifier.
                line.append(mod.getCommand()).append('\t').append(mod.getParms());
                writer.println(line);
            }
        }
    }

    @Override
    public Iterator<Modifier> iterator() {
        return this.modifiers.iterator();
    }

}
