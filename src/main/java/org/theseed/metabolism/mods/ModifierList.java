/**
 *
 */
package org.theseed.metabolism.mods;

import java.io.File;
import java.io.IOException;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.metabolism.MetaModel;
import org.theseed.metabolism.Reaction.ActiveDirections;

import com.github.cliftonlabs.json_simple.JsonArray;
import com.github.cliftonlabs.json_simple.JsonKey;
import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This object represents a list of modifiers.  The modifiers are applied to a metabolic model to
 * to suppress directional flow of certain reactions and to otherwise modify construction of pathways.
 * The modifier list can be stored in JSON form or as a tab-delimited flat file with headers. In the latter
 * case, the command is in the first column and the parameter string in the second.
 *
 * @author Bruce Parrello
 *
 */
public class ModifierList {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ModifierList.class);
    /** list of (command, flow modifier) pairs */
    private List<Map.Entry<String, Modifier>> modifiers;

    /**
     * This enum defines the JSON keys we use.
     */
    private static enum FlowKeys implements JsonKey {
        COMMAND("SUPPRESS"), PARMS("");

        private final Object m_value;

        private FlowKeys(final Object value) {
            this.m_value = value;
        }

        /** This is the string used as a key in the incoming JsonObject map.
         */
        @Override
        public String getKey() {
            return this.name().toLowerCase();
        }

        /** This is the default value used when the key is not found.
         */
        @Override
        public Object getValue() {
            return this.m_value;
        }

    }

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
     * Create an empty modifier list.
     */
    public ModifierList() {
        this.modifiers = Collections.emptyList();
    }

    /**
     * Create a modifier list from a JSON object.
     *
     * @param json		JSON array containing the flow modifier descriptors
     *
     * @throws IOException
     */
    public ModifierList(JsonArray json) throws IOException {
        this.modifiers = new ArrayList<>(json.size());
        for (Object modObj : json) {
            JsonObject mod = (JsonObject) modObj;
            String command = mod.getStringOrDefault(FlowKeys.COMMAND);
            String parms = mod.getStringOrDefault(FlowKeys.PARMS);
            this.addModifier(command, parms);
        }
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
            this.addModifier(command, parms);
        }
    }

    /**
     * Add a flow modifier build from a command and a parameter string.
     *
     * @param command	command (all upper case)
     * @param parms		parameter string
     *
     * @throws IOException
     */
    private void addModifier(String command, String parms) throws IOException {
        Command commandCode;
        // An invalid command code is rethrown as an IO exception.
        try {
            commandCode = Command.valueOf(command);
        } catch (IllegalArgumentException e) {
            throw new IOException("Invalid flow modifier command code \"" + command + "\".");
        }
        // Build the modifier from the command and parameter string.
        Modifier modifier = commandCode.create(parms);
        this.modifiers.add(new AbstractMap.SimpleEntry<String, Modifier>(command, modifier));
    }

    /**
     * This creates the JSON representation of the modifier list.  Each modifier is represented by an
     * object containing the command name as "command" and the modifier parms as "parms".
     *
     * @return a JSON array for this modifier list
     */
    public JsonArray toJson() {
        JsonArray retVal = new JsonArray();
        for (Map.Entry<String, Modifier> modEntry : this.modifiers) {
            JsonObject modJson = new JsonObject().putChain("command", modEntry.getKey())
                    .putChain("parms", modEntry.getValue().getParms());
            retVal.add(modJson);
        }
        return retVal;
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
        for (Map.Entry<String, Modifier> modEntry : this.modifiers) {
            Modifier mod = modEntry.getValue();
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

}
