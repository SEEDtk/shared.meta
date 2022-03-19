/**
 *
 */
package org.theseed.metabolism;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.apache.commons.lang3.StringUtils;
import org.theseed.counters.CountMap;
import org.theseed.excel.CustomWorkbook;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.utils.ParseFailureException;

import com.github.cliftonlabs.json_simple.JsonArray;
import com.github.cliftonlabs.json_simple.JsonException;
import com.github.cliftonlabs.json_simple.JsonKey;
import com.github.cliftonlabs.json_simple.JsonObject;
import com.github.cliftonlabs.json_simple.Jsoner;

/**
 * A pathway is an ordered set of reactions from one gene to another.  The rule is
 * that each reaction must connect to the next through a metabolite, and no reaction
 * can occur twice.  For each reaction we specify the direction and the output
 * metabolite.
 *
 * @author Bruce Parrello
 *
 */
public class Pathway implements Iterable<Pathway.Element>, Comparable<Pathway> {

    // FIELDS
    /** ordered list of pathway elements */
    private List<Element> elements;
    /** pathway goal compound */
    private String goal;
    /** pathway start compound */
    private String input;
    /** empty element list for Json */
    private static final JsonArray EMPTY_ARRAY = new JsonArray();
    /** filename extension for pathway files */
    public static final String FILE_EXT = ".path.json";

    /**
     * This enum defines the JSON keys we use.
     */
    private static enum PathwayKeys implements JsonKey {
        REVERSED(false), OUTPUT(null), INPUT(null), ELEMENTS(EMPTY_ARRAY);

        private final Object m_value;

        private PathwayKeys(final Object value) {
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
     * This represents a single reaction element of the pathway.  Elements are sorted by reaction ID, then
     * output metabolite, and finally with unreversed before reversed.
     */
    public class Element implements Comparable<Element>, IReactionSource {

        /** TRUE if the reaction is reversed */
        private boolean reversed;
        /** BiGG ID of the output metabolite */
        private String output;
        /** reaction at this point in the pathway */
        private Reaction reaction;
        /** sequence number of reaction */
        private int seqNum;

        /**
         * Create a new element for a pathway.
         *
         * @param reaction	reaction to use
         * @param node		desired output node
         * @param num 		sequence number of reaction
         */
        private Element(Reaction reaction, Reaction.Stoich node, int num) {
            this.reaction = reaction;
            this.output = node.getMetabolite();
            this.reversed = ! node.isProduct();
            this.seqNum = num;
        }

        /**
         * Create a new element for a pathway.
         *
         * @param reaction	reaction to use
         * @param compound	desired output compound
         * @param reverse	TRUE if the reaction should be reversed
         * @param num 		sequence number of reaction
         */
        private Element(Reaction reaction, String compound, boolean reverse, int num) {
            this.reaction = reaction;
            this.output = compound;
            this.reversed = reverse;
            this.seqNum = num;
        }

        /**
         * Create a pathway element from a JSON representation.
         *
         * @param json		JSON representation of the element
         * @param model		parent model to contain the pathway
         *
         * @throws ParseFailureException
         */
        private Element(JsonObject json, MetaModel model) throws ParseFailureException {
            JsonObject reactionO = (JsonObject) json.get("reaction");
            String reactionId = reactionO.getStringOrDefault(Reaction.ReactionKeys.BIGG_ID);
            this.reaction = model.getReaction(reactionId);
            if (this.reaction == null)
                throw new ParseFailureException("Pathway reaction " + reactionId
                        + " is not in model " + model.getMapName() + ".");
            this.output = json.getStringOrDefault(PathwayKeys.OUTPUT);
            this.reversed = json.getBooleanOrDefault(PathwayKeys.REVERSED);
        }

        /**
         * Create a new element by reversing an old one.
         *
         * @param old		element to reverse
         * @param output	desired output compound
         *
         * @throws IllegalArgumentException 	if the reaction is not reversible
         */
        protected Element(Element old, String output) {
            this.reaction = old.reaction;
            if (! this.reaction.isReversible())
                throw new IllegalArgumentException("Attempt to reverse irreversible reaction " + this.reaction.getBiggId()
                        + ".");
            this.reversed = ! old.reversed;
            this.output = output;
        }

        /**
         * @return TRUE if the reaction is reversed
         */
        @Override
        public boolean isReversed() {
            return this.reversed;
        }

        /**
         * @return the BiGG ID of the output metabolite
         */
        public String getOutput() {
            return this.output;
        }

        /**
         * @return the reaction
         */
        @Override
        public Reaction getReaction() {
            return this.reaction;
        }

        /**
         * @return the set of input compounds for this element
         */
        public Set<String> getInputs() {
            Collection<Reaction.Stoich> inputNodes = this.reaction.getOutputs(this.output);
            Set<String> retVal = inputNodes.stream().map(x -> x.getMetabolite()).collect(Collectors.toSet());
            return retVal;
        }

        /**
         * @return the previous element's output, or NULL if this is the first element
         */
        public String getMainInput() {
            int idx = Pathway.this.elements.indexOf(this);
            String retVal = null;
            if (idx > 0)
                retVal = Pathway.this.elements.get(idx - 1).output;
            return retVal;
        }

        @Override
        public String toString() {
            return "-(" + this.reaction.getBiggId() + ")" + this.output;
        }

        @Override
        public int compareTo(Element o) {
            int retVal = this.reaction.getBiggId().compareTo(o.reaction.getBiggId());
            if (retVal == 0) {
                retVal = this.output.compareTo(o.output);
                if (retVal == 0)
                    retVal = Boolean.compare(this.reversed, o.reversed);
            }
            return retVal;
        }

        /**
         * @return a JSON representation of this pathway element
         */
        public JsonObject toJson() {
            JsonObject retVal = new JsonObject();
            retVal.put("reversed", this.reversed);
            retVal.put("output", this.output);
            retVal.put("reaction", this.reaction.toJson());
            return retVal;
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((this.output == null) ? 0 : this.output.hashCode());
            result = prime * result + ((this.reaction == null) ? 0 : this.reaction.hashCode());
            result = prime * result + (this.reversed ? 1231 : 1237);
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof Element)) {
                return false;
            }
            Element other = (Element) obj;
            if (this.output == null) {
                if (other.output != null) {
                    return false;
                }
            } else if (!this.output.equals(other.output)) {
                return false;
            }
            if (this.reaction == null) {
                if (other.reaction != null) {
                    return false;
                }
            } else if (!this.reaction.equals(other.reaction)) {
                return false;
            }
            if (this.reversed != other.reversed) {
                return false;
            }
            return true;
        }

        /**
         * @return the sequence number
         */
        public int getSeqNum() {
            return this.seqNum;
        }

        /**
         * Specify the sequence number.
         * @param seqNum 	the sequence number to set
         */
        protected void setSeqNum(int seqNum) {
            this.seqNum = seqNum;
        }

        /**
         * @return the set of output compounds for this element
         */
        public Set<String> getOutputs() {
            var retVal = reaction.getMetabolites().stream().filter(x -> (x.isProduct() != this.reversed))
                    .map(x -> x.getMetabolite()).collect(Collectors.toSet());
            return retVal;
        }

        @Override
        public Set<String> getSpecial() {
            // For a pathway element, the input and output are special.
            var retVal = new TreeSet<String>();
            String input = this.getMainInput();
            if (input != null)
                retVal.add(input);
            retVal.add(this.output);
            return retVal;
        }

    }

    /**
     * This class implements a filename filter for pathway files.
     */
    public static class FileFilter implements java.io.FileFilter {

        @Override
        public boolean accept(File pathname) {
            return pathname.getName().endsWith(FILE_EXT);
        }

    }

    /**
     * This class represents triggering information for the Excel output.
     */
    protected static class Trigger implements Comparable<Trigger> {

        /** triggering gene ID */
        private String bigg;
        /** target of the trigger */
        private String target;
        /** TRUE for a branching trigger, else FALSE */
        private boolean branch;
        /** weight of trigger */
        private double weight;

        /**
         * Construct a trigger descriptor.
         *
         * @param bigg		BiGG ID of the triggering gene
         * @param target	name of the trigger target (either a compound or a reaction)
         * @param branch	TRUE for a branching trigger, FALSE for a normal one
         * @param weight	weight of the trigger; a higher weight indicates more importance
         */
        protected Trigger(String bigg, String target, boolean branch, double weight) {
            this.bigg = bigg;
            this.target = target;
            this.branch = branch;
            this.weight = weight;
        }

        @Override
        public int compareTo(Trigger o) {
            // Float more important triggers to the top.
            int retVal = Double.compare(o.weight, this.weight);
            if (retVal == 0) {
                // Branches are more important than mainline triggers. (Note that FALSE preceeds TRUE.)
                retVal = Boolean.compare(o.branch, this.branch);
                if (retVal == 0) {
                    // Sort by target.
                    retVal = this.target.compareTo(o.target);
                    if (retVal == 0) {
                        // Finally, sort by feature ID.
                        retVal = this.bigg.compareTo(o.bigg);
                    }
                }
            }
            return retVal;
        }
    }

    /**
     * Construct an empty pathway with a specified input.
     *
     * @param input		BiGG ID of the input compound
     */
    public Pathway(String input) {
        this.elements = new ArrayList<Element>();
        this.input = input;
    }

    /**
     * Construct a pathway from a single reaction element.
     *
     * @param input		input compound for the pathway
     * @param reaction	first reaction in pathway
     * @param stoich	stoichiometric element for computing output
     */
    public Pathway(String input, Reaction reaction, Reaction.Stoich node) {
        this.input = input;
        this.elements = new ArrayList<Element>();
        add(reaction, node);
    }

    /**
     * Construct a pathway from a JSON representation.
     *
     * @param json		JSON representation of the pathway
     * @param model		parent model of the pathway
     *
     * @throws ParseFailureException
     */
    public Pathway(JsonObject json, MetaModel model) throws ParseFailureException {
        this.fromJson(json, model);
    }

    /**
     * Construct a pathway from a json string.
     *
     * @param jsonString	JSON string describing the pathway
     * @param model			parent model of the pathway
     *
     * @throws ParseFailureException
     */
    public Pathway(String jsonString, MetaModel model) throws ParseFailureException {
        JsonObject json = Jsoner.deserialize(jsonString, (JsonObject) null);
        if (json == null)
            throw new ParseFailureException("Cannot deserialize JSON string \"" + StringUtils.abbreviate(jsonString, 30)
                    + "\".");
        this.fromJson(json, model);
    }

    /**
     * Construct a pathway from a json file.
     *
     * @param inFile	file containing the JSON representation of the pathway
     * @param model		parent model of the pathway
     *
     * @throws ParseFailureException
     * @throws JsonException
     * @throws IOException
     */
    public Pathway(File inFile, MetaModel model) throws ParseFailureException, JsonException, IOException {
        try (Reader reader = new FileReader(inFile)) {
            JsonObject json = (JsonObject) Jsoner.deserialize(reader);
            this.fromJson(json, model);
        }
    }

    /**
     * Construct a pathway from a single reaction element and specify a goal
     *
     * @param input		input compound
     * @param reaction	first reaction in pathway
     * @param stoich	stoichiometric element for computing output
     * @param goal		target output compound
     */
    public Pathway(String input, Reaction reaction, Reaction.Stoich node, String goal) {
        this.input = input;
        this.elements = new ArrayList<Element>();
        add(reaction, node);
        this.goal = goal;
    }

    /**
     * Initialize this pathway from a JSON representation.
     *
     * @param json		JSON representation of the pathway
     * @param model		parent model of the pathway
     *
     * @throws ParseFailureException
     */
    private void fromJson(JsonObject json, MetaModel model) throws ParseFailureException {
        this.input = json.getStringOrDefault(PathwayKeys.INPUT);
        JsonArray pathJson = json.getCollectionOrDefault(PathwayKeys.ELEMENTS);
        this.elements = new ArrayList<Element>(pathJson.size());
        for (Object obj : pathJson) {
            JsonObject elementO = (JsonObject) obj;
            Element element = new Element(elementO, model);
            this.elements.add(element);
            // We must set the sequence number manually-- the numbers are 1-based.
            element.setSeqNum(this.elements.size());
        }
    }

    /**
     * Add a new reaction to this pathway.
     *
     * @param reaction	reaction to add
     * @param node		stoichiometric element indicating the output
     */
    public Pathway add(Reaction reaction, Reaction.Stoich node) {
        Element element = new Element(reaction, node, this.elements.size() + 1);
        this.elements.add(element);
        return this;
    }

    /**
     * @return TRUE if the specified reaction is already in this pathway
     *
     * @param reaction	reaction to check
     */
    public boolean contains(Reaction reaction) {
        boolean retVal = this.elements.stream().anyMatch(x -> x.reaction.equals(reaction));
        return retVal;
    }

    /**
     * @return a copy of this pathway.
     */
    public Pathway clone() {
        Pathway retVal = new Pathway(this.input);
        this.elements.stream().forEach(x -> retVal.elements.add(x));
        retVal.goal = this.goal;
        return retVal;
    }

    /**
     * @return TRUE if this pathway is reversible
     */
    public boolean isReversible() {
        boolean retVal = this.elements.stream().allMatch(x -> x.getReaction().isReversible());
        return retVal;
    }

    /**
     * Attempt to reverse this pathway.  If the pathway is not reversible, an IllegalArgumentException
     * will be thrown.
     *
     * @return the reverse of this pathway
     *
     */
    public Pathway reverse() {
        // Only bother if this pathway has elements.  A null path reverses to itself.
        Pathway retVal;
        if (this.elements.size() <= 0) {
            retVal = this.clone();
        } else {
            // Create the new pathway.
            Element lastElement = this.getLast();
            retVal = new Pathway(lastElement.output);
            // Compute the new outputs for each reaction.
            String[] outputs = new String[this.size()];
            outputs[0] = this.input;
            final int n = this.size() - 1;
            for (int i = 1; i <= n; i++)
                outputs[i] = this.getElement(i-1).output;
            // Now assemble the reactions in reverse order.
            for (int i = n; i >= 0; i--) {
                Element reversed = new Element(this.getElement(i), outputs[i]);
                retVal.elements.add(reversed);
            }
        }
        return retVal;
    }

    /**
     * @return the last pathway element
     */
    public Element getLast() {
        final int n = this.elements.size();
        Element retVal = null;
        if (n > 0)
            retVal = this.elements.get(n - 1);
        return retVal;
    }

    /**
     * @return the first pathway element
     */
    public Element getFirst() {
        Element retVal = null;
        if (! this.elements.isEmpty())
            retVal = this.elements.get(0);
        return retVal;
    }

    @Override
    public Iterator<Element> iterator() {
        return this.elements.iterator();
    }

    /**
     * @return an iterator through the tail of the pathway
     *
     * @param n		number of reactions to iterate through
     */
    public Iterator<Element> tailIterator(int n) {
        int i2 = this.elements.size();
        int i1 = i2 - n;
        if (i1 < 0) i1 = 0;
        return this.elements.subList(i1, i2).iterator();
    }

    /**
     * @return the number of segments in this path
     */
    public int size() {
        return this.elements.size();
    }

    /**
     * @return the specified pathway element
     *
     * @param i		index of the element to return
     */
    public Element getElement(int i) {
        return this.elements.get(i);
    }

    /**
     * @return TRUE if this path includes all the reactions in the specified set, else FALSE
     *
     * @param includes	set of BiGG IDs for required reactions
     */
    public boolean includesAll(Set<String> includes) {
        Iterator<String> iter = includes.iterator();
        boolean retVal = true;
        while (iter.hasNext() && retVal) {
            String rID = iter.next();
            final int n = this.elements.size();
            boolean found = false;
            for (int i = 0; i < n && ! found; i++) {
                Reaction r = this.elements.get(i).reaction;
                if (r.getBiggId().equals(rID))
                    found = true;
            }
            retVal = found;
        }
        return retVal;
    }

    /**
     * Find the branches off each node in the pathway.
     *
     * @param model		the model in which the pathway occurs
     *
     * @return a map from each output compound to the set of branching reactions
     */
    public Map<String, Set<Reaction>> getBranches(MetaModel model) {
        var retVal = new HashMap<String, Set<Reaction>>(this.size() * 4 / 3);
        final int n = this.size() - 1;
        for (int i = 0; i < n; i++) {
            // Get the output compounds for this step and the next step.
            String intermediate = this.elements.get(i).output;
            String next = this.elements.get(i+1).output;
            // The branches are the successor reactions to this step's
            // output that do NOT lead to next step's output and are not
            // reverses of a reaction in the path.
            var reactions = model.getSuccessors(intermediate);
            for (Reaction reaction : reactions) {
                boolean isBranch = (reaction.getOutputs(intermediate).stream()
                        .allMatch(x -> ! x.getMetabolite().equals(next)) &&
                        ! this.contains(reaction));
                if (isBranch) {
                    var branchSet = retVal.computeIfAbsent(intermediate, x -> new TreeSet<Reaction>());
                    branchSet.add(reaction);
                }
            }
        }
        // Now, finally, we have the branches for the terminus.  We skip it if it is external.  Branching
        // doesn't make sense for an external compound.
        String terminus = this.getOutput();
        if (! terminus.endsWith("_e")) {
            var reactions = model.getSuccessors(terminus);
            Set<Reaction> branchSet = retVal.computeIfAbsent(terminus, x -> new TreeSet<Reaction>());
            branchSet.addAll(reactions);
        }
        return retVal;
    }

    /**
     * @return a stream of the pathway elements
     */
    public Stream<Element> stream() {
        return this.elements.stream();
    }

    /**
     * Specify the goal compound for this pathway.
     *
     * @param goal 	the desired output compound
     */
    public void setGoal(String goal) {
        this.goal = goal;
    }

    /**
     * @return TRUE if this pathway has reached its goal
     */
    public boolean isComplete() {
        String terminus = this.getOutput();
        return terminus.contentEquals(goal);
    }

    /**
     * @return the set of intermediate compounds for this pathway
     */
    public Set<String> getIntermediates() {
        var retVal = IntStream.range(0, this.size() - 1).mapToObj(i -> this.getElement(i).getOutput())
                .collect(Collectors.toSet());
        return retVal;
    }

    /**
     * @return the goal compound
     */
    public String getGoal() {
        return this.goal;
    }

    @Override
    public String toString() {
        StringBuffer retVal = new StringBuffer(80);
        if (this.input == null)
            retVal.append("X");
        else
            retVal.append(this.input);
        retVal.append(" ==> ");
        if (this.goal != null)
            retVal.append(this.goal);
        else if (this.elements.size() == 0)
            retVal.append("X");
        else
            retVal.append(this.getOutput());
        return retVal.toString();
    }

    @Override
    public int compareTo(Pathway o) {
        // We sort the shortest pathway first.
        int retVal = this.size() - o.size();
        if (retVal == 0) {
            // Here the pathways are the same number of elements.  Compare each element until we find a difference.
            final int n = this.size();
            for (int i = 0; retVal == 0 && i < n; i++)
                retVal = this.elements.get(i).compareTo(o.elements.get(i));
        }
        return retVal;
    }

    /**
     * @return a JSON representation of this pathway
     */
    public JsonObject toJson() {
        JsonArray elementJson = new JsonArray();
        for (Element element : this)
            elementJson.add(element.toJson());
        JsonObject retVal = new JsonObject().putChain("input", this.input).putChain("elements", elementJson);
        return retVal;
    }

    /**
     * @return a string representation of this pathway
     */
    public String toJsonString() {
        return Jsoner.serialize(this.toJson());
    }

    /**
     * @return a pathway formed by extending this one to a new goal
     *
     * @param model		parent model of the pathway
     * @param bigg2		BiGG ID of the output compound
     */
    public Pathway extend(MetaModel model, String bigg2) {
        return model.extendPathway(this, bigg2);
    }

    /**
     * @return a pathway formed by looping this one back to a starting compound
     *
     * @param model		parent model of the pathway
     */
    public Pathway loop(MetaModel model) {
        return model.loopPathway(this);
    }

    /**
     * Save this pathway to a file.
     *
     * @param outFile	output file for saving the pathway
     *
     * @throws IOException
     */
    public void save(File outFile) throws IOException {
        String jsonString = this.toJsonString();
        try (PrintWriter writer = new PrintWriter(outFile)) {
            writer.println(Jsoner.prettyPrint(jsonString));
        }
    }

    /**
     * @return the list of outputs for this pathway (in order)
     */
    public List<String> getOutputs() {
        return this.elements.stream().map(x -> x.getOutput()).collect(Collectors.toList());
    }


    /**
     * Compute the set of inputs for this pathway.  An input is a compound that does not appear on
     * the output side of any reaction and that is not uncommon.
     *
     * @param model			model containing this pathway
     * @param includeAll	TRUE to include common inputs
     *
     * @return a map containing the number of each input metabolite required
     */
    public CountMap<String> getInputs(MetaModel model, boolean includeAll) {
        // Create the exclude set.
        Set<String> commons;
        if (includeAll)
            commons = Collections.emptySet();
        else
            commons = model.getCommons();
        // Build the count map.
        CountMap<String> retVal = new CountMap<String>();
        Set<String> outputs = new TreeSet<String>();
        for (Element element : this.elements) {
            Reaction reaction = element.getReaction();
            boolean mode = element.isReversed();
            for (Reaction.Stoich stoich : reaction.getMetabolites()) {
                String metabolite = stoich.getMetabolite();
                // Note we only process the uncommon compounds.
                if (! commons.contains(metabolite)) {
                    if (mode == stoich.isProduct())
                        retVal.count(stoich.getMetabolite(), stoich.getCoeff());
                    else
                        outputs.add(stoich.getMetabolite());
                }
            }
        }
        // Now we must erase the metabolites that are also outputs.
        outputs.stream().forEach(x -> retVal.remove(x));
        // Return the remaining counts.
        return retVal;
    }

    /**
     * @return the input compound for this path
     */
    public String getInput() {
        return this.input;
    }

    /**
     * @return the current output metabolite for this path.
     */
    public String getOutput() {
        String retVal;
        if (this.elements.size() == 0)
            retVal = this.getInput();
        else
            retVal = this.getLast().output;
        return retVal;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.elements == null) ? 0 : this.elements.hashCode());
        result = prime * result + ((this.input == null) ? 0 : this.input.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof Pathway)) {
            return false;
        }
        Pathway other = (Pathway) obj;
        if (this.elements == null) {
            if (other.elements != null) {
                return false;
            }
        } else if (!this.elements.equals(other.elements)) {
            return false;
        }
        if (this.input == null) {
            if (other.input != null) {
                return false;
            }
        } else if (!this.input.equals(other.input)) {
            return false;
        }
        return true;
    }

    /**
     * Write this path to an excel file.
     *
     * @param excelFile		output excel file name
     * @param model			underlying metabolic model
     */
    public void saveToExcel(File excelFile, MetaModel model) {
        // Get the compound ratings.
        var weightMap = CompoundRating.getRatingMap(this, model);
        // This will collect the triggering genes.
        var genes = new TreeSet<Trigger>();
        // This will collect the input metabolites.
        CountMap<String> inputCounts = new CountMap<String>();
        // Get the common compounds.
        Set<String> commons = model.getCommons();
        // Get the branching reactions.
        var branches = this.getBranches(model);
        // This will hold the triggering protein ratings.
        Map<String, ProteinRating> inserts = new HashMap<String, ProteinRating>(this.size() * 5);
        // Each pathway element transmits a direct-line input to an output.
        // The inputs we want to count are the ones not in the direct line.
        // The first direct-line input is the main input.
        String oldInput = this.input;
        // Start the pathway report.
        try (CustomWorkbook workbook = CustomWorkbook.create(excelFile)) {
            workbook.addSheet("Pathway", true);
            workbook.setHeaders(Arrays.asList("reaction", "reaction_name", "triggering_rule", "output",
                    "formula"));
            for (Pathway.Element element : this) {
                Reaction reaction = element.getReaction();
                String intermediate = element.getOutput();
                String rule = reaction.getReactionRule();
                // Find out if the reaction is reversed.
                boolean reversed = element.isReversed();
                // Compute the weight.
                double weight = reaction.getWeight(weightMap, ! reversed);
                // Apply the weight to the reaction rule to get the protein weights.
                String translatedRule = Reaction.getTranslatedRule(rule, model);
                applyWeights(inserts, true, weight, translatedRule, reaction, reversed, intermediate);
                // Add the reaction row.
                workbook.addRow();
                workbook.storeCell(reaction.getBiggId());
                workbook.storeCell(reaction.getName());
                workbook.storeCell(translatedRule);
                workbook.storeCell(element.getOutput());
                workbook.storeCell(reaction.getLongFormula(element.isReversed(), model));
                // Add the triggering genes to the gene set.
                reaction.getTriggers().stream()
                        .forEach(x -> genes.add(new Trigger(x, reaction.getBiggId(), false, weight)));
                // Get the reaction inputs.
                var inputs = reaction.getOutputs(intermediate);
                for (Reaction.Stoich input : inputs) {
                    // Note we don't count the direct-line input, just the ancillaries.
                    if (! input.getMetabolite().equals(oldInput))
                        inputCounts.count(input.getMetabolite(), Math.abs(input.getCoeff()));
                }
                // Remember our output as the direct-line input for the next reaction.
                oldInput = intermediate;
            }
            workbook.autoSizeColumns();
            workbook.addSheet("Inputs", true);
            // Now list the inputs.
            workbook.setHeaders(Arrays.asList("metabolite", "needed", "name"));
            for (CountMap<String>.Count counter : inputCounts.sortedCounts()) {
                workbook.addRow();
                String compound = counter.getKey();
                workbook.storeCell(compound);
                workbook.storeCell(counter.getCount());
                workbook.storeCell(model.getCompoundName(compound));
            }
            workbook.autoSizeColumns();
            // Now we are working with genes, so we need the base genome.
            Genome baseGenome = model.getBaseGenome();
            // The next sheet is the branch reactions.  For each one we we want to show the input metabolite
            // and the details of the reaction itself.  We also track the triggering genes for the branches.
            Map<String, ProteinRating> deletes = new HashMap<String, ProteinRating>(branches.size() * 5);
            workbook.addSheet("Branches", true);
            workbook.setHeaders(Arrays.asList("input", "reaction", "reaction_name", "triggering_rule",
                    "formula"));
            for (Map.Entry<String, Set<Reaction>> branchList : branches.entrySet()) {
                String input = branchList.getKey();
                for (Reaction reaction : branchList.getValue()) {
                    // Check to see if there is a rule.  If there is no rule, we skip the branch.
                    String rule = reaction.getReactionRule();
                    if (! StringUtils.isBlank(rule)) {
                        // We need to see if the input requires reversing the reaction.
                        boolean reverse = reaction.isProduct(input);
                        // From that we compute the weight.
                        double weight = reaction.getWeight(weightMap, reverse);
                        String translatedRule = Reaction.getTranslatedRule(rule, model);
                        // If the compound is uncommon, get the protein weights for the reaction.
                        if (! commons.contains(input))
                            applyWeights(deletes, false, weight, translatedRule, reaction, reverse, input);
                        // Now create the output row.
                        workbook.addRow();
                        workbook.storeCell(input);
                        workbook.storeCell(reaction.getBiggId());
                        workbook.storeCell(reaction.getName());
                        workbook.storeCell(translatedRule);
                        workbook.storeCell(reaction.getLongFormula(reverse, model));
                        // Add the triggering genes to the gene set.
                        reaction.getTriggers().forEach(x -> genes.add(new Trigger(x, input, true, weight)));
                    }
                }
            }
            workbook.autoSizeColumns();
            // Write the triggering analysis. There are good ones that trigger the path, and bad
            // ones that bleed off metabolites into other reactions.
            workbook.addSheet("Triggers", true);
            workbook.setHeaders(Arrays.asList("weight", "gene", "target", "fid", "aliases", "type", "function"));
            this.writeGenes(workbook, genes, baseGenome);
            workbook.autoSizeColumns();
            // Finally, we write the protein analysis.
            List<ProteinRating> protList = new ArrayList<ProteinRating>(inserts.size() + deletes.size());
            protList.addAll(inserts.values());
            protList.addAll(deletes.values());
            Collections.sort(protList);
            workbook.addSheet("Genes", true);
            workbook.setHeaders(Arrays.asList("gene", "weight", "compound", "reaction"));
            for (ProteinRating rating : protList) {
                workbook.addRow();
                workbook.storeCell(rating.getProteinSpec());
                workbook.storeCell(rating.getWeight());
                String compound = model.getCompoundName(rating.getTarget());
                workbook.storeCell(compound);
                String formula = rating.getReaction().getLongFormula(rating.isReversed(), model);
                workbook.storeCell(formula);
            }
            workbook.autoSizeColumns();
        }

    }

    /**
     * Apply a reaction's weight to the proteins in its reaction rule.
     *
     * @param ratingMap		map of protein IDs to ratings
     * @param type			TRUE for insert, FALSE for delete
     * @param weight		weight of the reaction
     * @param rule			reaction rule string
     * @param reaction		target reaction
     * @param reverse 		TRUE if the reaction is reversed
     * @param compound		ID of the relevant compound
     */
    public static void applyWeights(Map<String, ProteinRating> ratingMap, boolean type, double weight, String rule,
            Reaction reaction, boolean reverse, String compound) {
        ReactionRule parsed = ReactionRule.parse(rule);
        var weightMap = (type ? parsed.getTriggerWeights() : parsed.getBranchWeights());
        for (Map.Entry<String, Double> weightEntry : weightMap.entrySet()) {
            String prot = weightEntry.getKey();
            ProteinRating rating = ratingMap.computeIfAbsent(prot, x -> new ProteinRating(x, type));
            rating.add(weightEntry.getValue() * weight, reaction, reverse, compound);
        }
    }

    /**
     * This method will write a set of genes to the triggering worksheet.
     *
     * @param workbook		output workbook
     * @param genes			sorted set of genes to write
     * @param baseGenome	base genome for the current model
     */
    private void writeGenes(CustomWorkbook workbook, Set<Trigger> genes, Genome baseGenome) {
        var aliasMap = baseGenome.getAliasMap();
        for (Trigger geneEntry : genes) {
            double weight = geneEntry.weight;
            String gene = geneEntry.bigg;
            var fids = aliasMap.get(gene);
            if (fids == null) {
                workbook.addRow();
                workbook.storeCell(weight);
                workbook.storeCell(gene);
                workbook.storeBlankCell();
                workbook.storeBlankCell();
                workbook.storeBlankCell();
                workbook.storeBlankCell();
                workbook.storeBlankCell();
            } else {
                for (String fid : fids) {
                    workbook.addRow();
                    Feature feat = baseGenome.getFeature(fid);
                    var aliases = feat.getAliases();
                    String aliasList = StringUtils.join(aliases, ", ");
                    workbook.storeCell(weight);
                    workbook.storeCell(gene);
                    workbook.storeCell(geneEntry.target);
                    workbook.storeCell(fid);
                    workbook.storeCell(aliasList);
                    workbook.storeCell(geneEntry.branch ? "branch" : "trigger");
                    workbook.storeCell(feat.getPegFunction());
                }
            }
        }
    }

    /**
     * @return the element preceding the specified one, or NULL if it is the first
     *
     * @param element	element of interest
     */
    public Element getPrevious(Element element) {
        Element retVal = null;
        if (element != this.elements.get(0)) {
            for (int i = 1; i < this.elements.size() && retVal == null; i++) {
                if (this.elements.get(i) == element)
                    retVal = this.elements.get(i-1);
            }
        }
        return retVal;
    }

    /**
     * Append the specified path to this one.  This will only work if the input of the
     * second path is the output of this path.
     *
     * @param path2		pathway to append
     */
    public void append(Pathway path2) {
        for (Pathway.Element element : path2) {
            Reaction react = element.getReaction();
            this.add(react, react.getStoich(element.getOutput()));
        }
    }

    /**
     * Add the specified reaction to this pathway.
     *
     * @param reaction		reaction to add
     * @param compoundId	proposed output compound
     * @param reverse		TRUE to reverse the reaction
     */
    public void add(Reaction reaction, String compoundId, boolean reverse) {
        int num = this.elements.size() + 1;
        Element element = new Element(reaction, compoundId, reverse, num);
        this.elements.add(element);
    }

}
