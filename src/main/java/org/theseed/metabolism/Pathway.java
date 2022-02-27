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
    public static class Element implements Comparable<Element> {

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
        // Only bother if this pathway has elements.
        Pathway retVal;
        if (this.elements.size() <= 0) {
            retVal = new Pathway(this.goal);
            retVal.goal = this.input;
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
            // The branch count is the number of successor reactions to this step's
            // output that do NOT lead to next step's output.
            var reactions = model.getSuccessors(intermediate);
            for (Reaction reaction : reactions) {
                boolean isBranch = reaction.getOutputs(intermediate).stream()
                        .allMatch(x -> ! x.getMetabolite().equals(next));
                if (isBranch) {
                    var branchSet = retVal.computeIfAbsent(intermediate, x -> new TreeSet<Reaction>());
                    branchSet.add(reaction);
                }
            }
        }
        // Now, finally, we have the branches for the terminus.
        String terminus = this.getLast().output;
        var reactions = model.getSuccessors(terminus);
        Set<Reaction> branchSet = retVal.computeIfAbsent(terminus, x -> new TreeSet<Reaction>());
        branchSet.addAll(reactions);
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
        String terminus = this.getLast().getOutput();
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
            retVal.append(this.getLast().getOutput());
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
     * @param filters	filters to apply
     */
    public Pathway extend(MetaModel model, String bigg2, PathwayFilter... filters) {
        return model.extendPathway(this, bigg2, filters);
    }

    /**
     * @return a pathway formed by looping this one back to a starting compound
     *
     * @param model		parent model of the pathway
     * @param filters	filters to apply
     */
    public Pathway loop(MetaModel model, PathwayFilter... filters) {
        return model.loopPathway(this, filters);
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


}
