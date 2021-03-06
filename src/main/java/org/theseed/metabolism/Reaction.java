/**
 *
 */
package org.theseed.metabolism;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import com.github.cliftonlabs.json_simple.JsonArray;
import com.github.cliftonlabs.json_simple.JsonKey;
import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This object represents a reaction present in a metabolic model.  The reaction contains
 * segments that describe the nodes to which it connects, the list of genes that trigger it,
 * various bits of identifiying information, and the stoichiometry of the compounds involved.
 *
 * The reaction also contains an active-directions status.  Depending on the state of the cell,
 * a reaction may be suppressed completely, or have one or both directions of flow shut down.
 * The active-directions status is reflected in the "isInput" and "isOutput" methods.
 *
 * @author Bruce Parrello
 *
 */
public class Reaction implements Comparable<Reaction> {

    // FIELDS
    /** ID number of this reaction */
    private int id;
    /** name of the reaction */
    private String name;
    /** BiGG identifier */
    private String biggId;
    /** rule for triggering the reaction (text based on BiGG IDs) */
    private String reactionRule;
    /** reversibility flag */
    private boolean reversible;
    /** list of gene aliases */
    private Map<String, String> aliases;
    /** metabolite list */
    private List<Stoich> metabolites;
    /** connections list */
    private List<Segment> segments;
    /** display location of reaction label */
    private Coordinate labelLoc;
    /** active direction indicator */
    private ActiveDirections active;
    /** connectors for reaction rules */
    private static final Set<String> CONNECTORS = Set.of("and", "or", "not");
    /** delimiters for reaction rules */
    private static final Pattern DELIMITERS = Pattern.compile("[ ()]+");
    /** delimiter set for reaction rules */
    private static final String DELIM_STRING = " ()";

    /**
     * This enumeration determines the active directions of a reaction.
     */
    public static enum ActiveDirections {
        BOTH {

            @Override
            public boolean isInput(Stoich stoich) {
                return true;
            }

            @Override
            public boolean isOutput(Stoich stoich) {
                return true;
            }

            @Override
            public String getCommand() {
                return "REVERSIBLE";
            }

        }, FORWARD {

            @Override
            public boolean isInput(Stoich stoich) {
                return ! stoich.isProduct();
            }

            @Override
            public boolean isOutput(Stoich stoich) {
                return stoich.isProduct();
            }

            @Override
            public String getCommand() {
                return "ONEWAY";
            }

        }, REVERSE {

            @Override
            public boolean isInput(Stoich stoich) {
                return stoich.isProduct();
            }

            @Override
            public boolean isOutput(Stoich stoich) {
                return ! stoich.isProduct();
            }

            @Override
            public String getCommand() {
                return "INVERTED";
            }

        }, NEITHER {

            @Override
            public boolean isInput(Stoich stoich) {
                return false;
            }

            @Override
            public boolean isOutput(Stoich stoich) {
                return false;
            }

            @Override
            public String getCommand() {
                return "SUPPRESS";
            }

        };

        /**
         * @return TRUE if the specified stoichiometric element is an eligible input in this active direction
         *
         * @param stoich	stoichiometric element of interest
         */
        public abstract boolean isInput(Stoich stoich);

        /**
         * @return TRUE if the specified stoichiometric element is an eligible output in this active direction
         *
         * @param stoich	stoichiometric element of interest
         */
        public abstract boolean isOutput(Stoich stoich);

        /**
         * Store the default active-direction indicator for a reaction.
         *
         * @param react		reaction to set up
         */
        protected static void setDefault(Reaction react) {
            react.active = (react.isReversible() ? BOTH : FORWARD);
        }

        /**
         * @return the flow modifier command code for this direction
         */
        public abstract String getCommand();

    }

    protected static enum ReactionKeys implements JsonKey {
        NAME("<unknown>"), BIGG_ID(""), REVERSIBILITY(false), LABEL_X(0.0), LABEL_Y(0.0),
        GENE_REACTION_RULE(""), COEFFICIENT(1);

        private final Object m_value;

        private ReactionKeys(final Object value) {
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

    private static enum SegmentKeys implements JsonKey {
        FROM_NODE_ID(0), TO_NODE_ID(0), B1(null), B2(null), X(0.0), Y(0.0);

        private final Object m_value;

        private SegmentKeys(final Object value) {
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
     * This is a simple object to represent stoichiometry.  The sort order puts reactants
     * before products.
     */
    public static class Stoich implements Comparable<Stoich> {

        /** stoichiometric coefficient */
        private int coefficient;
        /** BiGG identifier of metabolite */
        private String biggId;

        /**
         * Construct a new stoichiometric representation.
         *
         * @param coeff		coefficient
         * @param bigg		BiGG ID for metabolite
         */
        public Stoich(int coeff, String bigg) {
            this.coefficient = coeff;
            this.biggId = bigg;
        }

        @Override
        public int compareTo(Stoich o) {
            int retVal = this.coefficient - o.coefficient;
            if (retVal == 0)
                retVal = this.biggId.compareTo(o.biggId);
            return retVal;
        }

        /**
         * @return the coefficient (always positive)
         */
        public int getCoeff() {
            return (this.coefficient < 0 ? -this.coefficient : this.coefficient);
        }

        /**
         * @return TRUE for a product, FALSE for a reactant
         */
        public boolean isProduct() {
            return (this.coefficient > 0);
        }

        @Override
        public String toString() {
            int coeff = this.getCoeff();
            String retVal;
            if (coeff == 1)
                retVal = this.biggId;
            else
                retVal = String.format("%d*%s", coeff, this.biggId);
            return retVal;
        }

        public String getMetabolite() {
            return this.biggId;
        }

    }

    /**
     * This class represents a connection made by a reaction.  It contains the IDs of the
     * source and destination nodes, and the display coordinates (if any).
     * @author Bruce Parrello
     *
     */
    public static class Segment {

        /** ID of the segment */
        private int id;
        /** ID of origin node */
        private int fromNode;
        /** ID of terminal node */
        private int toNode;
        /** first location */
        private Coordinate b1;
        /** second location */
        private Coordinate b2;

        /**
         * Construct a segment from a JSON object.
         *
         * @param segId			segment ID;
         * @param segObject		source JSON object for the segment
         */
        public Segment(int segId, JsonObject segObject) {
            this.id = segId;
            this.fromNode = segObject.getIntegerOrDefault(SegmentKeys.FROM_NODE_ID);
            this.toNode = segObject.getIntegerOrDefault(SegmentKeys.TO_NODE_ID);
            this.b1 = this.getCoordinates(segObject, SegmentKeys.B1);
            this.b2 = this.getCoordinates(segObject, SegmentKeys.B2);
        }

        /**
         * @return the coordinates specified with the indicated key
         *
         * @param segObject		segment JSON object
         * @param string		segment key to use
         */
        private Coordinate getCoordinates(JsonObject segObject, SegmentKeys string) {
            JsonObject coord = (JsonObject) segObject.getOrDefault(string.getKey(), string.getValue());
            Coordinate retVal = null;
            if (coord != null)
                retVal = new Coordinate(coord.getDoubleOrDefault(SegmentKeys.X),
                        coord.getDoubleOrDefault(SegmentKeys.Y));
            return retVal;
        }

        /**
         * @return the from-node ID number
         */
        public int getFromNode() {
            return this.fromNode;
        }

        /**
         * @return the to-node ID number
         */
        public int getToNode() {
            return this.toNode;
        }

        /**
         * @return the ID of this segment
         */
        public int getId() {
            return this.id;
        }

        /**
         * @return a JSON object for this segment
         */
        public JsonObject toJson() {
            JsonObject retVal = new JsonObject();
            retVal.put("from_node_id", this.fromNode);
            retVal.put("to_node_id", this.toNode);
            if (this.b1 == null) {
                retVal.put("b1", null);
                retVal.put("b2", null);
            } else {
                retVal.put("b1", this.b1.toJson());
                retVal.put("b2", this.b2.toJson());
            }
            return retVal;
        }

        /**
         * @return the start coordinate
         */
        public Coordinate getB1() {
            return this.b1;
        }

        /**
         * @return the end coordinate
         */
        public Coordinate getB2() {
            return this.b2;
        }

    }

    /**
     * Construct a reaction from a JSON object.
     *
     * @param id				ID of this reaction
     * @param reactionObject	JSON object containing the reaction
     */
    public Reaction(int id, JsonObject reactionObject) {
        this.id = id;
        this.name = reactionObject.getStringOrDefault(ReactionKeys.NAME);
        this.biggId = reactionObject.getStringOrDefault(ReactionKeys.BIGG_ID);
        this.reversible = reactionObject.getBooleanOrDefault(ReactionKeys.REVERSIBILITY);
        this.labelLoc = new Coordinate(reactionObject.getDoubleOrDefault(ReactionKeys.LABEL_X),
                reactionObject.getDoubleOrDefault(ReactionKeys.LABEL_Y));
        this.reactionRule = reactionObject.getStringOrDefault(ReactionKeys.GENE_REACTION_RULE);
        ActiveDirections.setDefault(this);
        // For the gene list, we extract all the aliases and store them in a map.
        JsonArray geneList = (JsonArray) reactionObject.get("genes");
        this.aliases = new TreeMap<String, String>();
        if (geneList != null) {
            for (int i = 0; i < geneList.size(); i++) {
                JsonObject gene = (JsonObject) geneList.get(i);
                String bigg_id = (String) gene.get("bigg_id");
                String name = (String) gene.get("name");
                if (! StringUtils.isBlank(bigg_id))
                    this.addAlias(bigg_id, name);
            }
        }
        // For the metabolites list, we convert each one to a stoichiometry.
        List<Stoich> reactionParts;
        JsonArray metaList = (JsonArray) reactionObject.get("metabolites");
        if (metaList == null)
            reactionParts = Collections.emptyList();
        else {
            reactionParts = new ArrayList<Stoich>();
            // Get all the pieces.
            for (int i = 0; i < metaList.size(); i++) {
                JsonObject meta = (JsonObject) metaList.get(i);
                Stoich stoich = new Stoich(meta.getIntegerOrDefault(ReactionKeys.COEFFICIENT),
                        meta.getStringOrDefault(ReactionKeys.BIGG_ID));
                reactionParts.add(stoich);
            }
            // Sort them into a list.
            Collections.sort(reactionParts);
        }
        this.metabolites = reactionParts;
        // Finally, we process the segments.
        JsonObject segList = (JsonObject) reactionObject.get("segments");
        if (segList == null)
            this.segments = Collections.emptyList();
        else {
            this.segments = new ArrayList<Segment>();
            for (Map.Entry<String, Object> segEntry : segList.entrySet()) {
                int segId = Integer.valueOf(segEntry.getKey());
                JsonObject segItem = (JsonObject) segEntry.getValue();
                Segment seg = new Segment(segId, segItem);
                this.segments.add(seg);
            }
        }
    }

    /**
     * Create a new reaction with the specified name.
     *
     * @param reactionId		ID number of the reaction
     * @param biggId			BiGG ID of the reaction
     * @param reactionName		name of the reaction
     */
    public Reaction(int reactionId, String biggId, String reactionName) {
        this.id = reactionId;
        this.biggId = biggId;
        this.name = reactionName;
        this.aliases = new TreeMap<String, String>();
        this.labelLoc = null;
        this.metabolites = new ArrayList<Stoich>();
        this.reversible = false;
        this.segments = new ArrayList<Segment>();
        this.reactionRule = "";
        this.active = ActiveDirections.FORWARD;
    }

    /**
     * Add an alias with the specified label (BiGG ID) and the specified name.
     *
     * @param label2	BiGG ID
     * @param name2		name; if NULL, the BiGG ID will be used
     */
    public void addAlias(String label2, String name2) {
        if (StringUtils.isBlank(name2))
            name2 = label2;
        this.aliases.put(label2, name2);
    }

    /**
     * @return the reaction ID number
     */
    public int getId() {
        return this.id;
    }

    /**
     * @return the reaction's role name
     */
    public String getName() {
        return this.name;
    }

    /**
     * @return the BiGG ID for the reaction
     */
    public String getBiggId() {
        return this.biggId;
    }

    /**
     * @return the rule for triggering the reaction
     */
    public String getReactionRule() {
        return this.reactionRule;
    }

    /**
     * @return TRUE if this reaction is reversible
     */
    public boolean isReversible() {
        return this.reversible;
    }

    /**
     * Determine whether or not a compound is a product.  This only checks its position in the normal
     * form of the reaction.  To determine whether or not the compound is a potential input or output,
     * use "isInput" and "isOutput".
     *
     * @return TRUE if the specified metabolite is a product of this reaction
     *
     * @param compound	metabolite to check
     */
    public boolean isProduct(String compound) {
        boolean retVal = this.metabolites.stream().anyMatch(x -> x.isProduct() && x.getMetabolite().contentEquals(compound));
        return retVal;
    }

    /**
     * @return TRUE if the specified stoichiometric element is a possible output of this reaction
     *
     * @param stoich	stoichiometric element of interest
     */
    public boolean isOutput(Stoich stoich) {
        return this.active.isOutput(stoich);
    }

    /**
     * @return TRUE if the specified stoichiometric element is a possible input of this reaction
     *
     * @param stoich	stoichiometric element of interest
     */
    public boolean isInput(Stoich stoich) {
        return this.active.isInput(stoich);
    }

    /**
     * @return the aliases for the genes that relate to this reaction
     */
    public Set<String> getGenes() {
        Set<String> retVal = new TreeSet<String>(this.aliases.keySet());
        retVal.addAll(this.aliases.values());
        return retVal;
    }

    /**
     * @return the set of triggering genes for this reaction
     */
    public Set<String> getTriggers() {
        String[] parts = DELIMITERS.split(this.reactionRule);
        Set<String> retVal = new TreeSet<String>();
        for (String part : parts) {
            if (! part.isEmpty() && ! CONNECTORS.contains(part.toLowerCase()))
                retVal.add(part);
        }
        return retVal;
    }

    /**
     * @return a reaction rule with the gene IDs translated to gene names
     *
     * @param rule		reaction rule to translate
     * @param model		metabolic model containing the reaction
     */
    public static String getTranslatedRule(String rule, MetaModel model) {
        StringBuilder retVal = new StringBuilder(rule.length() + 20);
        // Tokenize the string.  Note we return delimiters as well as words.
        var tokenizer = new StringTokenizer(rule, DELIM_STRING, true);
        while (tokenizer.hasMoreTokens()) {
            String token = tokenizer.nextToken();
            if (DELIM_STRING.contains(token) || CONNECTORS.contains(token.toLowerCase())) {
                // Delimiters and connectors are copied unchanged.
                retVal.append(token);
            } else {
                // Here we have a trigger ID.
                String name = model.geneNameOf(token);
                retVal.append(name);
            }
        }
        return retVal.toString();
    }

    /**
     * @return the components of the reaction (reactants and products with stoichometric coefficients)
     */
    public List<Stoich> getMetabolites() {
        return this.metabolites;
    }

    /**
     * @return the connection segments of this reaction
     */
    public List<Segment> getSegments() {
        return this.segments;
    }

    /**
     * @return the coordinates of the reaction label
     */
    public Coordinate getLabelLoc() {
        return this.labelLoc;
    }

    /**
     * This rather arcane method returns the reaction formula string in a parsed format.  It will
     * consist of an array of strings.  The first string will be a connector, the second a compound
     * ID, the third another connector, and so on until the end.  The first string may be empty in
     * order to enforce the pattern.  This string list can then be formatted at will be the client
     * into HTML, JavaFX flowed text, or any number of crazy things.
     *
     * @param reverse	TRUE to reverse the reaction
     */
    public List<String> getParsedFormula(boolean reverse) {
        List<String> retVal = new ArrayList<String>(2 * this.metabolites.size());
        // Prime the left side with an empty string.
        retVal.add("");
        // Store the left-side metabolites.
        this.fillParsedFormula(retVal, reverse);
        // Add the connector.
        String connector = (this.reversible ? " <-> " : " --> ");
        retVal.add(connector);
        // Store the right-side metabolites.
        this.fillParsedFormula(retVal, ! reverse);
        // Return the list.
        return retVal;
    }
    /**
     * Add the metabolites on the current side of the formula.
     *
     * @param retVal	output list for the parsed formula
     * @param products	TRUE to add products, FALSE to add reactants
     */
    private void fillParsedFormula(List<String> retVal, boolean products) {
        // For the first metabolite, we don't need a connector.
        boolean first = true;
        // Loop through the metabolites.
        for (Stoich stoich : this.metabolites) {
            if (stoich.isProduct() == products) {
                // Insure we have a preceding connector.
                if (first)
                    first = false;
                else
                    retVal.add(" + ");
                int coeff = stoich.getCoeff();
                if (coeff != 1) {
                    // We have a non-unit coefficient, so add it to the preceding connector.
                    int n = retVal.size() - 1;
                    retVal.set(n, retVal.get(n) + Integer.toString(coeff) + "*");
                }
                // Now add the metabolite itself.
                retVal.add(stoich.biggId);
            }
        }
    }

    /**
     * @return the reaction formula as a string
     *
     * @param reverse	TRUE to reverse the reaction
     */
    public String getFormula(boolean reverse) {
        List<String> formList = this.getParsedFormula(reverse);
        String retVal = StringUtils.join(formList, "");
        return retVal;
    }

    /**
     * @return the reaction formula as a string with full compound names
     *
     * @param reverse	TRUE to reverse the reaction
     * @param model		parent metabolic model
     */
    public String getLongFormula(boolean reverse, MetaModel model) {
        List<String> formList = this.getParsedFormula(reverse);
        for (int i = 1; i < formList.size(); i += 2) {
            String name = "(" + model.getCompoundName(formList.get(i)) + ")";
            formList.set(i, name);
        }
        return StringUtils.join(formList, "");
    }

    @Override
    public int compareTo(Reaction o) {
        return (this.id - o.id);
    }

    @Override
    public String toString() {
        return "Reaction " + this.id + "(" + this.getName() + ")";
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + this.id;
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        Reaction other = (Reaction) obj;
        if (this.id != other.id)
            return false;
        return true;
    }

    /**
     * Here we want to determine what the possible outputs are for a reaction with
     * the specified metabolite as input.
     *
     * @param inputId	metabolite used as input
     *
     * @return a list of the eligible output metabolite elements for this reaction
     */
    public Collection<Stoich> getOutputs(String inputId) {
        // Find the stoichiometric coefficient for the input metabolite.
        Optional<Stoich> input = this.metabolites.stream()
                .filter(x -> x.biggId.equals(inputId)).findAny();
        Collection<Stoich> retVal;
        if (! input.isPresent())
            retVal = Collections.emptyList();
        else {
            double inCoeff = input.get().coefficient;
            retVal = this.metabolites.stream().filter(x -> x.coefficient * inCoeff < 0)
                    .collect(Collectors.toList());
        }
        return retVal;
    }

    /**
     * Store a new reaction rule for this reaction.
     *
     * @param reactionRule2		new reaction rule
     */
    protected void setRule(String reactionRule2) {
        this.reactionRule = reactionRule2;
        this.aliases.clear();
    }

    /**
     * Specify the reversibility of this reaction.
     *
     * @param reversible2	TRUE if the reaction is to be reversible, else FALSE
     */
    public void setReversible(boolean reversible2) {
        this.reversible = reversible2;
        ActiveDirections.setDefault(this);
    }

    /**
     * Add a new stoichiometric component to this reaction.
     *
     * @param coeff			coefficient
     * @param metabolite	BiGG ID of the compound
     */
    protected void addStoich(int coeff, String metabolite) {
        this.metabolites.add(new Stoich(coeff, metabolite));
    }

    /**
     * @return the stoichiometric node for the specified metabolite
     *
     * @param compound		BiGG ID of the relevant metabolite
     */
    public Stoich getStoich(String compound) {
        Optional<Stoich> retVal = this.metabolites.stream().filter(x -> x.getMetabolite().contentEquals(compound))
                .findFirst();
        if (retVal.isEmpty())
            throw new IllegalArgumentException("Could not find compound \"" + compound + "\" in reaction " +
                    this.biggId + ".");
        return retVal.get();
    }

    /**
     * @return a JSON object for this reaction
     */
    public JsonObject toJson() {
        JsonObject retVal = new JsonObject();
        retVal.put("name", this.name);
        retVal.put("bigg_id", this.biggId);
        retVal.put("reversibility", this.reversible);
        retVal.put("gene_reaction_rule", this.reactionRule);
        // Store the label coordinates.
        if (this.labelLoc != null) {
            retVal.put("label_x", this.labelLoc.getX());
            retVal.put("label_y", this.labelLoc.getY());
        }
        // Now we do an array of genes.
        JsonArray genes = new JsonArray();
        for (Map.Entry<String, String> geneEntry : this.aliases.entrySet()) {
            // Check for a missing name.
            String label = geneEntry.getKey();
            String name = geneEntry.getValue();
            JsonObject pair = new JsonObject();
            pair.put("bigg_id", label);
            if (! label.contentEquals(name))
                pair.put("name", name);
            genes.add(pair);
        }
        retVal.put("genes", genes);
        // Next, the reaction stoichiometry.
        JsonArray metabolites = new JsonArray();
        for (Stoich stoich : this.metabolites) {
            JsonObject stoichO = new JsonObject();
            stoichO.put("coefficient", stoich.coefficient);
            stoichO.put("bigg_id", stoich.biggId);
            metabolites.add(stoichO);
        }
        retVal.put("metabolites", metabolites);
        // Finally, the segments.
        JsonObject segments = new JsonObject();
        for (Segment segment : this.segments) {
            int segId = segment.getId();
            JsonObject segObject = segment.toJson();
            segments.put(Integer.toString(segId), segObject);
        }
        retVal.put("segments", segments);
        // Return the reaction.
        return retVal;
    }

    /**
     * @return the active-directions status
     */
    public ActiveDirections getActive() {
        return this.active;
    }

    /**
     * Specify a new active-directions status.
     *
     * @param active 	the active-directions status to set
     */
    public void setActive(ActiveDirections active) {
        this.active = active;
    }

    /**
     * This is a comparator for sorting reactions by name.
     */
    public static class ByName implements Comparator<Reaction> {

        @Override
        public int compare(Reaction o1, Reaction o2) {
            int retVal = o1.name.compareTo(o2.name);
            if (retVal == 0)
                retVal = o1.id - o2.id;
            return retVal;
        }

    }

    /**
     * This method computes the weight of a reaction.  Reactions are weighted by the sum of the
     * weights of their inputs (for a branching reaction) or outputs (for a main-line reaction).
     * Here, the user specifies the weighting hash and whether we are weighting reactants or
     * products.
     *
     * @param	weights		map of compound IDs to ratings
     * @param	products	TRUE if we are weighting by products, FALSE if by inputs
     *
     * @return the total weight of the reaction
     */
    public double getWeight(Map<String, CompoundRating> weights, boolean products) {
        double retVal = 0.0;
        for (Stoich stoich : this.metabolites) {
            if (stoich.isProduct() == products) {
                // Here we have a compound of interest.
                CompoundRating rating = weights.getOrDefault(stoich.biggId, CompoundRating.IRRELEVANT);
                retVal += stoich.getCoeff() * rating.getWeight();
            }
        }
        return retVal;
    }

}
