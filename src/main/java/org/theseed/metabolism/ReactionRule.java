/**
 *
 */
package org.theseed.metabolism;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;

import java.util.NoSuchElementException;

/**
 * A reaction rule is a structure that represents a reaction rule in disjunctive normal form.  The
 * rule consists of a collection of collections.  The high-level collection is a massive OR.  The
 * lower-level collections are all conjunctions (AND).  Note that the collections are internally
 * represented as sets, because AND and OR are idempotent (A and A == A; A or A == A).  The elements
 * of the inner collections are strings representing BiGG IDs of proteins.
 *
 * @author Bruce Parrello
 *
 */
public class ReactionRule {

    // FIELDS
    /** list of conjunctions */
    private NavigableSet<And> formula;

    /**
     * This class represents a conjunction of proteins
     */
    public static class And extends TreeSet<String> implements Comparable<And> {

        /** serialization version ID */
        private static final long serialVersionUID = -5711234107699161436L;

        /**
         * Create an empty conjunction.
         */
        public And() {
            super();
        }

        /**
         * Create a copy of an existing conjunction.
         *
         * @param source	conjunction to copy
         */
        public And(And source) {
            super(source);
        }

        /**
         * Create a singleton conjunction.
         *
         * @param protein	protein to include
         */
        public And(String protein) {
            super();
            this.add(protein);
        }

        @Override
        public int compareTo(And o) {
            Iterator<String> iter = this.iterator();
            Iterator<String> o_iter = o.iterator();
            int retVal = 0;
            while (retVal == 0 && iter.hasNext()) {
                if (! o_iter.hasNext())
                    retVal = 1;
                else
                    retVal = iter.next().compareTo(o_iter.next());
            }
            if (retVal == 0 && o_iter.hasNext())
                retVal = -1;
            return retVal;
        }

    }

    /**
     * Create a reaction rule for a single protein.
     *
     * @param protein		protein ID
     */
    protected ReactionRule(String protein) {
        this.formula = new TreeSet<And>();
        this.formula.add(new And(protein));
    }

    /**
     * Create a new, empty reaction rule.
     */
    protected ReactionRule() {
        this.formula = new TreeSet<And>();
    }

    /**
     * Construct a reaction rule from a rule string.  The rule string can only contain "and", "or", spaces,
     * parentheses, and protein IDs.
     *
     * @param rule		reaction rule string
     */
    public static ReactionRule parse(String rule) {
        ReactionRule retVal = null;
        // Create a stack of unprocessed rules.
        var outputs = new LinkedList<ReactionRule>();
        // Create the operator stack.
        var ops = new LinkedList<String>();
        // Break the rule into tokens.
        List<String> tokens = tokenize(rule);
        try {
            // Loop through the tokens.
            for (String token : tokens) {
                switch (token) {
                case "(" :
                    ops.push("(");
                    break;
                case "and" :
                    ops.push("and");
                    break;
                case "or" :
                    while (ops.peek() == "and") {
                        processOperator(ops.pop(), outputs);
                    }
                    ops.push("or");
                    break;
                case ")" :
                    // Here we have to unspool everything to the next left paren.
                    while (ops.peek() != "(") {
                        processOperator(ops.pop(), outputs);
                    }
                    // Pop off the left paren.
                    ops.pop();
                    break;
                default :
                    // Here we have a protein name.
                    ReactionRule tokenRule = new ReactionRule(token);
                    outputs.push(tokenRule);
                }
            }
            while (! ops.isEmpty())
                processOperator(ops.pop(), outputs);
            // Return the final rule.
            retVal = outputs.pop();
        } catch (NoSuchElementException e) {
            throw new RuntimeException("Invalid syntax in reaction rule \"" + rule + "\".", e);
        }
        return retVal;
    }

    /**
     * Apply the specified operator to the top two outputs and push the result back on the
     * output stack.
     *
     * @param pop		operator to apply
     * @param outputs	output stack
     */
    private static void processOperator(String pop, LinkedList<ReactionRule> outputs) {
        var right = outputs.pop();
        var left = outputs.pop();
        switch (pop) {
        case "or" :
            left.formula.addAll(right.formula);
            outputs.push(left);
            break;
        case "and" :
            // To AND two rules together, we distribute each conjunction of the right rule across each
            // conjunction on the left.  Start with an empty result.
            ReactionRule result = new ReactionRule();
            for (And clause : left.formula) {
                for (And clause2 : right.formula) {
                    And copy = new And(clause);
                    copy.addAll(clause2);
                    result.formula.add(copy);
                }
            }
            outputs.push(result);
            break;
        default :
            throw new NoSuchElementException();
        }
    }

    /**
     * Separate a reaction rule into tokens.  Each token is either a connector ("and", "or"), a parenthesis,
     * or a protein ID.
     *
     * @param rule		reaction rule to parse
     *
     * @return a list of the tokens in the rule
     */
    private static List<String> tokenize(String rule) {
        final int n = rule.length();
        List<String> retVal = new ArrayList<String>(n / 3);
        StringBuffer token = new StringBuffer(10);
        for (int i = 0; i < n; i++) {
            char c = rule.charAt(i);
            switch (c) {
            case '(' :
                closeToken(retVal, token);
                retVal.add("(");
                break;
            case ')' :
                closeToken(retVal, token);
                retVal.add(")");
                break;
            case ' ' :
                closeToken(retVal, token);
                break;
            default :
                token.append(c);
            }
        }
        // Make sure we have any residual token.
        closeToken(retVal, token);
        return retVal;
    }

    /**
     * This method adds the current token to the token list and prepares for the
     * next one.
     *
     * @param result	token list
     * @param token		current token
     */
    private static void closeToken(List<String> result, StringBuffer token) {
        if (token.length() > 0) {
            result.add(token.toString());
            token.setLength(0);
        }
    }

    /**
     * @return an empty weight map for this rule
     */
    private Map<String, Double> emptyMap() {
        var retVal = new TreeMap<String, Double>();
        this.formula.stream().flatMap(x -> x.stream()).forEach(x -> retVal.put(x, 0.0));
        return retVal;
    }

    /**
     * Compute trigger weights for this rule.  A protein's weight in an and-set is the reciprocal of the
     * set size.  The total weight is its maximum over all the and-sets.  The idea is we want to compute
     * the effect of over-expressing the protein.  Will it stimulate the reaction on its own (weight 1.0)
     * or have a negligible effect (weight close to 0.0)?
     *
     * @return a map of protein IDs to trigger weights for this reaction rule
     */
    public Map<String, Double> getTriggerWeights() {
        Map<String, Double> retVal = this.emptyMap();
        // Loop through the individual And-sets.  The weight of each protein in its set is 1 over
        // the number of proteins in the and-set.
        for (And andSet : this.formula) {
            double weight = 1.0 / andSet.size();
            for (String protein : andSet) {
                double current = retVal.get(protein);
                if (weight > current)
                    retVal.put(protein, weight);
            }
        }
        return retVal;
    }

    /**
     * Compute branching weights for this rule.  A protein's weight in the number of and-sets it occurs
     * in divided by the total number of and-sets.  The idea is we want to compute the effect of
     * knocking out the protein.  Will that suppress the reaction (weight 1.0) or have a negligible effect
     * (weight close to 0.0)?
     *
     * @return a map of protein IDs to branch weights for this reaction rule
     */
    public Map<String, Double> getBranchWeights() {
        Map<String, Double> retVal = this.emptyMap();
        // Compute the weihgt of each occurrence of the protein.
        double weight = 1.0 / this.formula.size();
        // Loop through the individual and-sets.
        for (And andSet : this.formula) {
            for (String protein : andSet)
                retVal.put(protein, retVal.get(protein) + weight);
        }
        return retVal;
    }

    /**
     * @return TRUE if this rule contains only a single protein
     */
    public boolean isSingleton() {
        boolean retVal = this.formula.size() == 1;
        if (retVal)
            retVal = this.formula.first().size() == 1;
        return retVal;
    }

    /**
     * @return the formula as a string
     */
    @Override
    public String toString() {
        String retVal = this.formula.stream().map(x -> StringUtils.join(x, " and "))
                .collect(Collectors.joining(" or "));
        return retVal;
    }

}
