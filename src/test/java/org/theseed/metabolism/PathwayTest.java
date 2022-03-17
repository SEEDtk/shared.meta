/**
 *
 */
package org.theseed.metabolism;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import org.junit.jupiter.api.Test;
import org.theseed.genome.Genome;
import org.theseed.metabolism.mods.ModifierList;
import org.theseed.utils.ParseFailureException;

import com.github.cliftonlabs.json_simple.JsonException;

/**
 * @author Bruce Parrello
 *
 */
class PathwayTest {

    @Test
    void testCompoundRatings() throws IOException, ParseFailureException, JsonException {
        File gFile = new File("data", "MG1655-wild.gto");
        File mFile = new File("data", "ecoli_all.json");
        File tFile = new File("data", "threonine.path.json");
        File nFile = new File("data", "normal.flow");
        Genome genome = new Genome(gFile);
        MetaModel model = new MetaModel(mFile, genome);
        ModifierList flowMods = new ModifierList(nFile);
        flowMods.apply(model);
        Pathway path = new Pathway(tFile, model);
        var ratings = CompoundRating.getRatingMap(path, model);
        var compounds = ratings.keySet();
        assertThat(compounds, hasItem("thr__L_e"));
        assertThat(compounds, hasItem("thr__L_c"));
        assertThat(compounds, hasItem("g6p_c"));
        assertThat(compounds, hasItem("glu__L_c"));
        var ratingList = new ArrayList<CompoundRating>(ratings.values());
        Collections.sort(ratingList);
        try (PrintWriter writer = new PrintWriter(new File("data", "ratings.ser"))) {
            writer.println("compound\tweight\tusage\toutput\tcommon\tgoal\teffort");
            for (CompoundRating rating : ratingList) {
                writer.format("%s\t%8.2f\t%d\t%d\t%b\t%b\t%d%n", rating.getCompoundId(), rating.getWeight(),
                        rating.getUsage(), rating.getOutput(), rating.isCommon(), rating.isGoal(),
                        rating.getEffort());
            }
        }
    }

    @Test
    void TestReactionRule() throws IOException, ParseFailureException, JsonException {
        File gFile = new File("data", "MG1655-wild.gto");
        File mFile = new File("data", "ecoli_all.json");
        File nFile = new File("data", "normal.flow");
        Genome genome = new Genome(gFile);
        MetaModel model = new MetaModel(mFile, genome);
        ModifierList flowMods = new ModifierList(nFile);
        flowMods.apply(model);
        var reactions = model.getAllReactions();
        // Verify that all reactions are parseable and that the singletons work.
        for (Reaction react : reactions) {
            String ruleString = react.getReactionRule();
            if (! ruleString.isBlank()) {
                ReactionRule rule = ReactionRule.parse(ruleString);
                if (! ruleString.contains(" ")) {
                    assertThat(ruleString, rule.isSingleton(), equalTo(true));
                    var weightMap = rule.getBranchWeights();
                    assertThat(ruleString, weightMap.size(), equalTo(1));
                    assertThat(ruleString, weightMap.get(ruleString), equalTo(1.0));
                }
            }
        }
        // Now check some complicated ones to insure they work.
        String ruleString = "(A or (B and C) and (E and F or D) or (B and C and D))";
        ReactionRule rule = ReactionRule.parse(ruleString);
        assertThat(rule.toString(), equalTo("A or B and C and D or B and C and E and F"));
        var branchMap = rule.getBranchWeights();
        assertThat(branchMap.get("A"), closeTo(0.333, 0.001));
        assertThat(branchMap.get("B"), closeTo(0.666, 0.001));
        assertThat(branchMap.get("C"), closeTo(0.666, 0.001));
        assertThat(branchMap.get("D"), closeTo(0.333, 0.001));
        assertThat(branchMap.get("E"), closeTo(0.333, 0.001));
        var triggerMap = rule.getTriggerWeights();
        assertThat(triggerMap.get("A"), closeTo(1.0, 0.001));
        assertThat(triggerMap.get("B"), closeTo(0.333, 0.001));
        assertThat(triggerMap.get("C"), closeTo(0.333, 0.001));
        assertThat(triggerMap.get("D"), closeTo(0.333, 0.001));
        assertThat(triggerMap.get("E"), closeTo(0.250, 0.001));
        assertThat(triggerMap.get("F"), closeTo(0.250, 0.001));
    }

}
