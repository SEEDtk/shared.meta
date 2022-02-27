/**
 *
 */
package org.theseed.metabolism;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.theseed.genome.Genome;

/**
 * @author Bruce Parrello
 *
 */
class TestModelLoad {

    @Test
    void testLoad() throws IOException {
        // Verify that all the reactions and metabolites are found.
        File gFile = new File("data", "MG1655-wild.gto");
        File mFile = new File("data", "ecoli_all.json");
        Genome genome = new Genome(gFile);
        MetaModel model = new MetaModel(mFile, genome);
        assertThat(model.getAllReactions().size(), equalTo(2719));
        assertThat(model.getBReactionMap().size(), equalTo(2713));
        assertThat(model.getMetaboliteCount(), equalTo(1877));
    }

}
