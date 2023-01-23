package test.rootfreqs;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;

import org.junit.jupiter.api.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.FilteredAlignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.tree.Tree;
import test.beast.BEASTTestCase;

public class TreeLikelihoodWithRootSequenceTest {

    protected GenericTreeLikelihood newTreeLikelihood() {
    	System.setProperty("java.only","true");
        return new rootfreqs.TreeLikelihood();
    }

    @Test
    public void testRootFreqsLikelihood() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        Alignment data0 = BEASTTestCase.getAlignment();
        FilteredAlignment data = new FilteredAlignment();
        data.initByName("filter", "1-10", "data", data0);
        Tree tree = BEASTTestCase.getTree(data);

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);

        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", "1.0 0.0 0.0 0.0",
                "estimate", false);

        GenericTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootFrequencies", freqs);
        double logP = 0;
        double logP0 = -42.99453211866883;
        logP = likelihood.calculateLogP();
        assertEquals(logP, logP0, BEASTTestCase.PRECISION);

        Sequence seq = new Sequence();
        seq.initByName("value", "AAAAAAAAAA", "taxon", "root", "totalcount", 4);
        likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq);
        logP = likelihood.calculateLogP();
        assertEquals(logP, logP0, BEASTTestCase.PRECISION);
        
        seq = new Sequence();
        seq.initByName("value", "1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;", "uncertain", true, "taxon", "root", "totalcount", 4);
        likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq);
        logP = likelihood.calculateLogP();
        assertEquals(logP, logP0, BEASTTestCase.PRECISION);
    }
    
    
    
}
