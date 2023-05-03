package test.rootfreqs;


import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.tree.Tree;
import rootfreqs.TreeLikelihood;

import test.beast.BEASTTestCase;
import test.beast.evolution.alignment.UncertainAlignmentTest;

/**
 * This test mimics the testLikelihood.xml file from Beast 1, which compares Beast 1 results to PAUP results.
 * So, it these tests succeed, then Beast II calculates the same for these simple models as Beast 1 and PAUP.
 * *
 */
public class TreeLikelihoodWithRootFreqsTest  {

    public TreeLikelihoodWithRootFreqsTest() {
        super();
    }

    protected GenericTreeLikelihood newTreeLikelihood() {
    	System.setProperty("java.only","true");
        return new rootfreqs.TreeLikelihood();
    }

    @Test
    public void testRootFreqsLikelihood() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        Alignment data = BEASTTestCase.getAlignment();
        Tree tree = BEASTTestCase.getTree(data);

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);

        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", "0.25 0.25 0.25 0.25",
                "estimate", false);

        GenericTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootFrequencies", freqs);
        double logP = 0;
        logP = likelihood.calculateLogP();
        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);



        likelihood.initByName("useAmbiguities", true, "data", data, "tree", tree, "siteModel", siteModel, "rootFrequencies", freqs);
        logP = likelihood.calculateLogP();
        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);



        freqs.initByName("frequencies", "0.35 0.15 0.25 0.25",
                "estimate", false);
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootFrequencies", freqs);
        logP = likelihood.calculateLogP();
        assertFalse(Math.abs(logP - -1992.2056440317247) < BEASTTestCase.PRECISION);
        assertEquals(logP, -1999.8294984023041, BEASTTestCase.PRECISION);



    }


    @Test
    public void testJC69LikelihoodWithUncertainCharacters() throws Exception {

    	Alignment data = UncertainAlignmentTest.getAlignment();
    	Alignment data2 = UncertainAlignmentTest.getUncertainAlignment();
    	double[] logL, logL_uncertain;

    	System.out.println("\nTree A:");
    	Tree tree = UncertainAlignmentTest.getTreeA(data2);
    	logL = testJC69Likelihood(data,tree);
    	logL_uncertain = testJC69Likelihood(data2,tree);
    	double x1 = -11.853202336328778;
    	double x2 = -12.069603116476458;
    	double x2f = -12.167015830399428;
    	assertEquals(logL[0], x1, BEASTTestCase.PRECISION);
    	assertEquals(logL[1], x1, BEASTTestCase.PRECISION);
    	assertEquals(logL_uncertain[0], x1, BEASTTestCase.PRECISION);
    	assertEquals(logL_uncertain[1], x2, BEASTTestCase.PRECISION);
    	assertEquals(logL_uncertain[2], x2, BEASTTestCase.PRECISION);
    	assertEquals(logL_uncertain[3], x2f, BEASTTestCase.PRECISION);

    	System.out.println("\nTree B:");
    	tree = UncertainAlignmentTest.getTreeB(data2);
    	logL = testJC69Likelihood(data,tree);
    	logL_uncertain = testJC69Likelihood(data2,tree);
    	double x3 = -12.421114302827698;
    	double x4 = -11.62105662310513;
    	double x4f = -11.703443703870239;
    	assertEquals(logL[0], x3, BEASTTestCase.PRECISION);
    	assertEquals(logL[1], x3, BEASTTestCase.PRECISION);
    	assertEquals(logL_uncertain[0], x3, BEASTTestCase.PRECISION);
    	assertEquals(logL_uncertain[1], x4, BEASTTestCase.PRECISION);
    	assertEquals(logL_uncertain[2], x4, BEASTTestCase.PRECISION);
    	assertEquals(logL_uncertain[3], x4f, BEASTTestCase.PRECISION);

    	System.out.println("\nTesting alignment doubling:");
    	Alignment data3 = UncertainAlignmentTest.getUncertainAlignmentDoubled();
    	logL_uncertain = testJC69Likelihood(data3,tree);
    	assertEquals(logL_uncertain[0], 2 * x3, BEASTTestCase.PRECISION);
    	assertEquals(logL_uncertain[1], 2 * x4, BEASTTestCase.PRECISION);
    	assertEquals(logL_uncertain[2], 2 * x4, BEASTTestCase.PRECISION);
    	assertEquals(logL_uncertain[3], 2 * x4f, BEASTTestCase.PRECISION);

    }

    public double[] testJC69Likelihood(Alignment data, Tree tree) throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "0.6", "substModel", JC);
        // NB The rate in the JC model used here is actually alpha * 3 in the usual sense, because
        // it's divided by 3 before multiplying in the exponent (not sure why)

        System.out.println("Without tip likelihoods:");
        GenericTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "scaling", TreeLikelihood.Scaling.none);
        double[] logP = new double[4];
        logP[0] = likelihood.calculateLogP();
        System.out.println(logP[0]);

        System.out.println("With tip likelihoods:");
        likelihood = newTreeLikelihood();
        likelihood.initByName("useTipLikelihoods", true, "data", data, "tree", tree, "siteModel", siteModel, "scaling", TreeLikelihood.Scaling.none);
        logP[1]= likelihood.calculateLogP();
        System.out.println(logP[1]);

        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", "0.25 0.25 0.25 0.25",
                "estimate", false);

        System.out.println("With tip likelihoods and equal root frequencies:");
        likelihood = newTreeLikelihood();
        likelihood.initByName("useTipLikelihoods", true, "data", data, "tree", tree, "siteModel", siteModel, "scaling", TreeLikelihood.Scaling.none, "rootFrequencies", freqs);
        logP[2]= likelihood.calculateLogP();
        System.out.println(logP[2]);


        freqs.initByName("frequencies", "0.35 0.15 0.25 0.25",
                "estimate", false);

        System.out.println("With tip likelihoods and unequal root frequencies:");
        likelihood = newTreeLikelihood();
        likelihood.initByName("useTipLikelihoods", true, "data", data, "tree", tree, "siteModel", siteModel, "scaling", TreeLikelihood.Scaling.none, "rootFrequencies", freqs);
        logP[3]= likelihood.calculateLogP();
        System.out.println(logP[3]);

        return logP;
    }

} // class TreeLikelihoodTest
