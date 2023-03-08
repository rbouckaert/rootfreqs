package test.rootfreqs;

import static org.junit.Assert.assertFalse;
import static org.junit.jupiter.api.Assertions.assertEquals;

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
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);

        double logP0 = -42.99453211866883;
        testLikelihood(siteModel, logP0);
    }
    
    @Test
    public void testRootFreqsLikelihoodWithRateHeterogeneity() throws Exception {
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 4, "shape", "1.0", "substModel", JC);
        double logP0 = -41.55016128066514;
        testLikelihood(siteModel, logP0);
    }
    	
    	
    private void testLikelihood(SiteModel siteModel, double logP0) throws Exception {
    	
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        Alignment data0 = BEASTTestCase.getAlignment();
        FilteredAlignment data = new FilteredAlignment();
        data.initByName("filter", "1-10", "data", data0);
        Tree tree = BEASTTestCase.getTree(data);

        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", "1.0 0.0 0.0 0.0",
                "estimate", false);

        // test with rootFrequencies input
        GenericTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootFrequencies", freqs);
        double logP = 0;
        logP = likelihood.calculateLogP();
        assertEquals(logP, logP0, BEASTTestCase.PRECISION);

        // test with rootfreqseq input, normal sequence
        Sequence seq = new Sequence();
        seq.initByName("value", "AAAAAAAAAA", "taxon", "root", "totalcount", 4);
        likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq);
        logP = likelihood.calculateLogP();
        assertEquals(logP, logP0, BEASTTestCase.PRECISION);

        // test with forced scaling
        likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq, "scaling", "always");
        logP = likelihood.calculateLogP();
        assertEquals(logP, logP0, BEASTTestCase.PRECISION);
        
        
        
        // test with rootfreqseq input, uncertain sequence
        seq = new Sequence();
        seq.initByName("value", "1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;", "uncertain", true, "taxon", "root", "totalcount", 4);
        likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq);
        logP = likelihood.calculateLogP();
        assertEquals(logP, logP0, BEASTTestCase.PRECISION);

        // test with forced scaling
        likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq, "scaling", "always");
        logP = likelihood.calculateLogP();
        assertEquals(logP, logP0, BEASTTestCase.PRECISION);
    }

    
    @Test
    public void testRootFreqsLikelihood2() throws Exception {
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);
        double logP0 = -20.30112275493719;
        testLikelihood2(siteModel, logP0);
    }
    
    @Test
    public void testRootFreqsLikelihoodWithRateHeterogeneity2() throws Exception {
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 4, "shape", "1.0", "substModel", JC);
        double logP0 = -20.1372796693738;
        testLikelihood2(siteModel, logP0);
    }

    private void testLikelihood2(SiteModel siteModel, double logP0) throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        Alignment data0 = BEASTTestCase.getAlignment();
        FilteredAlignment data = new FilteredAlignment();
        data.initByName("filter", "1-10", "data", data0);
        Tree tree = BEASTTestCase.getTree(data);


        // test with default root frequencies from site model, which for JC96 is [1/4,1/4,1/4,1/4]
        GenericTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);
        double logP = 0;
        logP = likelihood.calculateLogP();
        assertEquals(logP, logP0, BEASTTestCase.PRECISION);

        // test with rootfreqseq using ambiguous character sequence
        Sequence seq = new Sequence();
        seq.initByName("value", "?????-----", "taxon", "root", "totalcount", 4);
        likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq);
        logP = likelihood.calculateLogP();
        assertEquals(logP, logP0, BEASTTestCase.PRECISION);
        
        // test with rootfreqseq using uncertain sequence
        seq = new Sequence();
        seq.initByName("value", "0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;", "uncertain", true, "taxon", "root", "totalcount", 4);
        likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq);
        logP = likelihood.calculateLogP();
        assertEquals(logP, logP0, BEASTTestCase.PRECISION);
    }
    
    
    @Test
    public void testRootFreqsLikelihood3() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);
        double logP0 = -20.30112275493719;
        testLikelihood3(siteModel, logP0);
    }

    @Test
    public void testRootFreqsLikelihoodWithRateHeterogeneity3() throws Exception {
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 4, "shape", "1.0", "substModel", JC);
        double logP0 = -20.1372796693738;
        testLikelihood3(siteModel, logP0);
    }

    
    private void testLikelihood3(SiteModel siteModel, double logP0) throws Exception {

        Alignment data0 = BEASTTestCase.getAlignment();
        FilteredAlignment data = new FilteredAlignment();
        data.initByName("filter", "1-10", "data", data0);
        Tree tree = BEASTTestCase.getTree(data);


        // test with default root frequencies from site model, which for JC96 is [1/4,1/4,1/4,1/4]
        GenericTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);
        double logP = 0;
        logP = likelihood.calculateLogP();
        assertEquals(logP, logP0, BEASTTestCase.PRECISION);

        // test with rootfreqseq using ambiguous character sequence
        Sequence seq = new Sequence();
        seq.initByName("value", "T????-----", "taxon", "root", "totalcount", 4);
        likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq);
        logP = likelihood.calculateLogP();
        assertFalse(Math.abs(logP - logP0)< BEASTTestCase.PRECISION);
        
        // test with rootfreqseq using uncertain sequence
        seq = new Sequence();
        seq.initByName("value", "0,0,0,1;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;", "uncertain", true, "taxon", "root", "totalcount", 4);
        likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq);
        double logP2 = likelihood.calculateLogP();
        assertFalse(Math.abs(logP2 - logP0)< BEASTTestCase.PRECISION);
        assertEquals(logP, logP2, BEASTTestCase.PRECISION);
    }
    
    @Test
    public void testRootFreqsLikelihoodWithInvariableCategory() throws Exception {
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", 
        		"gammaCategoryCount", 1, 
        		"shape", "1.0",
        		"proportionInvariant", "0.25",
        		"substModel", JC);
        // NB: there are constant sites in the alignment
        double logP0 = -20.203302531716012;
        testLikelihood3(siteModel, logP0);
    }
     
}
