package test.rootfreqs;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.FilteredAlignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.Nucleotide;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.Assert;
import org.junit.jupiter.api.Test;
import phylonco.beast.evolution.errormodel.ErrorModel;
import phylonco.beast.evolution.errormodel.ErrorModelBase;
import phylonco.beast.evolution.likelihood.BeagleTreeLikelihoodWithError;
import test.beast.BEASTTestCase;

import static org.junit.Assert.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class BeagleTreeLikelihoodWithErrorTest extends TreeLikelihoodWithRootSequenceTest {

    @Override
    protected GenericTreeLikelihood newTreeLikelihood() {
        System.setProperty("java.only","false");
        return new rootfreqs.TreeLikelihoodWithError();
    }

    @Test
    public void testRootFreqsWithErrorModelEpsilonZero() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);

        Nucleotide datatype = new Nucleotide();
        ErrorModelBase errorModel = new ErrorModelBase();
        errorModel.initByName("epsilon", "0.0", "datatype", datatype);
        errorModel.initAndValidate();

        // test data
        Alignment data0 = BEASTTestCase.getAlignment();
        FilteredAlignment data = new FilteredAlignment();
        data.initByName("filter", "1-10", "data", data0);
        Tree tree = BEASTTestCase.getTree(data);

        double logP = 0;
        Sequence seq = new Sequence();
        seq.initByName("value", "AAAAAAAAAA", "taxon", "root", "totalcount", 4);
        GenericTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data,
                "tree", tree,
                "siteModel", siteModel,
                "rootfreqseq", seq,
                "errorModel", errorModel);
        logP = likelihood.calculateLogP();

        double logP0 = -42.99453211866883;
        assertEquals(logP0, logP, BEASTTestCase.PRECISION);

        System.out.println("logP with ErrorModel = " + logP);
    }

    @Test
    public void testRootFreqsWithErrorModelEpsilons() throws Exception {
        // test for varying levels of epsilon
        double[] epsilons = {0.01, 0.1, 0.2, 0.3, 0.4, 0.5};
        double[] logP0 = {-43.072386980977484, -45.192132962314105, -49.33853004133921,
                -54.786093240560405, -61.33213107081169, -68.40276833565605};
        for (int i = 0; i < epsilons.length; i++) {
            Double epsilon = Double.valueOf(epsilons[i]);
            testErrorModel(epsilon, logP0[i]);
        }
    }

    private void testErrorModel(Double epsilon, double logP0) throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);

        Nucleotide datatype = new Nucleotide();
        ErrorModelBase errorModel = new ErrorModelBase();
        errorModel.initByName("epsilon", epsilon.toString(), "datatype", datatype);
        errorModel.initAndValidate();

        // test data
        Alignment data0 = BEASTTestCase.getAlignment();
        FilteredAlignment data = new FilteredAlignment();
        data.initByName("filter", "1-10", "data", data0);
        Tree tree = BEASTTestCase.getTree(data);

        double logP = 0;
        Sequence seq = new Sequence();
        seq.initByName("value", "AAAAAAAAAA", "taxon", "root", "totalcount", 4);
        GenericTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data,
                "tree", tree,
                "siteModel", siteModel,
                "rootfreqseq", seq,
                "errorModel", errorModel);
        logP = likelihood.calculateLogP();

        assertEquals(logP0, logP, BEASTTestCase.PRECISION);

//        System.out.println("logP with ErrorModel = " + logP);
    }

    @Test
    public void testJCLikelihoodSmallWithErrorNoRootseq() {
        // test error model likelihood without rootseq
        Alignment data = new Alignment();
        Sequence seqA = new Sequence("a", "A");
        Sequence seqB = new Sequence("b", "A");
        data.initByName(
                "sequence", seqA,
                "sequence", seqB,
                "dataType", "nucleotide"
        );

        TreeParser tree = getTree(data);

        JukesCantor subsModel = new JukesCantor();
        subsModel.initAndValidate();

        SiteModel siteModel = getSiteModel(subsModel);

        Nucleotide datatype = new Nucleotide();

        ErrorModelBase errorModel = new ErrorModelBase();
        errorModel.initByName("epsilon", "0.1", "datatype", datatype);
        errorModel.initAndValidate();

        double logP = getLogLikelihood(data, tree, siteModel, errorModel);
        double expectedLogP = -2.3063595712034233;
        Assert.assertEquals(expectedLogP, logP, BEASTTestCase.PRECISION);
    }

    private TreeParser getTree(Alignment data) {
        TreeParser tree = new TreeParser();
        tree.initByName(
                "taxa", data,
                "newick", "(a: 0.5, b: 0.5);",
                "IsLabelledNewick", true
        );
        return tree;
    }


    private SiteModel getSiteModel(SubstitutionModel subsModel) {
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", subsModel);
        siteModel.initAndValidate();
        return siteModel;
    }


    private double getLogLikelihood(Alignment data, TreeParser tree,
                                           SiteModel siteModel, ErrorModel errorModel) {
        GenericTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName(
                "data", data,
                "tree", tree,
                "siteModel", siteModel,
                "useAmbiguities", true,
                "useTipLikelihoods", true,
                "errorModel", errorModel);

        return likelihood.calculateLogP();
    }


    private void testLikelihood(SiteModel siteModel, double logP0) throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        Alignment data0 = BEASTTestCase.getAlignment();
        FilteredAlignment data = new FilteredAlignment();
        data.initByName("filter", "1-10", "data", data0);
        Tree tree = BEASTTestCase.getTree(data);

        // test with rootfreqseq input, normal sequence
        double logP = 0;
        Sequence seq = new Sequence();
        seq.initByName("value", "AAAAAAAAAA", "taxon", "root", "totalcount", 4);
        GenericTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq);
        logP = likelihood.calculateLogP();
        assertEquals(logP0, logP, BEASTTestCase.PRECISION);

        // test with forced scaling
        likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq, "scaling", "always");
        logP = likelihood.calculateLogP();
        assertEquals(logP0, logP, BEASTTestCase.PRECISION);

        // test with rootfreqseq input, uncertain sequence
        seq = new Sequence();
        seq.initByName("value", "1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;", "uncertain", true, "taxon", "root", "totalcount", 4);
        likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq);
        logP = likelihood.calculateLogP();
        assertEquals(logP0, logP, BEASTTestCase.PRECISION);

        // test with forced scaling
        likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "rootfreqseq", seq, "scaling", "always");
        logP = likelihood.calculateLogP();
        assertEquals(logP0, logP, BEASTTestCase.PRECISION);
    }

}
