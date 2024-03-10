package test.rootfreqs;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.FilteredAlignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.tree.Tree;
import org.junit.jupiter.api.Test;
import test.beast.BEASTTestCase;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class BeagleTreeLikelihoodTest  extends TreeLikelihoodWithRootSequenceTest {

    @Override
    protected GenericTreeLikelihood newTreeLikelihood() {
        System.setProperty("java.only","false");
        return new rootfreqs.TreeLikelihood();
    }

}
