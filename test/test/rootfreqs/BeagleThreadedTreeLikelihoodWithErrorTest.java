package test.rootfreqs;

import beast.base.evolution.likelihood.GenericTreeLikelihood;

public class BeagleThreadedTreeLikelihoodWithErrorTest extends BeagleTreeLikelihoodWithErrorTest {

    @Override
    protected GenericTreeLikelihood newTreeLikelihood() {
        System.setProperty("java.only","false");
        return new rootfreqs.ThreadedTreeLikelihoodWithError();
    }


}
