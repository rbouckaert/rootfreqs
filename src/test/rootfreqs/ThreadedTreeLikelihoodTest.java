package test.rootfreqs;

import beast.base.evolution.likelihood.GenericTreeLikelihood;

public class ThreadedTreeLikelihoodTest extends TreeLikelihoodTest {

	@Override
	protected GenericTreeLikelihood newTreeLikelihood() {
		return new rootfreqs.ThreadedTreeLikelihood();
	}
	
}
