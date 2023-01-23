
package test.rootfreqs;

import beast.base.evolution.likelihood.GenericTreeLikelihood;

public class ThreadedTreeLikelihoodWithRootFreqsTest extends TreeLikelihoodWithRootFreqsTest {

	@Override
	protected GenericTreeLikelihood newTreeLikelihood() {
		return new rootfreqs.ThreadedTreeLikelihood();
	}
	
}
