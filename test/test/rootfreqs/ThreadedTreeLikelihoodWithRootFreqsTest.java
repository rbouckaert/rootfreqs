
package rootfreqs;

import beast.base.evolution.likelihood.GenericTreeLikelihood;

public class ThreadedTreeLikelihoodWithRootFreqsTest extends rootfreqs.TreeLikelihoodWithRootFreqsTest {

	@Override
	protected GenericTreeLikelihood newTreeLikelihood() {
		return new rootfreqs.ThreadedTreeLikelihood();
	}
	
}
