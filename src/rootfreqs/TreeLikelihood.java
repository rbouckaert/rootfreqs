package rootfreqs;

import java.util.List;

import beagle.Beagle;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.sitemodel.SiteModelInterface.Base;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

@Description("Tree-likelihood that allows site specific root frequencies")
public class TreeLikelihood extends beast.base.evolution.likelihood.TreeLikelihood {

	final public Input<Sequence> rootFrequenciesSequenceInput = new Input<>("rootfreqseq", 
			  "If defined, specifies site specific root frequencies instead of uniform root frequencies. "
			+ "If it is a standard sequence, root frequencies will be set at 1 for the observed value in the sequence "
			+ "(or uniform if an ambiguous value) and all other frequencies will be set at 0. "
			+ "If the sequence is uncertain, each site represents a root frequencies distribution over all states.");

	private double [][] rootFrequenciesSequence;
	private int siteCount, stateCount, categoryCount, patternCount;


	@Override
	public void initAndValidate() {
		super.initAndValidate();
			
		if (rootFrequenciesSequenceInput.get() != null) {
			
			// sanity check
			if (rootFrequenciesInput.get() != null) {
				throw new IllegalArgumentException("Either rootFrequencies or rootfreqseq can be specified, but not both");
			}
			
			
			initRootFrequencies();
			
			// sanity check
			if (siteCount != rootFrequenciesSequence.length) {
				throw new IllegalArgumentException("root sequence length (" + rootFrequenciesSequence.length + ") differs from alignment length("+ siteCount + ")");
			}
			
			patternLogLikelihoods = new double[siteCount];
			m_siteModel = (Base) siteModelInput.get();
			substitutionModel = m_siteModel.getSubstitutionModel();

		} else {
			rootFrequenciesSequence = null;
		}
	}
	
	

	protected void initRootFrequencies() {
		Sequence seq = rootFrequenciesSequenceInput.get();
		stateCount = seq.totalCountInput.get();
		siteCount = dataInput.get().getSiteCount();
        // account for alignment weights here
		patternCount = dataInput.get().getPatternCount();
		categoryCount = ((SiteModelInterface.Base)siteModelInput.get()).getCategoryCount();
		if (categoryCount <= 0) {
			categoryCount = 1;
		}

		// deal with uncertain root sequence first
		if (seq.uncertainInput.get() != null && seq.uncertainInput.get()) {
			rootFrequenciesSequence = seq.getLikelihoods();
			return;
		}
		
		// it is an ordinairy sequence, not an uncertain one
		// todo: test for ambiguous codes
		rootFrequenciesSequence = new double[siteCount][stateCount];
		DataType dataType = dataInput.get().getDataType();
		List<Integer> values = seq.getSequence(dataType);
		for (int i = 0; i < siteCount; i++) {
			int [] state = dataType.getStatesForCode(values.get(i));
			for (int j : state) {
				rootFrequenciesSequence[i][j] = 1.0/state.length;
			}	
		}
	}

	
	@Override
	public double calculateLogP() {
        if (beagle != null) {
            logP = beagle.calculateLogP();
            if (rootFrequenciesSequence != null) {
            	logP = recalculateBeagleLogPWithRootFrequences();
            }
            return logP;
        }
        logP = super.calculateLogP();;
		return logP;
	}	


	private double [] rootpartials2 = null;
	private double [] rootpartials = null;
	private double recalculateBeagleLogPWithRootFrequences() {
		if (rootpartials == null) {
			rootpartials = new double[patternCount * stateCount];
			if (categoryCount > 1) {
				rootpartials2 = new double[patternCount * stateCount * categoryCount];
			}
		} else if (rootpartials.length != patternCount * stateCount) {
            System.out.println("root partials not equal to num patterns " + rootpartials.length + " != " + (patternCount * stateCount));
			rootpartials = new double[patternCount * stateCount];
			if (categoryCount > 1) {
				rootpartials2 = new double[patternCount * stateCount * categoryCount];
			}
		}
				
		int number = treeInput.get().getRoot().getNr();
		int node = beagle.getPartialBufferHelper().getOffsetIndex(number);
		if (categoryCount <= 1) {
			beagle.getBeagle().getPartials(node, Beagle.NONE, rootpartials);
		} else {
			beagle.getBeagle().getPartials(node, Beagle.NONE, rootpartials2);
			// integrate categories first
            final double[] proportions = ((SiteModelInterface.Base)siteModelInput.get()).getCategoryProportions(treeInput.get().getRoot());
            calculateIntegratePartials(rootpartials2, proportions, rootpartials);
		}
		
		double [] beaglePatternLogLikelihoods = beagle.getPatternLogLikelihoods();
		double [] scaleFactor = new double[patternCount];
		double [] freqs = substitutionModel.getFrequencies();
        int u = 0;
		for (int i = 0; i < patternCount; i++) {
            double sum = 0.0;
            for (int j = 0; j < stateCount; j++) {
                sum += freqs[j] * rootpartials[u];
                u++;
            }
            scaleFactor[i] = beaglePatternLogLikelihoods[i] - Math.log(sum);
		}
		
        for (int k = 0; k < siteCount; k++) {
        	int j = dataInput.get().getPatternIndex(k);
        	int v = j * stateCount;
            double sum = 0.0;
            for (int i = 0; i < stateCount; i++) {

                sum += rootFrequenciesSequence[k][i] * rootpartials[v];
                v++;
            }
            patternLogLikelihoods[k] = Math.log(sum) + scaleFactor[j];
        }
		calcLogP();
		return logP;
	}

	
	protected void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials) {

        int u = 0;
        int v = 0;
        for (int k = 0; k < patternCount; k++) {

            for (int i = 0; i < stateCount; i++) {

                outPartials[u] = inPartials[v] * proportions[0];
                u++;
                v++;
            }
        }


        for (int l = 1; l < categoryCount; l++) {
            u = 0;

            for (int k = 0; k < patternCount; k++) {

                for (int i = 0; i < stateCount; i++) {

                    outPartials[u] += inPartials[v] * proportions[l];
                    u++;
                    v++;
                }
            }
        }
    }
	
    protected void calcLogP() {
    	if (rootFrequenciesSequence != null) {
            logP = 0.0;

            if (dataInput.get().siteWeightsInput.get() != null) {
                // handle alignment with weights
                for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                    logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
                }
            } else {
                // full alignment without weights
                for (int i = 0; i < siteCount; i++) {
                    logP += patternLogLikelihoods[i];
                }
            }
            if (useAscertainedSitePatterns) {
                // at this point, logP contains contributions of ascertained sites, which it should not have.
                final double ascertainmentCorrection = dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
                
                int totalPatternWeight = 0;
                for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                	totalPatternWeight += dataInput.get().getPatternWeight(i);
                }
            	// TODO: test this
                // here, we subtract ascertainmentCorrection once for each site that should contribute to the likelihood (i.e. non-ascertaining sites)
                // and subtract the ascertainmentCorrection for the ascertaining sites already accumulated in logP above 
                logP +=  -ascertainmentCorrection * (totalPatternWeight + 1);
            }
    		
    	} else {
    		super.calcLogP();
    	}
    }

	

	@Override
	protected int traverse(Node node) {
        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        // First update the transition probability matrix(ices) for this branch
        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
            m_branchLengths[nodeIndex] = branchTime;
            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse(child1);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
                }

                if (m_siteModel.integrateAcrossCategories()) {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods

                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

                    if (getConstantPattern() != null) { // && !SiteModel.g_bUseOriginal) {
                        proportionInvariant = m_siteModel.getProportionInvariant();
                        // some portion of sites is invariant, so adjust root partials for this
                        for (final int i : getConstantPattern()) {
                            m_fRootPartials[i] += proportionInvariant;
                        }
                    }

                    if (this.rootFrequenciesSequence != null) {
                    	// use site specific root frequencies
                    	calculateLogLikelihoods(m_fRootPartials, this.rootFrequenciesSequence, patternLogLikelihoods);
                    } else {
                    	double[] rootFrequencies = substitutionModel.getFrequencies();
                    	if (rootFrequenciesInput.get() != null) {
                    		rootFrequencies = rootFrequenciesInput.get().getFreqs();
                    	}
                    	likelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);
                    }
                }

            }
        }
        return update;
    } // traverse
	

	private void calculateLogLikelihoods(
			double[] partials, 
			double[][] rootFrequencies,
			double[] outLogLikelihoods) {

        if (dataInput.get().siteWeightsInput.get() != null) {
            // handle site weights
            int v = 0;
            for (int k = 0; k < patternCount; k++) {
                double sum = 0.0;
                for (int i = 0; i < stateCount; i++) {
                    sum += rootFrequencies[k][i] * partials[v];
                    v++;
                }
                outLogLikelihoods[k] = Math.log(sum) + getLikelihoodCore().getLogScalingFactor(k);
            }
        } else {
            // alignment and root are unweighted
            for (int k = 0; k < siteCount; k++) {
                int j = dataInput.get().getPatternIndex(k);
                int v = j * stateCount;
                double sum = 0.0;
                for (int i = 0; i < stateCount; i++) {

                    sum += rootFrequencies[k][i] * partials[v];
                    v++;
                }
                outLogLikelihoods[k] = Math.log(sum) + getLikelihoodCore().getLogScalingFactor(j);
            }
        }
	}



	@Override
	// defined here to provide access for rootfreqs.ThreaderTreeLikelihood
	protected boolean requiresRecalculation() {
		return super.requiresRecalculation();
	}
}
