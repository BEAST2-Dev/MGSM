package mgsm.evolution.sitemodel;


import java.io.PrintStream;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import mgsm.evolution.likelihood.MGSMBeagleTreeLikelihood;
import beast.base.core.Log;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.distribution.ParametricDistribution;

@Description("Uncorrelated Relaxed Gamma site model that allows different shape parameters for each branche")
public class RelaxedGammaSiteModel extends SiteModel implements Loggable, Function  {
	//public Input<RealParameter> shapesParameterInput = new Input<RealParameter>("shapes", "array of shapes, one for each branch", Validate.REQUIRED);

	public Input<ParametricDistribution> rateDistInput = new Input<ParametricDistribution>("distr", 
			"the distribution governing the rates among branches. " +
			"Unlike UC Relaxed clock model distributions, does not need to have a mean of 1.", Input.Validate.REQUIRED);
    public Input<IntegerParameter> categoryInput = new Input<IntegerParameter>("rateCategories", 
    		"the rate categories associated with nodes in the tree for sampling of individual " +
    		"rates among branches.", Input.Validate.REQUIRED);
    public Input<Tree> treeInput = new Input<Tree>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);

    double [][] categoryRates;
	
	// member vars for relaxed part 
    ParametricDistribution distribution;
    IntegerParameter categories;
    Tree tree;

    private boolean recompute = true;

    private double[] rates;
    private double[] storedRates;
	
	
	public RelaxedGammaSiteModel() {
		shapeParameterInput.setRule(Validate.FORBIDDEN);
	}
	
	@Override
	public void initAndValidate() {
        //shapesParameter = shapesParameterInput.get();
        super.initAndValidate();
        
        categoryCount = gammaCategoryCount.get();
        if (categoryCount <= 1) {
            throw new RuntimeException("RelaxedGammaSiteModel: Invalid category count (" + categoryCount + "). Should be at least 2");
        }

        boolean hasMGSMBeagleLikelihood = false;
        for (Object o : getOutputs()) {
        	if (o instanceof MGSMBeagleTreeLikelihood) {
        		hasMGSMBeagleLikelihood = true;
        		break;
        	}
        }
        if (!hasMGSMBeagleLikelihood) {
	        if (!Boolean.valueOf(System.getProperty("java.only"))) {
	        	Log.warning.println("MultiGammaSiteModel is java-only -- unless you use MGSMBeagleTreeLikelihood");
	        	System.setProperty("java.only", "true");
	        }
        }
	
        tree = treeInput.get();

        categories = categoryInput.get();
        int nCategoryCount = tree.getNodeCount() - 1;
        categories.setDimension(nCategoryCount);
        Integer[] iCategories = new Integer[nCategoryCount];
        for (int i = 0; i < nCategoryCount; i++) {
            iCategories[i] = i;
        }
        IntegerParameter other = new IntegerParameter(iCategories);
        categories.assignFromWithoutID(other);
        categories.setLower(0);
        categories.setUpper(categories.getDimension() - 1);

        distribution = rateDistInput.get();

        rates = new double[categories.getDimension()];
        storedRates = new double[categories.getDimension()];
        for (int i = 0; i < rates.length; i++) {
            try {
				rates[i] = distribution.inverseCumulativeProbability((i + 0.5) / rates.length);
			} catch (MathException e) {
				throw new IllegalArgumentException(e);
			}
        }
        System.arraycopy(rates, 0, storedRates, 0, rates.length);
        //normalize = normalizeInput.get();

        hasPropInvariantCategory = false;
        if (invarParameter.getValue() > 0) {
        	categoryCount += 1;
        	hasPropInvariantCategory = true;
        	//throw new RuntimeException("MultiGammaSiteModel: Cannot handle invariant sites, please set prop invariant to zero");
        }

	}
	
	@Override
    /**
     * discretisation of gamma distribution with equal proportions in each
     * category
     * @param node
     */
    protected void calculateCategoryRates(final Node nodeX) {
        double propVariable = 1.0;
        int cat = 0;

        // 	deal with invarParameter.getValue() > 0
        if (hasPropInvariantCategory) {
            // categoryRates[0][x] = 0.0;
            categoryProportions[0] = invarParameter.getValue();
            cat = 1;
        }
        propVariable = 1.0 - invarParameter.getValue();


        final int gammaCatCount = categoryCount - cat;

        // calc proportions
        for (int i = 0; i < gammaCatCount; i++) {
            categoryProportions[i + cat] = propVariable / gammaCatCount;
        }
        
        if (nodeX == null || nodeX.isRoot()) {
        	return;
        }

    	Tree tree = nodeX.getTree();
        if (categoryRates == null) {
        	int n = tree.getNodeCount();
        	categoryRates = new double[n][cat + gammaCatCount];
        }
        
//        if (node.getParent().isRoot()) {
//        	Node root = node.getParent();
//        	a = (shapesParameter.getValue(root.getLeft().getNr()) + shapesParameter.getValue(root.getLeft().getNr()))/ 2.0;
//        } else {
//        	a = shapesParameter.getValue(node.getNr());
//        }
        
        for (int k= 0; k < categoryRates.length; k++) {
            double mean = 0.0;
            double a = 0.0;
        	Node node = tree.getNode(k);
        	if (!node.isRoot()) {
       			a = getShapeParameterForBranch(node);//shapesParameter.getValue(node.getNr());
        		final GammaDistribution g = new GammaDistributionImpl(a, 1.0 / a);
        	
		        for (int i = 0; i < gammaCatCount; i++) {
		            try {
		                // RRB: alternative implementation that seems equally good in
		                // the first 5 significant digits, but uses a standard distribution object
		            	if (useBeast1StyleGamma) {
		                    categoryRates[k][i + cat] = GammaDistributionQuantile((2.0 * i + 1.0) / (2.0 * gammaCatCount), a, 1.0 / a);
		            	} else {
		            		categoryRates[k][i + cat] = g.inverseCumulativeProbability((2.0 * i + 1.0) / (2.0 * gammaCatCount));
		            	}
		
		            } catch (Exception e) {
		                e.printStackTrace();
		                System.err.println("Something went wrong with the gamma distribution calculation");
		                System.exit(-1);
		            }
		            mean += categoryRates[k][i + cat];
		
		        }
		        mean = (propVariable * mean) / gammaCatCount;
		        
		        for (int i = 0; i < gammaCatCount; i++) {
		            categoryRates[k][i + cat] /= mean;
		        }
        	}
        }

        ratesKnown = true;
    }
	
	
    public double getShapeParameterForBranch(Node node) {
        if (node.isRoot()) {
            // root has no rate
            return 1;
        }
        if (recompute) {
            prepare();
            recompute = false;
        }

        int nodeNumber = node.getNr();

        if (nodeNumber == categories.getDimension()) {
            // root node has nr less than #categories, so use that nr
            nodeNumber = node.getTree().getRoot().getNr();
        }

        int rateCategory = categories.getValue(nodeNumber);

        return rates[rateCategory];// * scaleFactor * 1.0;//meanRate.getValue();
    }


    private void prepare() {
        //System.out.println("prepare");

        categories = (IntegerParameter) categoryInput.get();

        distribution = rateDistInput.get();

        tree = treeInput.get();

        rates = new double[categories.getDimension()];
        try {
            for (int i = 0; i < rates.length; i++) {
                rates[i] = distribution.inverseCumulativeProbability((i + 0.5) / rates.length);
            }
        } catch (Exception e) {
            // Exception due to distribution not having  inverseCumulativeProbability implemented.
            // This should already been caught at initAndValidate()
            e.printStackTrace();
            System.exit(0);
        }
    }
	
	
    @Override
    protected void refresh() {
    	categoryCount = gammaCategoryCount.get();

        if (/*invarParameter != null && */invarParameter.getValue() > 0) {
            if (hasPropInvariantCategory) {
                categoryCount += 1;
            }
        }

        categoryRates = null;
        categoryProportions = new double[categoryCount];
        //calculateCategoryRates(null);
    }


    @Override
    public double getRateForCategory(final int category, final Node node) {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates(node);
            }
        }

        //final double mu = (muParameter != null) ? muParameter.getValue() : 1.0;

        return categoryRates[node.getNr()][category] * muParameter.getValue();
    }

    @Override
    public double[] getCategoryRates(final Node node) {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates(node);
            }
        }

        final double mu = muParameter.getValue();//(muParameter != null) ? muParameter.getValue() : 1.0;

        final double[] rates = new double[categoryCount];
        if (ratesKnown) {
        	for (int i = 0; i < rates.length; i++) {
        		rates[i] = categoryRates[node.getNr()][i] * mu;
        	}
        }

        return rates;
    }

    @Override
    protected boolean requiresRecalculation() {
        // we only get here if something is dirty in its inputs, so always return true
        ratesKnown = false;

        
        recompute = false;
        //renormalize = true;

//        if (treeInput.get().somethingIsDirty()) {
//        	recompute = true;
//            return true;
//        }
        // rateDistInput cannot be dirty?!?
        if (rateDistInput.get().isDirtyCalculation()) {
            recompute = true;
            return true;
        }
        // NOT processed as trait on the tree, so DO mark as dirty
        if (categoryInput.get().somethingIsDirty()) {
            //recompute = true;
            return true;
        }
//        if (meanRate.somethingIsDirty()) {
//            return true;
//        }

        return true;
    }

    @Override
    public void store() {
        System.arraycopy(rates, 0, storedRates, 0, rates.length);
        super.store();
    }

    @Override
    public void restore() {
        double[] tmp = rates;
        rates = storedRates;
        storedRates = tmp;
        super.restore();
    }

	@Override
	public void init(PrintStream out) {
		for (int i = 0; i < tree.getNodeCount() - 1; i++) {
			out.append("shape" + (i+1) + "\t");
		}
	}

	@Override
	public void log(long nSample, PrintStream out) {
		Node [] node = tree.getNodesAsArray();
		for (int i = 0; i < tree.getNodeCount() - 1; i++) {
			out.append(getShapeParameterForBranch(node[i]) + "\t");
		}
	}

	@Override
	public void close(PrintStream out) {
		// nothing to do
	}

	@Override
	public int getDimension() {
		return tree.getNodeCount();
	}

	@Override
	public double getArrayValue() {
		return 0;
	}

	@Override
	public double getArrayValue(int iDim) {
		Node [] node = tree.getNodesAsArray();
		return getShapeParameterForBranch(node[iDim]);
	}
	
    //@Override
    public boolean hasNodeIndependentCategories() {
    	return false;
    }
}
