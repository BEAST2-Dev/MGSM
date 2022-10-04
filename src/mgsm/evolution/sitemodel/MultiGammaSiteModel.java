package mgsm.evolution.sitemodel;

import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import mgsm.evolution.likelihood.MGSMBeagleTreeLikelihood;
import beast.base.core.Log;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

@Description("Gamma site model that allows different shape parameters for each branche")
public class MultiGammaSiteModel extends SiteModel {
	public Input<RealParameter> shapesParameterInput = new Input<RealParameter>("shapes", "array of shapes, one for each bracnh", Validate.REQUIRED);
	
	RealParameter shapesParameter;
	double [][] categoryRates;
	
	
	public MultiGammaSiteModel() {
		shapeParameterInput.setRule(Validate.FORBIDDEN);
	}
	
	@Override
	public void initAndValidate() {
        shapesParameter = shapesParameterInput.get();
        super.initAndValidate();
        
        hasPropInvariantCategory = false;
        categoryCount = gammaCategoryCount.get();
        if (categoryCount <= 1) {
            throw new RuntimeException("MultiGammaSiteModel: Invalid category count (" + categoryCount + "). Should be at least 2");
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
    	Tree tree = nodeX.getTree();		
		
        double propVariable = 1.0;
        int cat = 0;

        // 	deal with invariant sites
        if (hasPropInvariantCategory) {
        	// categoryRates may not be create here, so do not initialise
        	// When created (usin new double[][]) its values will be 0 anyway
            // categoryRates[0] = 0.0;
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
       			a = shapesParameter.getValue(node.getNr());
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
        for (int i = 0; i < rates.length; i++) {
            rates[i] = categoryRates[node.getNr()][i] * mu;
        }

        return rates;
    }

    @Override
    protected boolean requiresRecalculation() {
        // we only get here if something is dirty in its inputs, so always return true
        ratesKnown = false;
        return true;
    }

    
    //@Override
    public boolean hasNodeIndependentCategories() {
    	return false;
    }
}
