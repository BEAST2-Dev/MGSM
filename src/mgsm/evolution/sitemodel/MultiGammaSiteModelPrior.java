package mgsm.evolution.sitemodel;

import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

@Description("Prior on gammashape parameters for Multi-Gamma Site Model to ensure branches from root "
		+ "have identical gamma shape (this is required in order to prevent identifiability issues)")
public class MultiGammaSiteModelPrior extends Distribution {
	public Input<Tree> treeInput = new Input<Tree>("tree", "tree associated with gamma shape site model", Validate.REQUIRED);
	public Input<Parameter> shapesParameterInput = new Input<Parameter>("shapes", "array of shapes, one for each bracnh", Validate.REQUIRED);
	//public Input<IntegerParameter> indicatorParameterInput = new Input<IntegerParameter>("indicator", "array of indicator, one for each bracnh", Validate.XOR, shapesParameterInput);

	Tree tree;
	Parameter<?> shapesParameter;
	
	@Override
	public void initAndValidate() {
		shapesParameter = shapesParameterInput.get();
		tree = treeInput.get();
    	if (shapesParameter.getDimension() != tree.getNodeCount()) {
            if (shapesParameter instanceof RealParameter) {
            	((RealParameter)shapesParameter).setDimension(tree.getNodeCount());
            } else {
            	((IntegerParameter)shapesParameter).setDimension(tree.getNodeCount());
            }
    	}

	}
	
	
    public double calculateLogP() {
        logP = 0;
        Node root = tree.getRoot();
        final double x;
        if (shapesParameter instanceof RealParameter) {
	        x = ((RealParameter)shapesParameter).getValue(root.getLeft().getNr()) - ((RealParameter)shapesParameter).getValue(root.getRight().getNr());
        } else {
	        x = ((IntegerParameter)shapesParameter).getValue(root.getLeft().getNr()) - ((IntegerParameter)shapesParameter).getValue(root.getRight().getNr());
        }
        final double standardDeviation = 0.01;
        logP = -(x) * (x) / (2.0 * standardDeviation * standardDeviation);
        return logP;
    }

	
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub

	}

}
