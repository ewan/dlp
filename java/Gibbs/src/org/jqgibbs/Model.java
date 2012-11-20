package org.jqgibbs;

import java.util.List;
import java.util.Map;

import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Integer1D;
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.mathstat.RandomVar;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;

public abstract class Model {
	protected Map<String, RandomVar<? extends Numeric<?>>> params;
	protected Map<String, Numeric<? extends Numeric<?>>> hypers;
	protected int dims;

	public Model() {
		super();
	}
	
	public Model(Map<String, Numeric<? extends Numeric<?>>> hypers, int dims) throws ProbDistParmException {
		this.hypers = hypers;
		this.dims = dims;
	}
	
	public RandomVar<? extends Numeric<?>> getParam(String paramName) {
		return this.params.get(paramName);
	}
	
	public RandomVar<? extends Numeric<?>> getParam(String paramName, Double1D v) {
		RandomVar<? extends Numeric<?>> rv = this.params.get(paramName); // FIXME - check!
		return rv.cloneFromVector(v);
	}

	public Numeric<? extends Numeric<?>> getHyper(String hyperName) {
		return this.hypers.get(hyperName);
	}

	public ChainLink getChainLink(Double1D v, Map<String, Integer1D> paramIndices) throws GibbsException {
		ChainLink cl = new ChainLink();
		int[] indices;
		Double1D values;
		for (String paramName : paramIndices.keySet()) {
			indices = paramIndices.get(paramName).value();
			values = v.getAll(indices);
			cl.add(this.getParam(paramName, values));
		}
		return cl;
	}
	
	public abstract ChainLink getInitialLink();

	public List<String> getVars() {
		return this.getInitialLink().getOrder();
	}
}
