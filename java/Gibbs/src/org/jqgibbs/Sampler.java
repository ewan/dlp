package org.jqgibbs;

import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;

public abstract class Sampler {
	public static double MIN_VALUE = Double.MIN_VALUE;
	
	protected Model model;
	private Double2D data;
	protected ChainLink current;
	
	public Sampler(Model m, Double2D d) throws ProbDistParmException {
		this.model = m;
		this.setData(d);
		this.setChainLink(m.getInitialLink());
	}
	
	protected void setChainLink(ChainLink initialLink) {
		this.current = initialLink;
	}

	public Model getModel() {
		return this.model;
	}
	
	protected void setData(Double2D d) {
		this.data = d;
	}

	protected Double2D getData() {
		return data;
	}

	public abstract ChainLink variate();

	public ChainLink variateFast() {
		return this.variate();
	}
}
