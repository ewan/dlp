package org.jqgibbs;

public abstract class Sampler {
	public static double MIN_VALUE = Double.MIN_VALUE;
	
	protected Model model;
	protected ChainLink current;
	
	public Sampler(Model m) {
		this.model = m;
		this.setChainLink(m.getInitialLink());
	}
	
	protected void setChainLink(ChainLink initialLink) {
		this.current = initialLink;
	}

	public Model getModel() {
		return this.model;
	}
	
	public abstract ChainLink variate();

	public ChainLink variateFast() {
		return this.variate();
	}
}
