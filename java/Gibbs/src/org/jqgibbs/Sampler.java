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

	public abstract ChainLink variate();
}
