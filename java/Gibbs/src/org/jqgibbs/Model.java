package org.jqgibbs;

import java.util.Map;

public abstract class Model {
	protected Map<String, RandomVar<? extends Flattenable>> params;
	protected Map<String, Flattenable> hypers;
	protected int dims;

	public Model() {
		super();
	}
	
	public Model(Map<String, Flattenable> hypers, int dims) {
		this.hypers = hypers;
		this.dims = dims;
	}
	
	public abstract ChainLink getInitialLink();

}
