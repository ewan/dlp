package org.jqgibbs;

import java.util.Map;

public abstract class Model {
	protected Map<String, RandomVar<? extends Flattenable>> params;
	protected Map<String, Flattenable> hypers;
	protected Map<String, Object> extended_hypers;

	public Model() {
		super();
	}
	
	public Model(Map<String, Flattenable> hypers) {
		this.hypers = hypers;
	}
	
	public abstract ChainLink getInitialLink();

}
