package org.jqgibbs;

import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

public class RandomEngineSelector {
	private static final boolean fixedSeed = true;
	private static RandomEngine engine = null;
	
	public static RandomEngine getEngine() {
		if(engine == null) {
			if(fixedSeed) {
				engine = new MersenneTwister();
			} else {
				engine = RandomEngineSelector.getEngine();
			}
		}
		return engine;
	}
}