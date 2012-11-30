package org.jqgibbs;

import cern.jet.random.engine.RandomEngine;
import cern.jet.random.engine.MersenneTwister;

public class RandomEngineSelector {
	private static int seed = 100;
	private static RandomEngine re;
	
	public static RandomEngine getEngine() {
		if (re == null) {
			//re = new MersenneTwister(seed);
			re = new MersenneTwister();
		}
		return re;
	}
}
