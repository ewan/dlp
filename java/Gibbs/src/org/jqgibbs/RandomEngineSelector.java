package org.jqgibbs;

import java.util.Date;

import umontreal.iro.lecuyer.rng.RandomStream;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

public class RandomEngineSelector {
	private static boolean fixedSeed = false;
	
	private static RandomEngine engine = null;
	private static RandomStream stream = null;
	
	public static void setFixedSeed(boolean fixedSeed) {
		if (fixedSeed != RandomEngineSelector.fixedSeed) {
			RandomEngineSelector.fixedSeed = fixedSeed;
			if (RandomEngineSelector.engine != null) {
				RandomEngineSelector.engine = null;
				RandomEngineSelector.stream = null;
			}
		}
	}

	public static RandomEngine getEngine() {
		if (RandomEngineSelector.engine == null) {
			if (RandomEngineSelector.fixedSeed) {
				RandomEngineSelector.engine = new MersenneTwister();
			} else {
				RandomEngineSelector.engine = new MersenneTwister(new Date());
			}
		}
		return RandomEngineSelector.engine;
	}

	public static RandomStream getStream() {
		if (RandomEngineSelector.stream == null) {
			RandomEngine e = RandomEngineSelector.getEngine();
			RandomEngineSelector.stream = new RandomEngineStream(e);
		}
		return RandomEngineSelector.stream;
	}

	protected static class RandomEngineStream implements RandomStream {
		private RandomEngine engine = null;

		public RandomEngineStream(RandomEngine engine) {
			this.engine = engine;
		}

		public void nextArrayOfDouble(double[] u, int start, int n) {
			boolean notPastEnd = true;
			for (int k=0; k<n; k++) {
				double rand = this.nextDouble();
				if (notPastEnd) {
					if (start+k < u.length) {
						u[start+k] = rand;
					} else {
						notPastEnd = false;
					}
				}
			}			
		}

		public void nextArrayOfInt(int i, int j, int[] u, int start, int n) {
			boolean notPastEnd = true;
			for (int k=0; k<n; k++) {
				int randInt = this.nextInt(i, j);
				if (notPastEnd) {
					if (start+k < u.length) {
						u[start+k] = randInt;
					} else {
						notPastEnd = false;
					}
				}
			}
		}

		public double nextDouble() {
			return this.engine.nextDouble();
		}

		public int nextInt(int i, int j) {
			double unif = this.engine.nextDouble();
			double unifMax = (j - i + 1) * unif;
			int unifMaxInt = (int) Math.floor(unifMax);
			return (i + unifMaxInt);
		}

		public void resetNextSubstream() {
			throw new UnsupportedOperationException();
		}

		public void resetStartStream() {
			throw new UnsupportedOperationException();
		}

		public void resetStartSubstream() {
			throw new UnsupportedOperationException();
		}

	}
}