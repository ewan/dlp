package org.jqgibbs;

import umontreal.iro.lecuyer.rng.RandomStream;
import cern.jet.random.engine.RandomEngine;
import cern.jet.random.engine.MersenneTwister;

public class RandomEngineSelector {
	private static int seed = 100;
	private static RandomEngine re;
	private static RandomStream rs;
	
	public static RandomEngine getEngine() {
		if (re == null) {
			//re = new MersenneTwister(seed);
			re = new MersenneTwister();
		}
		return re;
	}
	
	
	public static RandomStream getStream() {
		if (rs == null) {
			rs = new RandomEngineStream(RandomEngineSelector.getEngine());
		}
		return rs;
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
