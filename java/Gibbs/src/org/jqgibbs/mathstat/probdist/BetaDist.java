package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.random.Beta;

public class BetaDist extends ProbDist<Double0D> {
	private Double0D shape1;
	private Double0D shape2;

	private Beta betaGen;

	public BetaDist(Double0D shape1, Double0D shape2, boolean checkParms)
			throws ProbDistParmException {
		this.shape1 = shape1;
		this.shape2 = shape2;
		if (checkParms) {
			checkParms();
		}
		setUpFromParms();
	}

	/**
	 * @error Need to check to see if the distribution is initialized on a
	 *        0-args call
	 */
	protected void checkInitialized(Numeric... parms) {
		if (parms.length == 0)
			return;
		this.shape1 = (Double0D) parms[0];
		this.shape2 = (Double0D) parms[1];
		setUpFromParms();
	}


	/**
	 * Beta distribution.
	 * <br>
	 * <pre>
	 * f(x) \propto x^(shape1 - 1) * (1-x)^(shape2 - 1)
	 * </pre>
	 * 
	 * @param shape1 The "x" parameter (see above). Must be positive.
	 * @param shape2 The "1-x" parameter (see above). Must be positive.
	 * @throws ProbDistParmException if either parameter is non-positive
	 */	
	public BetaDist(Double0D shape1, Double0D shape2)
			throws ProbDistParmException {
		this(shape1, shape2, CHECK_PARMS);
	}

	private void checkParms() throws ProbDistParmException {
		if (this.shape1.value() <= 0) {
			throw new ProbDistParmException("shape1 must be positive");
		}
		if (this.shape2.value() <= 0) {
			throw new ProbDistParmException("shape2 must be positive");
		}
	}

	protected Beta getBetaGen() {
		return this.betaGen;
	}

	private Double0D getShape1() {
		return shape1;
	}

	private Double0D getShape2() {
		return shape2;
	}

	private void setUpFromParms() {
		if (this.getBetaGen() == null) {
			this.betaGen = new Beta(this.getShape1().value(), this.getShape2()
					.value(), RandomEngineSelector.getEngine());
		} else {
			this.getBetaGen().setState(this.getShape1().value(),
					this.getShape2().value());
		}
	}

	@Override
	protected Double0D genVariate() throws ProbDistParmException {
		assert this.initialized;
		return new Double0D(this.getBetaGen().nextDouble());
	}

	@Override
	protected double getDensity(Double0D pt) {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}
}
