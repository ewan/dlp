package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;

/**
 * 
 * Note: creating a NormalDist object advances the pseudo-random number sequence
 * of RandomEngine.getSequence() by one
 * 
 * @author ewan
 * @author josh
 */
public class NormalDist extends ProbDist<Double0D> {

	private RandomEngine randomEngine;
	private Normal normalGen;

	private Double0D mu;
	private Double0D sg;

	/**
	 * @review note that this function (which is already supposed to be
	 *         reworked/removed anyway) never uses the checkParms flag
	 */
	protected void checkInitialized(Numeric... parms) {
		if (parms.length < 2) {
			if (!this.initialized) {
				throw new IllegalStateException(
						"use of uninitialized probability distribution");
			}
		} else {
			this.setParms((Double0D) parms[0], (Double0D) parms[1]);
		}
	}

	public NormalDist(Double0D mu, Double0D sg, boolean checkParms) {
		this.setParms(mu, sg, checkParms);
	}

	public NormalDist(Double0D mu, Double0D sg) {
		this(mu, sg, CHECK_PARMS);
	}

	public void setParms(Double0D mu, Double0D sg, boolean checkParms) {
		this.mu = mu;
		this.sg = sg;
		this.setUpFromParms(checkParms);
		this.initialized = true;
	}

	public void setParms(Double0D mu, Double0D sg) {
		this.setParms(mu, sg, CHECK_PARMS);
	}

	private void checkParms() {
		if (this.sg.value() < 0) {
			throw new IllegalArgumentException(
					"Expected non-negative value for sg");
		}
	}

	private void setUpFromParms(boolean checkParms) {
		if (checkParms) {
			this.checkParms();
		}
		if (this.normalGen == null) {
			this.randomEngine = RandomEngineSelector.getEngine();
			this.normalGen = new Normal(this.mu.value(), this.sg.value(),
					this.randomEngine);
		} else {
			assert this.randomEngine != null;
			this.normalGen.setState(this.mu.value(), this.sg.value());
		}
	}

	@Override
	protected Double0D genVariate() {
		return new Double0D(this.normalGen.nextDouble());
	}

	@Override
	protected double getDensity(Double0D pt) {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}
}
