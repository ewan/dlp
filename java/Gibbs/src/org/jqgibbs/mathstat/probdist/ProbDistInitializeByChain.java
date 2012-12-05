package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.ChainLink;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Numeric;

/**
 * 
 * A probability distribution whose parameters
 * depend on an object containing the
 * state of a Markov chain.
 * <br> 
 * The only arguments passed in to
 * checkInitialized() (and consequently to those
 * methods that invoke checkInitialized) are
 * the data and a chain link. All other information
 * is passed as arguments to the constructor.
 * <br>
 * The motivation for this division of labour
 * is that the chain parameters will be changing
 * with every sample in general, while the other information
 * (the "fixed parameters") only needs to be
 * set once. It would be clunky to have to pass
 * it in repeatedly. As a result, this supports a scheme
 * where highly abstracted code is used to do
 * the sampling, calling variate() without ever
 * needing access to any information other than
 * the previous state of the chain and the data.
 * 
 * @design The data is a fixed parameter!
 * @design These distributions will *definitely*
 * never take advantage of the zero-args call.
 * @review These distributions used to override
 * variateFast. Rather than simply
 * calling variate (i.e., doing nothing of
 * interest, as the parent method does) they used to call a faster version
 * of checkInitialized that would skip the value
 * checks. Now this mechanism has been replaced
 * with ProbDist.CHECK_PARMS (see docs at
 * ProbDist).
 * 
 * Initialization sequence in checkInitialized():
 * <ol>
 * <li>Check to see if the distribution is already initialized on a zero-args
 * call</li>
 * <li>Set the data (the second argument to checkInitialized)</li>
 * <li>Extract the information contained in the chain link (the first argument
 * to checkInitialized) by calling initializeChainParms()</li>
 * <li>Set the distribution parameters by calling initializeChainParms()</li>
 * </ol>
 * 
 * @author josh
 * @author ewan
 * @param <T>
 *            Type of atomic event for this distribution.
 */
public abstract class ProbDistInitializeByChain<T extends Numeric> extends
		ProbDist<T> {
	private Double2D data;

	protected Double2D getSamplerData() {
		return this.data;
	}

	protected abstract void initializeChainParms(ChainLink l);

	protected abstract void setUpFromChainParms();

	/**
	 * @error Need to check to see if the distribution is initialized on a
	 *         0-args call
	 */
	@Override
	protected void checkInitialized(Numeric... parms) {
		if (parms.length == 0)
			return;
		this.data = (Double2D) parms[1];
		this.initializeChainParms((ChainLink) parms[0]);
		this.setUpFromChainParms();
	}
}

