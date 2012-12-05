package org.jqgibbs.mathstat.probdist;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.jqgibbs.mathstat.Numeric;

/**
 * Generate variates and obtain density values for events under a given
 * probability distribution.
 * 
 * <br>
 * When a variate is sampled using variate(), (as opposed to variateFast(), the
 * other public sampling method, which may or may not go through the same
 * initialization sequence), a sequence of variates sampled using variatesIID(),
 * or the density of a point obtained using density(), the distribution will
 * invoke the abstract checkInitialized() method before performing the
 * operation. Extending classes should use the checkInitialized() method to do
 * two different things:
 * 
 * <ul>
 * <li>Use the arguments passed the invoking method (variate(), etc) and passed
 * on to checkInitialized() to set the distribution's parameters; any value
 * checks should also be performed here</li>
 * <li>Ensure that zero-args calls to variate(), etc, are only made after the
 * distribution's parameters have been set, either by a previous call to
 * checkInitialized(), or at class initialization</li>
 * </ul>
 * 
 * @design The value of using the varargs is that this class can exist to
 *         specify the interface. If all the calls were made zero-args,
 *         checkInitialized would become trivial, we would avoid the varargs
 *         array and the accompanying need to check the number and type of the
 *         arguments, as well as the explicit casts. This is essentially what
 *         the COLT distributions do: they have a setState() method to change
 *         the parameters without destroying the object which is not specified
 *         as part of the interface. However, although this solution avoids the
 *         awkwardness of the varargs, it makes the use of the class clunky, as
 *         we can no longer make most of our calls to variate() without adding a
 *         separate line to update the parameter values (we use virtually no
 *         zero-args calls). The COLT solution is to add another variate(...)
 *         method by convention with the full list of arguments -- but not to do
 *         the same for cdf() and pdf(), presumably because these get used in
 *         this way less often. My usage is similar, and so I would be happy to
 *         go with the same design. Note that the COLT people decided not to
 *         have the extra variate(...) call change the internal state, so that
 *         it is essentially being used as just a static method. In fact, by
 *         another convention, there is also a method (which we would call)
 *         staticVariate(...), and I would assume variate(...) is probably just
 *         a wrapper for this method). It's useful to have both for the same
 *         reason it's useful that we didn't make all the methods here flat-out
 *         static: so that we can do things like saving an object as a ProbDist
 *         and then drawing a variate from it without knowing the extending
 *         class -- which is the whole basis for the design of the sampler.
 *         However, this would unfortunately not kick in in any useful way for
 *         us, because the only case in which we can make use of the ProbDist
 *         interface is unavailable in the sampler situation, where the
 *         parameters keep changing with the state of the Markov chain. The best
 *         solution seems to be to make the above changes, and then subclass
 *         ProbDist with a class that will be the parent class for all the model
 *         inner classes (currently, ProbDistInitializeByChain). The difference
 *         would be that we could still drop the varargs and simply add a
 *         variate(ChainLink) method (we probably ought to make the data part of
 *         the initialization; it's just another observed variable, and not all
 *         variables will ultimately make use of it).
 * 
 * @review The old design had a varargs constructor as part of the interface, as
 *         well. This constructor was not really to be messed with, and it
 *         worked by first running the code that set up all the type and value
 *         checks for the varargs, (standardized as enforced by checkParm), and
 *         then doing some pre-initialization which was later to be presupposed
 *         by initializeParms. In practice, what this meant was something for
 *         ProbDistInitializeByChain (which preInitialized by setting the fixed
 *         parameters) and nothing for ProbDistInitializeDirectly (which just
 *         called initializeParms). Thus, in effect, the varargs in the constructor
 *         actually referred to two totally different things - in one case,
 *         something totally different than the varargs
 *         in the ordinary methods, namely the fixed parameters used as part
 *         of a ProbDistInitializeByChain; in the other, the regular parameters
 *         with their regular checks. We can safely remove this constructor
 *         and the pre-initialization it does and in so doing just remove all
 *         vestiges of any parameter initialization from the constructor, and that's
 *         what's been done here. This means that we now have available for
 *         actual use a zero-args constructor (namely, the one inherited from
 *         Object) whereas before the parameter initialization could in
 *         principle have barfed on a zero-args call (because it would have
 *         gone to the varargs). Perhaps somewhat confusingly, this was the
 *         case (barfing, that is) for ProbDistInitializeByChain but not
 *         ProbDistInitializeDirectly. One possible explanation for this is
 *         that the constructor was simply the only method I set up for setting
 *         the fixed parameters, as I was already putting the variate(...)
 *         args to a dedicated use passing in the data/chain. At any rate,
 *         I believe there is now a use for the zero-arg constructor, because
 *         there is a reason to create an uninitialized object: it allows
 *         us to clean up the branching logic which is used throughout the
 *         model inner classes whereby one either creates, or simply resets
 *         the parameters of, a distribution object (such as the
 *         MVNormalDist mentioned in a review comment in PostMDist). Oddly,
 *         this is not a case where I actually *needed* to do this, as far
 *         as I can tell, as a zero-args constructor call should have had the
 *         desired effect for a ProbDistInitializeDirectly all along. This
 *         should be changed with a very careful eye to things that might
 *         be breaking (perhaps try changing first in the ewan branch, which
 *         preserves the varargs constructor).
 * 
 * @author josh
 * @author ewan
 * 
 * @param <T>
 *            Type of atomic event for this distribution.
 */
public abstract class ProbDist<T extends Numeric> {
	/**
	 * Toggles value checks on parameters during initialization. It is up to the
	 * extending class to ensure that the value of CHECK_PARMS is respected.
	 */
	protected static final boolean CHECK_PARMS = true;
	protected boolean initialized;

	protected abstract double getDensity(T pt) throws ProbDistParmException;

	protected double getLogDensity(T pt) throws ProbDistParmException {
		return Math.log(this.getDensity(pt));
	}

	protected abstract T genVariate() throws ProbDistParmException;

	/*
	 * Initializes chain parms for ProbDistInitializeByChain, Fixed parms for
	 * other distributions
	 */
	protected abstract void checkInitialized(Numeric... parms);

	/*
	 * This is declared as final in order to guarantee that
	 * ProbDistParmException is never thrown on a zero-args call.
	 */
	public final double density(T pt, Numeric... parms)
			throws ProbDistParmException {
		this.checkInitialized(parms);
		return this.getDensity(pt);
	}

	public final double logDensity(T pt) throws ProbDistParmException {
		return this.getLogDensity(pt);
	}

	/*
	 * This is declared as final in order to guarantee that
	 * ProbDistParmException is never thrown on a zero-args call.
	 */
	public final T variate(Numeric... parms) throws ProbDistParmException {
		this.checkInitialized(parms);
		return this.genVariate();
	}

	/**
	 * @design Do we need this class?
	 * 
	 * @param n
	 *            Number of variates to return.
	 * @param parms
	 *            Varargs determining the distribution parameters.
	 * @return A sequence of variates from a distribution with parameters
	 *         determined by parms.
	 * @throws ProbDistParmException
	 *             if there are value errors or other problems with parms
	 */
	/*
	 * This is declared as final in order to guarantee that
	 * ProbDistParmException is never thrown on a zero-args call.
	 */
	@SuppressWarnings("unchecked")
	public final List<T> variatesIID(int n, Numeric... parms)
			throws ProbDistParmException {
		if (n < 1) {
			return new ArrayList<T>();
		}
		this.checkInitialized(parms);
		List<T> l = new LinkedList<T>();
		for (int i = 0; i < n; i++) {
			l.add(this.genVariate());
		}
		List<T> s = new ArrayList<T>();
		for (Numeric elem : l.get(0).rowVec()) {
			s.add((T) elem);
		}
		if (n > 1) {
			l.remove(0);
			s.addAll(l);
		}
		return s;
	}

	public T variateFast(Numeric... parms) throws ProbDistParmException {
		return this.variate(parms);
	}

}
