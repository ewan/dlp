package org.jqgibbs.mathstat;

import org.jqgibbs.ChainLink;
import org.jqgibbs.SamplerData;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;

/*
 * It would be *really* nice if we could implement value() such that it returns
 * this.getNumericValue().value(). But sadly, we can't, because Numeric.value() returns
 * an unknown primitive type - so Numeric.value() per se doesn't exist (isn't guaranteed
 * by the contract for Numeric). Similarly, there's no way to state its return type here,
 * either.
 * 
 * That way, we could treat RandomVar as an ordinary Numeric, and not have to worry about
 * the distinction in places where we don't care.
 * 
 *  							   >>>AAARRGGGHH!!<<<
 * 
 * There are two alternatives. One would be to use autoboxing/unboxing, but this would be
 * an unnecessary speed hit (we are using COLT matrices for speed, after all), and wouldn't help
 * with the problem of primitive array values, which means it's out before it hits the gate.
 * The other would be to declare this class as abstract and get subclasses to implement value
 * by convention. It's not clear to me that I gain too much by doing that, but it might be
 * useful/necessary in the future. 
 *   
 */

public class RandomVar<T extends Numeric> implements Numeric {
	private String name;
	private ProbDist<T> prior;
	private ProbDist<T> posterior;
	private T numericValue;

	public RandomVar(String name, ProbDist<T> prior,
			ProbDist<T> posterior, T t) {
		this.setName(name);
		this.setPrior(prior);
		this.setPosterior(posterior);
		this.setNumericValue(t);
	}

	public void updatePrior() {
		try {
			this.setNumericValue(this.getPrior().variate());
		} catch (ProbDistParmException e) {
			throw new IllegalStateException(
					"Encountered bug: Unexpected ProbDistParmException "
							+ "from 0-arg call to ProbDist.variate()", e);
		}
	}

	public void updatePosterior(ChainLink l, Double2D d)
			throws ProbDistParmException {
		this.setNumericValue(this.getPosterior().variate(l, d));
	}	
	
	public RandomVar<T> samplePrior() {
		try {
			T t = this.getPrior().variate();
			return this.cloneWith(t);
		} catch (ProbDistParmException e) {
			throw new IllegalStateException(
					"Encountered bug: Unexpected ProbDistParmException "
							+ "from 0-arg call to ProbDist.variate()", e);
		}
	}

	public RandomVar<T> samplePosterior(ChainLink l, SamplerData d)
			throws ProbDistParmException {
		T t = this.getPosterior().variate(l, d);
		return this.cloneWith(t);
	}

	public RandomVar<T> cloneWith(T t) {
		return new RandomVar<T>(this.getName(), this.getPrior(), this
				.getPosterior(), t);
	}

	public String getName() {
		return this.name;
	}

	protected void setNumericValue(T numericValue) {
		this.numericValue = numericValue;
	}

	public T getNumericValue() {
		return numericValue;
	}

	protected void setName(String name) {
		this.name = name;
	}

	protected void setPrior(ProbDist<T> prior) {
		this.prior = prior;
	}

	public ProbDist<T> getPrior() {
		return prior;
	}

	protected void setPosterior(ProbDist<T> posterior) {
		this.posterior = posterior;
	}

	public ProbDist<T> getPosterior() {
		return posterior;
	}

	@SuppressWarnings("unchecked")
	public Object clone() throws CloneNotSupportedException {
		RandomVar<T> clone = new RandomVar<T>(this.getName(), this.getPrior(),
				this.getPosterior(), (T) this.getNumericValue().clone());
		return clone;
	}

	public RandomVar<T> cloneFromVector(Double1D v) {
		T clone = (T) this.getNumericValue().cloneFromVector(v);
		return this.cloneWith(clone);
	}
	
	//@Override
	public String toString() {
		return this.getNumericValue().toString();
//		return this.getName() + ": " + this.getNumericValue().toString();
	}
	
	@SuppressWarnings("unchecked")
	//@Override
	public boolean equals(Object o) {
		if (!(this.getClass().isAssignableFrom(o.getClass()))) {
			return false;
		}
		if (!(this.getName().equals(((RandomVar<T>) o).getName()))) {
			// FIXME - how the hell do I do this?
			return false;
		}
		return this.getNumericValue().equals(((RandomVar<T>) o).getNumericValue());
	}
	
	//@Override
	public int hashCode() {
		return this.getName().hashCode() + this.getNumericValue().hashCode();
	}
	
	//@Override
	public Double1D rowVec() {
		return this.getNumericValue().rowVec();
	}

	//@Override
	public int length1D() {
		return this.getNumericValue().length1D();
	}

	public void updatePosteriorFast(ChainLink l, Double2D d) throws ProbDistParmException {
		this.setNumericValue(this.getPosterior().variateFast(l, d));		
	}
}
