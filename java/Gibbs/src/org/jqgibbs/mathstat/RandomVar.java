package org.jqgibbs.mathstat;

import org.jqgibbs.ChainLink;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistMC;

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

public class RandomVar<T extends Numeric> {
	private String name;
	private ProbDistMC<T> posterior;
	private T numericValue;

	public RandomVar(String name, ProbDistMC<T> posterior, T t) {
		this.name = name;
		this.posterior = posterior;
		this.numericValue = t;
	}

	public void updatePosterior(ChainLink l) {
		this.numericValue = this.posterior.variate(l);
	}

	public RandomVar<T> samplePosterior(ChainLink l) {
		T t = this.posterior.variate(l);
		return this.cloneWith(t);
	}

	public RandomVar<T> cloneWith(T t) {
		return new RandomVar<T>(this.name, this.posterior, t);
	}

	public String getName() {
		return this.name;
	}
	
	public T getNumericValue() {
		return numericValue;
	}

	public ProbDist<T> getPosterior() {
		return posterior;
	}

	@SuppressWarnings("unchecked")
	public Object clone() throws CloneNotSupportedException {
		RandomVar<T> clone = new RandomVar<T>(this.name,
				this.posterior, (T) this.numericValue.clone());
		return clone;
	}

	public RandomVar<T> cloneFromVector(Double1D v) {
		T clone = (T) this.numericValue.cloneFromVector(v);
		return this.cloneWith(clone);
	}

	// @Override
	public String toString() {
		return this.numericValue.toString();
		// return this.getName() + ": " + this.getNumericValue().toString();
	}

	@Override
	public boolean equals(Object o) {
		if (!(this.getClass().isAssignableFrom(o.getClass()))) {
			return false;
		}
		if (!(this.name.equals(((RandomVar<T>) o).name))) {
			return false;
		}
		return this.numericValue.equals(((RandomVar<T>) o).numericValue);
	}

	@Override
	public int hashCode() {
		return this.name.hashCode() + this.numericValue.hashCode();
	}

	public Double1D rowVec() {
		return this.numericValue.rowVec();
	}

	public int length1D() {
		return this.numericValue.length1D();
	}
}
