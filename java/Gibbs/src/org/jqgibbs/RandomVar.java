package org.jqgibbs;

import java.lang.reflect.InvocationTargetException;

import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistMC;

public class RandomVar<T extends Flattenable & Cloneable> {
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

	@SuppressWarnings("unchecked")
	public RandomVar<T> cloneShallow() {
		T cloneNumValue = null;
		try {
			cloneNumValue = (T) this.numericValue.getClass()
					.getMethod("clone").invoke(this.numericValue);
		} catch (IllegalAccessException e) {
			assert false;
		} catch (InvocationTargetException e) {
			assert false;
		} catch (NoSuchMethodException e) {
			assert false;
		}
		return new RandomVar<T>(this.name, this.posterior, cloneNumValue);
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

	public String toString() {
		return this.numericValue.toString();
	}

}
