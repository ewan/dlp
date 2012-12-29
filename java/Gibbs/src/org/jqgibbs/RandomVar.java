package org.jqgibbs;

import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistMC;

public class RandomVar<T extends Flattenable> {
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

	public RandomVar<T> cloneShallow() {
		return new RandomVar<T>(this.name, this.posterior, this.numericValue);
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
