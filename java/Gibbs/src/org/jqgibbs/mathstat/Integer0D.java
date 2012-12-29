package org.jqgibbs.mathstat;

import org.jqgibbs.Flattenable;

public class Integer0D implements Flattenable {

	private int i;
	
	public Integer0D(int i) {
		this.i = i;
	}
	
	public int value() {
		return this.i;
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof Integer0D)) {
			return false;
		}
		return this.value() == ((Integer0D) o).value();
	}

	@Override
	public int hashCode() {
		return this.value();
	}
	
	public Double1D rowVec() {
		return new Double1D((double) this.value()); 
	}
	
	public int length1D() {
		return 1;
	}
	
	// Convenience methods
	
	public Integer0D plus(int k) {
		return new Integer0D(this.value() + k);
	}

	public Double0D pow(double a) {
		return new Double0D(Math.pow(this.value(), a));
	}

	public Integer0D plus(Integer0D k) {
		return this.plus(k.value());
	}

	public Double0D divide(Double0D k) {
		return (new Double0D(this.value()).divide(k));
	}
	

}
