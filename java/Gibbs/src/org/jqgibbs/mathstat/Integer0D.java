package org.jqgibbs.mathstat;

import java.math.BigDecimal;

public class Integer0D implements Numeric {

	private int i;
	
	public Integer0D(int i) {
		this.i = i;
	}
	
	public int value() {
		return this.i;
	}

	//@Override
	public Integer1D sequence() {
		return new Integer1D(this);
	}

	public int compareTo(Integer0D o) {
		return this.value() - o.value();
	}
	
	public int compareTo(int o) {
		return this.value() - o;
	}

	//@Override
	public Object clone() throws CloneNotSupportedException {
		return new Integer0D(i);
	}

	//@Override
	public Integer0D cloneFromVector(Double1D v) {
		if (v.size() == 0) {
			return new Integer0D(0); // FIXME
		}
		return new Integer0D((int) v.get(0).value());
	}

	//@Override
	public boolean equals(Object o) {
		if (!(o instanceof Integer0D)) {
			return false;
		}
		return this.value() == ((Integer0D) o).value();
	}

	//@Override
	public int hashCode() {
		return this.value();
	}

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
	
	//@Override
	public Double1D rowVec() {
		return new Double1D((double) this.value()); 
	}
	
	//@Override
	public int length1D() {
		return 1;
	}
}
