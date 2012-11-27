package org.jqgibbs.mathstat;

import java.math.BigDecimal;
import java.math.RoundingMode;

public class Double0D implements Numeric {
	private static final int SCALE = 35;
	private BigDecimal d;
	private boolean dirty;
	private double cached;
	
	public Double0D(double d) {
		this.d = new BigDecimal(d);
		this.cached = d;
		this.dirty = false;
	}
	
	protected Double0D(BigDecimal d) {
		this.d = d;
		this.cached = this.d.doubleValue();
		this.dirty = false;
	}
	
	public double value() {
		if (this.dirty) {
			this.cached = this.d.doubleValue(); 
		}
		return this.cached;
	}
	
	private BigDecimal getBigDecimal() {
		return d;
	}

	
	public Double1D sequence() {
		return new Double1D(this.value());
	}

	
	public Object clone() throws CloneNotSupportedException {
		return new Double0D(this.value());
	}

	
	public Double0D cloneFromVector(Double1D v) {
		if (v.size() == 0) {
			return new Double0D(0); // FIXME
		}
		return new Double0D(v.get(0).value());
	}

	
	public boolean equals(Object o) {
		if (!(o instanceof Double0D)) {
			return false;
		}
		return this.value() == ((Double0D) o).value();
	}

	
	public int hashCode() {
		return (int) this.value();
	}

	public Double0D plus(double d) {
		return new Double0D(this.getBigDecimal().add(new BigDecimal(d))); // FIXME
	}
	
	
	public String toString() {
		return String.valueOf(this.value());
	}

	public Double0D mult(double e) {
		return new Double0D(this.value()*e); // FIXME
	}

	public Double0D plus(Integer0D k) {
		return new Double0D(this.getBigDecimal().add(new BigDecimal(k.value())));
	}

	public Double0D divide(Double0D k) {
		return new Double0D(this.getBigDecimal().divide(k.getBigDecimal(), Double0D.SCALE, RoundingMode.HALF_EVEN));
	}

	public Double0D minus(Double0D k) {
		return new Double0D(this.getBigDecimal().subtract(k.getBigDecimal()));
	}

	public Double0D mult(Double0D k) {
		return new Double0D(this.getBigDecimal().multiply(k.getBigDecimal()));
	}

	public Double0D plus(Double0D k) {
		return new Double0D(this.getBigDecimal().add(k.getBigDecimal()));
	}
	
	public Double0D pow(int a) {
		return new Double0D(this.getBigDecimal().pow(a));// FIXME
	}

	public Double0D recip() {
		return new Double0D((new BigDecimal(1)).divide(this.getBigDecimal(), Double0D.SCALE, RoundingMode.HALF_EVEN));
	}

	
	public Double1D rowVec() {
		return new Double1D(this.value());
	}
	
	
	public int length1D() {
		return 1;
	}

	public Double0D sqrt() {
		return new Double0D(Math.sqrt(this.value()));
	}

	public Double0D exp() {
		return new Double0D(Math.exp(this.value()));
	}

	public void set(double d) {
		this.d = new BigDecimal(d);
		this.cached = d;
		this.dirty = false;
	}
}
