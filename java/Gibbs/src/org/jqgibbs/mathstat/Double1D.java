package org.jqgibbs.mathstat;

import java.util.Arrays;
import java.util.Collection;

import org.jqgibbs.Flattenable;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;

public class Double1D implements Flattenable,Cloneable {

	private DoubleMatrix1D dm;
	
	public Double1D(Double0D[] d0Ds) {
		this.dm = new DenseDoubleMatrix1D(new double[0]);
		this.addAll(Arrays.asList(d0Ds));
	}
	
	public Double1D(double... ds) {
		this.dm = new DenseDoubleMatrix1D(ds);
	}

	public Double1D(DoubleMatrix1D dm) {
		this.dm = dm;
	}

	public int length1D() {
		return this.size();
	}

	public Double1D rowVec() {
		return this;
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof Double1D)) {
			return false;
		}
		return this.dm.equals(((Double1D) o).dm);
	}
	
	@Override
	public int hashCode() {
		return this.dm.hashCode();
	}	
	
	@Override
	public String toString() {
		String s = "";
		String prefix = "";
		for (int i=0; i < this.size(); i++) {
			s = s + prefix + String.valueOf(this.dm.get(i));
			prefix = " ";
		}
		return s;
	}	
	
	@Override
	public Object clone() {
		return new Double1D(this.value());
	}
	
	public DoubleMatrix1D getDm() {
		return this.dm;
	}

	// FIXME
	//@Override
	public boolean addAll(Collection<? extends Double0D> c) {
		double[] ds = Arrays.copyOf(this.dm.toArray(),
				this.size()+c.size());
		int i = this.size();
		for (Double0D d : c) {
			if (d == null) {
				throw new NullPointerException();
			}
			ds[i] = d.value();
			i = i + 1;
		}
		this.dm = new DenseDoubleMatrix1D(ds);
		return true;
	}

	public int size() {
		return this.dm.size();
	}

	public double[] value() {
		return this.dm.toArray();
	}
	
	// FIXME - to copy or not???
	public DoubleMatrix1D toColt() {
		return this.dm.copy();
	}

	public Double0D get(int i) {
		return new Double0D(this.dm.get(i));
	}
	
	public Double2D toDouble2D(int h) {
		return new Double2D(AlgebraStatic.unvec(this.dm, h));
	}

	// Convenience methods
	
	public Double1D plus(Double1D o) {
		return new Double1D(AlgebraStatic.plusStatic(this.dm, o.dm));		
	}

	public Double2D outer(Double1D o) {
		DoubleMatrix2D r = AlgebraStatic.multOuterStatic(this.dm, o.dm);
		return new Double2D(r);
	}
	
	public Double1D mult(double k) {
		return new Double1D(AlgebraStatic.multStatic(k, this.dm));
	}

	public Double1D minus(Double1D v) {
		return new Double1D(AlgebraStatic.minusStatic(this.dm, v.dm));
	}

	public Double1D mult(Double2D other) {
		return new Double1D(AlgebraStatic.multStatic(this.dm, other
				.getDm()));
	}

	public double mult(Double1D other) {
		return AlgebraStatic.multStatic(this.dm, other.dm);
	}
	
	public Double1D minus(double k) {
		return new Double1D(AlgebraStatic.minusStatic(this.dm, k));
	}
	
	public Double1D exp() {
		double[] ed = Arrays.copyOf(this.value(), this.size());
		for (int i=0; i<this.size(); i++) {
			if (Double.isInfinite(ed[i])) {
				ed[i] = -Double.MAX_VALUE;
			}
			ed[i] = Math.exp(ed[i]);
		}
		return new Double1D(ed);
	}

	public Double1D plus(double k) {
		return new Double1D(AlgebraStatic.plusStatic(this.dm, k));
	}

	public Double0D sum() {
		double sum = 0;
		for (int i=0; i<this.size(); i++) {
			sum = sum + this.value()[i];
		}
		return new Double0D(sum);
	}

	public Double1D mult(Double0D k) {
		return this.mult(k.value());
	}

	public Double1D divide(Double0D k) {
		return this.divide(k.value());
	}

	public Double1D divide(double k) {
		return new Double1D(AlgebraStatic.divideStatic(k, this.dm));
	}

	public Double2D outer(Integer1D o) {
		double[] od = new double[o.size()];
		int[] ov = o.value();
		for (int i=0; i<ov.length; i++) {
			od[i] = ov[i];
		}
		DoubleMatrix2D r = AlgebraStatic.multOuterStatic(this.dm,
														 new DenseDoubleMatrix1D(od));
		return new Double2D(r);
	}

	public Double1D sqrt() {
		double[] v = this.value();
		for (int i=0; i<this.size(); i++) {
			v[i] = Math.sqrt(v[i]);
		}
		return new Double1D(v);
	}
}
