package org.jqgibbs.mathstat;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;

public class Double1D implements Numeric, Iterable<Double0D> {

	private DoubleMatrix1D dm;
	
	public Double1D(Double0D[] d0Ds) {
		this.setDm(new double[0]);
		this.addAll(Arrays.asList(d0Ds));
	}
	
	public Double1D(Double[] doubles) {
		double[] ds = new double[doubles.length];
		for (int i=0; i<doubles.length; i++) {
			if (doubles[i] == null) {
				throw new NullPointerException();
			}
			ds[i] = doubles[i].doubleValue();
		}
		this.setDm(ds);
	}
	
	public Double1D(double... ds) {
		this.setDm(ds);
	}

	public Double1D(DoubleMatrix1D dm) {
		this.setDm(dm);
	}

	private synchronized void setDm(DoubleMatrix1D dm) {
		this.dm = dm;
	}
	
	private void setDm(double[] ds) {
		this.setDm(new DenseDoubleMatrix1D(ds));
	}

	synchronized DoubleMatrix1D getDm() {
		return this.dm;
	}

	//@Override
	public Double2D sequence() {
		return new Double2D(this);
	}

	//@Override
	public Double0D set(int i, Double0D t) throws IndexOutOfBoundsException {
		if (i > this.size() || i < 0) {
			throw new IndexOutOfBoundsException("Tried to set index out of range"); // FIXME
		}
		if (i == this.size()) {
			this.add(t);
		} else {
			this.getDm().set(i, t.value());
		}
		return t;
	}	
	
	//@Override
	public boolean add(Double0D d0D) {
		double[] ds = Arrays.copyOf(this.getDm().toArray(), this.size()+1);
		if (d0D == null) {
			throw new NullPointerException();
		}
		ds[this.size()] = d0D.value();
		this.setDm(ds);
		return true;
	}

	//@Override
	public boolean addAll(Collection<? extends Double0D> c) {
		double[] ds = Arrays.copyOf(this.getDm().toArray(),
				this.size()+c.size());
		int i = this.size();
		for (Double0D d : c) {
			if (d == null) {
				throw new NullPointerException();
			}
			ds[i] = d.value();
			i = i + 1;
		}
		this.setDm(ds);
		return true;
	}

	//@Override
	public boolean contains(Object o) {
		double d;
		if (o == null) {
			throw new NullPointerException();
		} else if (o instanceof Double0D) {
			d = ((Double0D) o).value();
		} else if (o instanceof Double) {
			d = ((Double) o).doubleValue();
		} else {
			return false;
		}
		for (int i=0;i<this.size();i++) {
			if (this.get(i).value() == d) {
				return true;
			}
		}
		return false;
	}

	//@Override
	public int size() {
		return this.getDm().size();
	}

	public double[] value() {
		return this.getDm().toArray();
	}
	
	public DoubleMatrix1D toColt() {
		return this.getDm();
	}

	//@Override
	public Double0D get(int i) {
		return new Double0D(this.getDm().get(i));
	}

	//@Override
	public Object clone() throws CloneNotSupportedException {
		return new Double1D(this.getDm().copy());
	}
	
	public Double1D plus(Double1D o) {
		return new Double1D(AlgebraStatic.plusStatic(this.getDm(), o.getDm()));		
	}

	public boolean nonnegative() {
		for (int i=0; i<this.size(); i++) {
			if (this.getDm().get(i) < 0) {
				return false;
			}
		}
		return true;
	}

	//@Override
	public Double1D getAll(int... dis) {
		double[] d = this.value();
		double[] all = new double[dis.length];
		int ai = 0;
		for (int di : dis) {
			all[ai] = d[di]; 
			ai++;
		}
		return new Double1D(all);
	}

	public Double2D outer(Double1D o) {
		DoubleMatrix2D r = AlgebraStatic.multOuterStatic(this.getDm(), o.getDm());
		return new Double2D(r);
	}
	
	public Integer1D toInteger1D() {
		int[] is = new int[this.size()];
		for (int i=0; i<this.size(); i++) {
			is[i] = (int) this.value()[i];
		}
		return new Integer1D(is);
	}

	//@Override
	public Double1D cloneFromVector(Double1D v) {
		double[] ds = new double[this.size()];
		for (int i=0; i<this.size(); i++) {
			if (i >= v.size()) {
				return new Double1D(ds);
			}
			ds[i] = v.get(i).value();
		}
		return new Double1D(ds);
	}
	
	//@Override
	public String toString() {
		String s = "";
		String prefix = "";
		for (int i=0; i < this.size(); i++) {
			s = s + prefix + String.valueOf(this.getDm().get(i));
			prefix = " ";
		}
		return s;
	}	
	

	//@Override
	public boolean equals(Object o) {
		if (!(o instanceof Double1D)) {
			return false;
		}
		return this.value() == ((Double1D) o).value();
	}
	
	//@Override
	public int hashCode() {
		int h = 0;
		for (int i=0; i<this.size(); i++) {
			h = h + (int) this.value()[i];
		}
		return h;
	}

	public Double1D mult(double k) {
		return new Double1D(AlgebraStatic.multStatic(k, this.getDm()));
	}

	public Double1D minus(Double1D v) {
		return new Double1D(AlgebraStatic.minusStatic(this.getDm(), v.getDm()));
	}

	public Double1D mult(Double2D other) {
		return new Double1D(AlgebraStatic.multStatic(this.getDm(), other
				.getDm()));
	}

	public double mult(Double1D other) {
		return AlgebraStatic.multStatic(this.getDm(), other.getDm());
	}
	
	public Double1D minus(double k) {
		return new Double1D(AlgebraStatic.minusStatic(this.getDm(), k));
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
		return new Double1D(AlgebraStatic.plusStatic(this.getDm(), k));
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

	public Double2D toDouble2D(int h) {
		return new Double2D(AlgebraStatic.unvec(this.getDm(), h));
	}

	//@Override
	public Double1D rowVec() {
		return this;
	}
	
	public Iterator<Double0D> iterator() {
		return new Iterator<Double0D>() {
			private int curr = 0;

			public boolean hasNext() {
				return (this.curr < Double1D.this.length1D());
			}

			public Double0D next() {
				Double0D d = Double1D.this.rowVec().get(this.curr);
				this.curr++;
				return d;
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}
	
	//@Override
	public int length1D() {
		return this.size();
	}

	public Double1D divide(Double0D k) {
		return this.divide(k.value());
	}

	public Double1D divide(double k) {
		return new Double1D(AlgebraStatic.divideStatic(k, this.getDm()));
	}

	public Double0D mode() {
		Map<Double0D,Integer0D> counts = new HashMap<Double0D,Integer0D>();
		for (Double0D d : this) {
			if (!counts.containsKey(d)) {
				counts.put(d, new Integer0D(0));
			}
			counts.put(d, counts.get(d).plus(1));
		}
		Double0D mode = this.get(0);
		int max = 0;
		for (Double0D d : counts.keySet()) {
			if (counts.get(d).value() > max) {
				mode = d;
				max = counts.get(d).value();
			}
		}
		return mode;
	}

	public Double2D outer(Integer1D o) {
		double[] od = new double[o.size()];
		int[] ov = o.value();
		for (int i=0; i<ov.length; i++) {
			od[i] = ov[i];
		}
		DoubleMatrix2D r = AlgebraStatic.multOuterStatic(this.getDm(),
														 new DenseDoubleMatrix1D(od));
		return new Double2D(r);
	}
	
	public boolean hasBad() {
		for (int i=0; i<this.size(); i++) {
			if (Double.isNaN(this.value()[i]) || Double.isInfinite(this.value()[i])) {
				return true;
			}
		}
		return false;
	}

	public Double1D sqrt() {
		double[] v = this.value();
		for (int i=0; i<this.size(); i++) {
			v[i] = Math.sqrt(v[i]);
		}
		return new Double1D(v);
	}
}
