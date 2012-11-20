package org.jqgibbs.mathstat;

import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import cern.colt.matrix.DoubleMatrix3D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix3D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;

public class Double3D extends AbstractSequence<Double3D, Double2D> {
	private DoubleMatrix3D dm;
	
	public Double3D(Double2D... d2Ds) {
		this.setDm(new double[0][0][0]);
		this.addAll(Arrays.asList(d2Ds));
	}
	
	public Double3D(double[][][] ds) {
		this.setDm(ds);
	}

	public Double3D(DoubleMatrix3D dm) {
		this.setDm(dm);
	}
	
	public Double3D(int l, int m, int n, double... ds) {
		if (!(ds.length == l*m*n || ds.length == 0)) {
			throw new IllegalArgumentException(String.format(
					"Incorrect number of values provided"
							+ " (expected %d, got %d)", l*m*n, ds.length));
		}
		double[][][] dm = new double[l][m][n];
		if (ds.length > 0) {
		int di = 0;
		for (int i = 0; i < l; i++) {
			for (int j = 0; j < m; j++) {
					for (int k = 0; k < n; k++) {
						dm[i][j][k] = ds[di];
						di++;
					}
				}
			}
		}
		this.setDm(new DenseDoubleMatrix3D(dm));
	}

	private synchronized void setDm(DoubleMatrix3D dm) {
		this.dm = dm;
	}
	
	private void setDm(double[][][] ds) {
		this.setDm(new DenseDoubleMatrix3D(ds));
	}

	private synchronized DoubleMatrix3D getDm() {
		return this.dm;
	}

	public boolean add(Double2D d2D) {
		double[][][] ds = Arrays.copyOf(this.getDm().toArray(), this.size()+1);
		ds[this.size()] = d2D.value();
		this.setDm(ds);
		return true;
	}

	public boolean addAll(Collection<? extends Double2D> c) {
		double[][][] ds = Arrays.copyOf(this.getDm().toArray(),
				this.size()+c.size());
		int i = this.size();
		for (Double2D d : c) {
			if (d == null) {
				throw new NullPointerException();
			}			
			ds[i] = d.value();
			i = i + 1;
		}
		this.setDm(ds);
		return true;
	}

	public boolean contains(Object o) {
		double[][] d;
		if (o == null) {
			throw new NullPointerException();
		} else if (o instanceof Double1D) {
			d = ((Double2D) o).value();
		} else if (o instanceof Double[][]) {
			d = new Double2D((Double[][]) o).value();
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

	public boolean containsAll(Collection<?> c) {
		for (Object o : c) {
			if (!(this.contains(o))) {
				return false;
			}
		}
		return true;
	}

	public boolean isEmpty() {
		return this.size() == 0;
	}

	
	public Iterator<Double2D> iterator() {
		return new Iterator<Double2D>() {
			private int curr = 0;
			public boolean hasNext() {
				return (this.curr < Double3D.this.size());
			}
			public Double2D next() {
				Double2D d = Double3D.this.get(this.curr);
				this.curr++;
				return d;
			}
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	public int size() {
		return this.getDm().slices();
	}
	
	public int numRows() {
		return this.getDm().rows();
	}
	
	public int numCols() {
		return this.getDm().columns();
	}

	public Object[] toArray() {
		Double2D[] ds = new Double2D[this.size()];
		int i=0;
		for (Double2D d : this) {
			ds[i] = d;
			i++;
		}
		return ds;
	}

	@SuppressWarnings("unchecked")
	public <T> T[] toArray(T[] a) {
		if (!(a.getClass().isAssignableFrom(Double2D.class))) {
			throw new ArrayStoreException();
		}
		if (!(a.length < this.size())) {
			a = (T[]) new Double2D[this.size()];
		}
		int i=0;
		for (Double2D d : this) {
			a[i] = (T) d;
			i++;
		}
		if (a.length > this.size()) {
			a[i] = null;
		}		
		return a;
	}
	
	public void clear() {
		throw new UnsupportedOperationException();
	}
	
	public boolean remove(Object o) {
		throw new UnsupportedOperationException();
	}

	public boolean removeAll(Collection<?> c) {
		throw new UnsupportedOperationException();
	}

	public boolean retainAll(Collection<?> c) {
		throw new UnsupportedOperationException();
	}

	public double[][][] value() {
		return this.getDm().toArray();
	}

	public Double2D get(int i) {
		return new Double2D(this.getDm().viewSlice(i));
	}
	
	public DoubleMatrix3D toColt() {
		return this.getDm();
	}

	@Override
	public Object clone() throws CloneNotSupportedException {
		return new Double3D(this.toColt().copy());
	}

	@Override
	public Double3D getAll(int... dis) {
		double[][][] d = this.value();
		double[][][] all = new double[dis.length][d[0].length][d[0][0].length];
		int ai = 0;
		for (int di : dis) {
			for (int j=0; j<d[0].length; j++) {
				for (int k=0; k<d[0][0].length; k++) {
					all[ai][j][k] = d[di][j][k];
				}
			}
			ai++;
		}
		return new Double3D(all);
	}

	@Override
	public Double3D cloneFromVector(Double1D v) {
		double[][][] ds = new double[this.size()][this.numRows()][this.numCols()];
		int l = 0;
		for (int i=0; i<this.size(); i++) {
			for (int j=0; j<this.numRows(); j++) {
				for (int k=0; k<this.numCols(); k++) {
					if (l >= v.size()) {
						// FIXME - should warn or throw here (and of course in the
						// other case)...
						return new Double3D(ds);
					}
					double a = v.get(l).value();
					ds[i][j][k] = a; 
					l++;
				}
			}
		}
		return new Double3D(ds);
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof Double3D)) {
			return false;
		}
		return this.value() == ((Double3D) o).value();
	}

	@Override
	public int hashCode() {
		int h = 0;
		for (int i=0; i<this.size(); i++) {
			h = h + this.get(i).hashCode();
		}
		return h;
	}

	@Override
	public Double2D set(int i, Double2D t) {
		if (i > this.size() || i < 0) {
			throw new IndexOutOfBoundsException("Tried to set index out of range"); // FIXME
		}
		if (i == this.size()) {
			this.add(t); // FIXME - no dims checking??
		} else {
			for (int j=0; j<this.numRows(); j++) {
				for (int k=0; k<this.numCols(); k++) {
					this.getDm().set(i, j, k, t.toColt().get(j, k));
				}
			}
		}
		return t;
	}
	

	@Override
	public String toString() {
		String s = "";
		String prefix = "";
		for (int i=0; i < this.size(); i++) {
			for (int j=0; j < this.numRows(); j++) {
				for (int k=0; k < this.numCols(); k++) {
					s = s + prefix + String.valueOf(this.getDm().get(i, j, k));
					prefix = " ";
				}
			}
		}
		return s;
	}

	@Override
	public Double1D rowVec() {
		double dd[] = new double[this.size()*this.numRows()*this.numCols()];
		int k = 0;
		for (int i=0; i<this.size(); i++) {
			double dv[] = this.get(i).rowVec().value();
			for (int j=0; j<dv.length; j++) {
				dd[k] = dv[j];
				k++;
			}
		}
		return new Double1D(dd);
	}	
	
	@Override
	public int length1D() {
		return this.size()*this.numCols()*this.numRows();
	}

	public Double3D cycle() {
		return this;
	}
}
