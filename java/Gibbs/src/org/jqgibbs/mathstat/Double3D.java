package org.jqgibbs.mathstat;

import java.util.Arrays;
import java.util.Collection;

import org.jqgibbs.Flattenable;

import cern.colt.matrix.DoubleMatrix3D;
import cern.colt.matrix.impl.DenseDoubleMatrix3D;

public class Double3D implements Flattenable,Cloneable {
	private DoubleMatrix3D dm;
	
	public Double3D(Double2D... d2Ds) {
		this.dm = new DenseDoubleMatrix3D(new double[0][0][0]);
		this.addAll(Arrays.asList(d2Ds));
	}
	
	public Double3D(double[][][] ds) {
		this.dm = new DenseDoubleMatrix3D(ds);
	}

	public Double3D(DoubleMatrix3D dm) {
		this.dm = dm;
	}
	
	public Double3D(int l, int m, int n) {
		double[][][] dm = new double[l][m][n];
		this.dm = new DenseDoubleMatrix3D(dm);
	}

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
	
	public int length1D() {
		return this.size()*this.numCols()*this.numRows();
	}	
	
	@Override
	public String toString() {
		String s = "";
		String prefix = "";
		for (int i=0; i < this.size(); i++) {
			for (int j=0; j < this.numRows(); j++) {
				for (int k=0; k < this.numCols(); k++) {
					s = s + prefix + String.valueOf(this.dm.get(i, j, k));
					prefix = " ";
				}
			}
		}
		return s;
	}	

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof Double3D)) {
			return false;
		}
		return this.dm.equals(((Double3D) o).dm);
	}

	@Override
	public int hashCode() {
		return this.dm.hashCode();
	}	
	
	@Override
	public Object clone() {
		return new Double3D(this.value());
	}
	
	public boolean add(Double2D d2D) {
		double[][][] ds = Arrays.copyOf(this.dm.toArray(), this.size()+1);
		ds[this.size()] = d2D.value();
		this.dm = new DenseDoubleMatrix3D(ds);
		return true;
	}

	// FIXME
	public boolean addAll(Collection<? extends Double2D> c) {
		double[][][] ds = Arrays.copyOf(this.dm.toArray(),
				this.size()+c.size());
		int i = this.size();
		for (Double2D d : c) {
			if (d == null) {
				throw new NullPointerException();
			}			
			ds[i] = d.value();
			i = i + 1;
		}
		this.dm = new DenseDoubleMatrix3D(ds);
		return true;
	}

	public int size() {
		return this.dm.slices();
	}
	
	public int numRows() {
		return this.dm.rows();
	}
	
	public int numCols() {
		return this.dm.columns();
	}

	public double[][][] value() {
		return this.dm.toArray();
	}

	public Double2D get(int i) {
		return new Double2D(this.dm.viewSlice(i));
	}
	
	// FIXME - to copy or not???
	public DoubleMatrix3D toColt() {
		return this.dm.copy();
	}

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

	// FIXME - is this the smartest way to do this?
	public Double2D set(int i, Double2D t) {
		if (i > this.size() || i < 0) {
			throw new IndexOutOfBoundsException("Tried to set index out of range"); // FIXME
		}
		if (i == this.size()) {
			this.add(t); // FIXME - no dims checking??
		} else {
			for (int j=0; j<this.numRows(); j++) {
				for (int k=0; k<this.numCols(); k++) {
					this.dm.set(i, j, k, t.toColt().get(j, k));
				}
			}
		}
		return t;
	}

}
