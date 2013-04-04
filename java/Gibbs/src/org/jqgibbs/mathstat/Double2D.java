package org.jqgibbs.mathstat;

import java.util.Arrays;
import java.util.List;

import org.jqgibbs.Flattenable;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.CholeskyDecomposition;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import cern.jet.math.Functions;

public class Double2D implements Flattenable,Cloneable {
	public static final double MAX_MATRIX_COND_NUM = 1e16;
	
	private DoubleMatrix2D dm;

	public Double2D(double[][] ds) {
		this.dm = new DenseDoubleMatrix2D(ds);
	}

	public Double2D(DoubleMatrix2D dm) {
		if (dm instanceof DenseDoubleMatrix2D) { // FIXME - this colt class (minimally) has its own
			this.dm = dm.viewSelection(null,null);    // multiplication procedure that's incompatible w
												 // the other matrix being a selection view. seems dumb..
												 // what we're doing here is basically casting to a selection
												 // view so that we never use that
		} else {
			this.dm = dm;
		} 
	}

	public Double2D(int m, int n, List<int[]> indices,
			Double1D values) {
		if (!(indices.size() == values.size())) {
			throw new IllegalArgumentException(
					"Length of indices was different from "
							+ "length of values");
		}
		this.dm = new SparseDoubleMatrix2D(m, n);
		// Could construct the double[][] matrix then make the SDM out of it -
		// FIXME
		for (int i = 0; i < indices.size(); i++) {
			try {
				this.dm.set(indices.get(i)[0], indices.get(i)[1],
						values.get(i).value());
			} catch (ArrayIndexOutOfBoundsException e) {
				throw new IllegalArgumentException("Index " + String.valueOf(i)
						+ " " + "did not have two elements");
			}
		}
	}

	public Double2D(int m, int n) {
		double[][] dm = new double[m][n];
		this.dm = (new DenseDoubleMatrix2D(dm)).viewSelection(null,null); // FIXME - see above
	}
	
	public int length1D() {
		return this.numCols() * this.numRows();
	}

	@Override
	public String toString() {
		String s = "";
		String prefix = "";
		for (int i = 0; i < this.numRows(); i++) {
			for (int j = 0; j < this.numCols(); j++) {
				s = s + prefix + String.valueOf(this.dm.get(i, j));
				prefix = " ";
			}
		}
		return s;
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof Double2D)) {
			return false;
		}
		return this.dm.equals(((Double2D) o).dm);		
	}

	@Override
	public int hashCode() {
		return this.dm.hashCode();
	}	
	
	@Override
	public Object clone() {
		return this.copy();
	}
	
	public DoubleMatrix2D getDm() {
		return this.dm;
	}
	
	// FIXME
	public DoubleMatrix2D toColt() {
		return this.dm.copy();
	}
	
	public boolean add(Double1D d1D) {
		double[][] ds = Arrays.copyOf(this.dm.toArray(), this.size() + 1);
		ds[this.size()] = d1D.value();
		this.dm = new DenseDoubleMatrix2D(ds);
		return true;
	}

	public int size() {
		return this.numRows();
	}

	public double[][] value() {
		return this.dm.toArray();
	}

//	public Double2D getAll(int... dis) {
//		double[][] d = this.value();
//		double[][] all = new double[dis.length][this.numCols()];
//		int ai = 0;
//		for (int di : dis) {
//			for (int j = 0; j < this.numCols(); j++) {
//				all[ai][j] = d[di][j];
//			}
//			ai++;
//		}
//		DoubleMatrix2D dm = new DenseDoubleMatrix2D(dis.length, this.numCols());
//		dm.assign(all);
//		return new Double2D(dm);
//	}

	// FIXME - check pour la seguridade
	public Double2D getAll(int... dis) {
		DoubleMatrix2D dm = this.dm.viewSelection(dis, null);
		return new Double2D(dm);
	}
	
	public Double1D colVec() {
		double[] ds = new double[this.numRows() * this.numCols()];
		int k = 0;
		for (int i = 0; i < this.numCols(); i++) {
			for (int j = 0; j < this.numRows(); j++) {
				ds[k] = this.value()[j][i];
				k++;
			}
		}
		return new Double1D(ds);
	}

	public Double1D rowVec() {
		double[] ds = new double[this.numRows() * this.numCols()];
		int k = 0;
		for (int i = 0; i < this.numRows(); i++) {
			for (int j = 0; j < this.numCols(); j++) {
				ds[k] = this.value()[i][j];
				k++;
			}
		}
		return new Double1D(ds);
	}
	
	public int numCols() {
		return this.dm.columns();
	}

	public int numRows() {
		return this.dm.rows();
	}

	public Double1D get(int i) {
		return new Double1D(this.dm.viewRow(i));
	}

	public Double1D getCol(int i) {
		return new Double1D(this.dm.viewColumn(i));
	}
	
	public Double2D getColAll(int... is) {
		double[][] d = new double[this.numRows()][is.length];
		for (int n = 0; n < is.length; n++) {
			for (int i = 0; i < this.numRows(); i++) {
				d[i][n] = this.dm.get(i, is[n]);
			}
		}
		return new Double2D(d);
	}
	
	// Convenience methods
	
	public Double1D mean() {
		double[] mean = new double[this.numCols()];
		double[] d;
		double sum;
		for (int c = 0; c < this.numCols(); c++) {
			sum = 0;
			d = this.dm.viewColumn(c).toArray();
			for (int i = 0; i < d.length; i++) {
				sum = sum + d[i];
			}
			mean[c] = sum / d.length;
		}
		return new Double1D(mean);
	}

	public boolean square() {
		return this.dm.columns() == this.dm.rows();
	}
	
	public Double2D inverse() {
		return new Double2D(AlgebraStatic.inverseStatic(this.dm));
	}

	public Double2D cholesky() {
		if (!this.square()) {
			throw new UnsupportedOperationException(
					"Cannot perform Cholesky decomposition for non-square matrix");
		}
		CholeskyDecomposition cd = new CholeskyDecomposition(this.dm);
		// FIXME - find another test - this one gives (nearly-)false positives
		// that
		// seem to be due to sensitivity to roundoff error.
		// if (!cd.isSymmetricPositiveDefinite()) {
		// throw new UnsupportedOperationException(
		// "Cannot perform Cholesky decomposition for non-positive "
		// + "definite matrix");
		// }
		return new Double2D(cd.getL());
	}

	public Double2D mult(Double2D other) {
		return new Double2D(this.dm.zMult(other.dm, null));
//		return new Double2D(AlgebraStatic.multStatic(this.dm, other
//				.dm));
	}

	public Double1D mult(Double1D other) {
		return new Double1D(AlgebraStatic.multStatic(this.dm, other
				.getDm()));
	}

	public Double2D transpose() {
		return new Double2D(AlgebraStatic.transposeStatic(this.dm).copy());
	}

	public Double2D mult(double k) {
		return new Double2D(AlgebraStatic.multStatic(k, this.dm));
	}

	public Double2D plus(Double2D o) {
		DoubleMatrix2D dm_copy = this.dm.copy().viewSelection(null,null);
		return new Double2D(dm_copy.assign(o.dm, Functions.plus));
//		return new Double2D(AlgebraStatic.plusStatic(this.dm, o.dm));
	}

	public Double2D minus(Double1D v) {
		return new Double2D(AlgebraStatic.minusStatic(this.dm, v.getDm()));
	}

	public Double2D minus(Double2D o) {
		DoubleMatrix2D dm_copy = this.dm.copy().viewSelection(null,null);
		return new Double2D(dm_copy.assign(o.dm, Functions.minus));
//		return new Double2D(AlgebraStatic.minusStatic(this.dm, o.dm));
	}

	public Double1D sum() {
		double[] sum = new double[this.numCols()];
		double[] d;
		for (int c = 0; c < this.numCols(); c++) {
			sum[c] = 0;
			d = this.dm.viewColumn(c).toArray();
			for (int i = 0; i < d.length; i++) {
				sum[c] = sum[c] + d[i];
			}
		}
		return new Double1D(sum);
	}
	
	public double det() {
		return AlgebraStatic.detStatic(this.dm);
	}

	public Double2D cov() {
		Double2D cthis = this.minus(this.mean());
		return cthis.transpose().mult(cthis).mult(
				1 / ((double) this.size() - 1));
	}

	public Double0D max() {
		double max = -Double.MAX_VALUE;
		for (int i = 0; i < this.numRows(); i++) {
			for (int j = 0; j < this.numCols(); j++) {
				if (this.dm.get(i, j) > max) {
					max = this.dm.get(i, j);
				}
			}
		}
		return new Double0D(max);
	}

	public Double2D mult(Double0D k) {
		return this.mult(k.value());
	}

	public Double2D divide(Double0D k) {
		return this.divide(k.value());
	}

	public Double2D divide(double k) {
		return new Double2D(AlgebraStatic.divideStatic(k, this.dm));
	}

	public Double2D diag() {
		// FIXME - fail if not square
		double[][] ds = new double[this.numRows()][this.numCols()];
		for (int i = 0; i < this.numRows(); i++) {
			ds[i][i] = this.value()[i][i];
		}
		Double2D d = new Double2D(ds);
		return d;
	}


	public Double2D kron(Double2D r) {
		return new Double2D(AlgebraStatic.kronStatic(this.dm, r.dm));
	}

	public static Double2D ident(int h) {
		Double2D n = new Double2D(h, h);
		for (int i = 0; i < h; i++) {
			n.dm.setQuick(i, i, 1);
		}
		return n;
	}
	
	public Double0D trace() {
		// FIXME - fail if not square
		double tr = 0;
		for (int i=0; i<this.numRows(); i++) {
			tr += this.value()[i][i];
		}
		return new Double0D(tr);
	}

	public Double1D mult(Integer1D other) {
		double[] od = new double[other.size()];
		for (int i = 0; i < other.size(); i++) {
			od[i] = other.value()[i];
		}
		return this.mult(new Double1D(od));
	}
	
	public boolean isWellConditioned() {
		SingularValueDecomposition svd = new SingularValueDecomposition(this.toColt());
		return svd.cond() <= Double2D.MAX_MATRIX_COND_NUM;
	}

	public Double2D divide(Double1D v) {
		double[][] td = this.value();
		double[][] sd = new double[this.numRows()][this.numCols()];
		for (int i=0; i<this.numRows(); i++) {
			for (int j=0; j<this.numCols(); j++) {
				sd[i][j] = td[i][j]/v.get(j).value();
			}
		}
		return new Double2D(sd);
	}

	public Double2D transposeMult(Double2D other) {
		return new Double2D(this.dm.zMult(other.dm, null, 1, 1, true, false));
	}

	public void plusEquals(Double2D o) {
		this.dm.assign(o.dm, Functions.plus);
	}

	public Double2D copy() {
		return new Double2D(this.dm.copy());
	}

	public Double1D transposeMult(Double1D other) {
		return new Double1D(this.dm.zMult(other.getDm(), null, 1, 1, true));
	}

	public Double2D multTranspose(Double2D other) {
		return new Double2D(this.dm.zMult(other.dm, null, 1, 1, false, true));
	}

	public Double2D subMatrix(int[] row_indices, int[] col_indices) {
		return new Double2D(this.dm.viewSelection(row_indices, col_indices));
	}

	public Double2D sparsifyRows(int[] row_indices) {
		int p = this.dm.rows();
		DoubleMatrix2D reduction_matrix = new SparseDoubleMatrix2D(p,p).viewSelection(null, null);
		for (int i=0;i<row_indices.length;i++) {
			reduction_matrix.setQuick(row_indices[i], row_indices[i], 1);
		}
		return new Double2D(reduction_matrix.zMult(this.dm, null));
	}

	public Double2D blockBDInv(int[] indices_lower, int[] indices_upper) {
		DoubleMatrix2D B = this.dm.viewSelection(indices_upper, indices_lower); 
		DoubleMatrix2D D = this.dm.viewSelection(indices_lower, indices_lower); 
		return new Double2D(B.zMult(AlgebraStatic.inverseStatic(D),null));
	}
	
	public Double2D schur(int[] indices_lower, int[] indices_upper) {
		DoubleMatrix2D A = this.dm.viewSelection(indices_upper, indices_upper).copy(); 
		DoubleMatrix2D B = this.dm.viewSelection(indices_upper, indices_lower); 
		DoubleMatrix2D C = this.dm.viewSelection(indices_lower, indices_upper); 
		DoubleMatrix2D D = this.dm.viewSelection(indices_lower, indices_lower); 
		DoubleMatrix2D BDI = B.zMult(AlgebraStatic.inverseStatic(D),null).viewSelection(null,null);
		A.assign(BDI.zMult(C, null),Functions.minus);
		return new Double2D(A);
	}

	public Double2D interleaveRows(Double2D o, int[] this_indices, int[] other_indices) {
		int bigRows = this.numRows()+o.numRows();
		int bigCols = this.numCols();
		DenseDoubleMatrix2D big = new DenseDoubleMatrix2D(bigRows,bigCols);
		int i_this = 0;
		for (int i : this_indices) {
			for (int j=0; j<bigCols; j++) {
				big.setQuick(i, j, this.dm.getQuick(i_this, j));
			}
			i_this++;
		}
		int i_o = 0;
		for (int i : other_indices) {
			for (int j=0; j<bigCols; j++) {
				big.setQuick(i, j, o.dm.getQuick(i_o, j));
			}
			i_o++;
		}
		return new Double2D(big);
	}	

}
