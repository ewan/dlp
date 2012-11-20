package org.jqgibbs.mathstat;

import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.CholeskyDecomposition;
import cern.jet.math.Functions;

public class Integer2D extends AbstractSequence<Integer2D, Integer1D> {

	int[][] ints;

	public Integer2D(Integer1D... i1Ds) {
		this.ints = new int[i1Ds.length][0];
		this.addAll(Arrays.asList(i1Ds));
	}

	public Integer2D(int[][] is) {
		this.ints = is;
	}

	public Integer2D(Integer[][] ints) {
		int[][] is = new int[ints.length][];
		for (int i = 0; i < ints.length; i++) {
			if (ints[i] == null) {
				throw new NullPointerException();
			}
			is[i] = new Integer1D((Integer[]) ints[i]).value();
		}
		this.ints = is;
	}

	public Integer2D(int m, int n, List<int[]> indices,
			AbstractSequence<?, Integer0D> values) {
		if (!(indices.size() == values.size())) {
			throw new IllegalArgumentException(
					"Length of indices was different from "
							+ "length of values");
		}
		this.ints = new int[m][n];
		for (int i = 0; i < indices.size(); i++) {
			try {
				this.ints[indices.get(i)[0]][indices.get(i)[1]] = values.get(i)
						.value();
			} catch (ArrayIndexOutOfBoundsException e) {
				throw new IllegalArgumentException("Index " + String.valueOf(i)
						+ " " + "did not have two elements");
			}
		}
	}

	public Integer2D(int m, int n, int... is) {
		if (!(is.length == m * n || is.length == 0)) {
			throw new IllegalArgumentException(String.format(
					"Incorrect number of values provided"
							+ " (expected %d, got %d)", m * n, is.length));
		}
		int[][] im = new int[m][n];
		if (is.length > 0) {
			int ii = 0;
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					im[i][j] = is[ii];
					ii++;
				}
			}
		}
		this.ints = im;
	}

	public boolean add(Integer1D i1D) {
		int[][] is = Arrays.copyOf(this.value(), this.size() + 1);
		is[this.size()] = i1D.value();
		this.ints = is;
		return true;
	}

	public boolean addAll(Collection<? extends Integer1D> c) {
		int[][] is = Arrays.copyOf(this.value(), this.size() + c.size());
		int i = this.size();
		for (Integer1D j : c) {
			if (j == null) {
				throw new NullPointerException();
			}
			is[i] = j.value();
			i = i + 1;
		}
		this.ints = is;
		return true;
	}

	public boolean contains(Object o) {
		int[] j;
		if (o == null) {
			throw new NullPointerException();
		} else if (o instanceof Integer1D) {
			j = ((Integer1D) o).value();
		} else if (o instanceof Integer[]) {
			j = new Integer1D((Integer[]) o).value();
		} else {
			return false;
		}
		for (int i = 0; i < this.size(); i++) {
			if (this.value()[i] == j) {
				return true;
			}
		}
		return false;
	}

	public Integer1D get(int i) {
		int[] ti = this.value()[i];
		int[] v = new int[this.numCols()];
		for (int j = 0; j < this.numCols(); j++) {
			v[j] = ti[j];
		}
		return new Integer1D(v);
	}

	@Override
	public Integer2D getAll(int... iis) {
		int[][] is = this.value();
		int[][] all = new int[iis.length][this.numCols()];
		int ai = 0;
		for (int ii : iis) {
			for (int j = 0; j < this.numCols(); j++) {
				all[ai][j] = is[ii][j];
			}
			ai++;
		}
		return new Integer2D(all);
	}

	@Override
	public Integer1D set(int i, Integer1D t) {
		if (i > this.size() || i < 0) {
			throw new IndexOutOfBoundsException(
					"Tried to set index out of range"); // FIXME
		}
		if (i == this.size()) {
			this.add(t); // FIXME - no dims checking??
		} else {
			this.value()[i] = t.value();
		}
		return t;
	}

	public int size() {
		return this.numRows();
	}

	@Override
	public Object clone() throws CloneNotSupportedException {
		return new Integer2D(Arrays.copyOf(this.value(), this.numRows()));
	}

	@Override
	public Integer2D cloneFromVector(Double1D v) {
		int[][] is = new int[this.numRows()][this.numCols()];
		int k = 0;
		for (int j = 0; j < this.numCols(); j++) {
			for (int i = 0; i < this.numRows(); i++) {
				if (k >= v.size()) {
					// FIXME - should warn or throw here (and of course in the
					// other case)...
					return new Integer2D(is);
				}
				is[i][j] = (int) v.get(k).value();
				k++;
			}
		}
		return new Integer2D(is);
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof Integer2D)) {
			return false;
		}
		return this.value() == ((Integer2D) o).value();
	}

	@Override
	public int hashCode() {
		int h = 0;
		for (int i = 0; i < this.size(); i++) {
			h = h + this.get(i).hashCode(); // FIXME - probably slow
		}
		return h;
	}

	public boolean isEmpty() {
		return this.size() == 0;
	}

	public Iterator<Integer1D> iterator() {
		return new Iterator<Integer1D>() {
			private int curr = 0;

			public boolean hasNext() {
				return (this.curr < Integer2D.this.size());
			}

			public Integer1D next() {
				Integer1D i = Integer2D.this.get(this.curr);
				this.curr++;
				return i;
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	public int[][] value() {
		return this.ints;
	}

	public Double2D toDouble2D() {
		double[][] ds = new double[this.numRows()][this.numCols()];
		for (int i = 0; i < this.numRows(); i++) {
			for (int j = 0; j < this.numCols(); j++) {
				ds[i][j] = (double) this.ints[i][j];
			}
		}
		return new Double2D(ds);
	}

	public Double1D mean() {
		double[] mean = new double[this.numCols()];
		double sum;
		for (int c = 0; c < this.numCols(); c++) {
			sum = 0;
			for (int i = 0; i < this.numRows(); i++) {
				sum = sum + this.value()[i][c];
			}
			mean[c] = sum / this.numRows();
		}
		return new Double1D(mean);
	}

	public Integer1D sum() {
		int[] sum = new int[this.numCols()];
		for (int c = 0; c < this.numCols(); c++) {
			sum[c] = 0;
			for (int i = 0; i < this.numRows(); i++) {
				sum[c] = sum[c] + this.value()[i][c];
			}
		}
		return new Integer1D(sum);
	}

	public boolean square() {
		return this.numCols() == this.numRows();
	}

	public int numCols() {
		if (this.numRows() == 0) {
			return 0;
		}
		return this.ints[0].length; // FIXME - no consistency checks in certain
									// places...
	}

	public int numRows() {
		return this.ints.length;
	}

	public Integer2D transpose() {
		int[][] is = new int[this.numCols()][this.numRows()];
		for (int i = 0; i < this.numRows(); i++) {
			for (int j = 0; j < this.numCols(); j++) {
				is[j][i] = this.ints[i][j];
			}
		}
		return new Integer2D(is);
	}

	public Double2D mult(Integer2D o) {
		return this.mult(o.toDouble2D());
	}

	public Double2D mult(Double2D o) {
		return new Double2D(AlgebraStatic.multStatic(this.toDouble2D().getDm(),
				o.getDm()));
	}

	public Double2D inverse() {
		return new Double2D(AlgebraStatic.inverseStatic(this.toDouble2D()
				.getDm()));
	}

	// public Double2D cholesky() {
	// if (!this.square()) {
	// throw new UnsupportedOperationException(
	// "Cannot perform Cholesky decomposition for non-square matrix");
	// }
	// CholeskyDecomposition cd = new CholeskyDecomposition(this.getDm());
	// // FIXME - find another test - this one gives (nearly-)false positives
	// // that
	// // seem to be due to sensitivity to roundoff error.
	// // if (!cd.isSymmetricPositiveDefinite()) {
	// // throw new UnsupportedOperationException(
	// // "Cannot perform Cholesky decomposition for non-positive "
	// // + "definite matrix");
	// // }
	// return new Double2D(cd.getL());
	// }
	//
	// public Double2D mult(Double2D other) {
	// return new Double2D(AlgebraStatic.multStatic(this.getDm(), other
	// .getDm()));
	// }
	//
	// public Double1D mult(Double1D other) {
	// return new Double1D(AlgebraStatic.multStatic(this.getDm(), other
	// .getDm()));
	// }
	//
	// public Double2D transpose() {
	// return new Double2D(AlgebraStatic.transposeStatic(this.getDm()).copy());
	// }
	//
	// public Double2D sumOuter() {
	// // Could iterate over rows of our matrix and use AlgebraStatic for speed
	// // As it is we're mixing levels...
	// // FIXME
	// DenseDoubleMatrix2D dm = new DenseDoubleMatrix2D(this.numCols(), this
	// .numCols());
	// Double2D outer;
	// for (Double1D xi : this) {
	// outer = xi.outer(xi);
	// dm.assign(outer.getDm(), Functions.plus);
	// }
	// return new Double2D(dm);
	// }
	//
	// public DoubleMatrix2D toColt() {
	// return this.getDm().copy();
	// }
	//
	// public Double2D mult(double k) {
	// return new Double2D(AlgebraStatic.multStatic(k, this.getDm()));
	// }
	//
	// public Double2D plus(Double2D o) {
	// return new Double2D(AlgebraStatic.plusStatic(this.getDm(), o.getDm()));
	// }
	//
	// public Double2D minus(Double1D v) {
	// return new Double2D(AlgebraStatic
	// .minusStatic(this.getDm(), v.getDm()));
	// }
	//
	// public Double2D minus(Double2D o) {
	// return new Double2D(AlgebraStatic.minusStatic(this.getDm(), o.getDm()));
	// }

	// public Double1D sum() {
	// double[] sum = new double[this.numCols()];
	// double[] d;
	// for (int c = 0; c < this.numCols(); c++) {
	// sum[c] = 0;
	// d = this.getDm().viewColumn(c).toArray();
	// for (int i = 0; i < d.length; i++) {
	// sum[c] = sum[c] + d[i];
	// }
	// }
	// return new Double1D(sum);
	// }

	@Override
	public String toString() {
		String s = "";
		String prefix = "";
		for (int i = 0; i < this.numRows(); i++) {
			for (int j = 0; j < this.numCols(); j++) {
				s = s + prefix + String.valueOf(this.ints[i][j]);
				prefix = " ";
			}
		}
		return s;
	}

	@Override
	public Double1D rowVec() {
		double[] ds = new double[this.numRows() * this.numCols()];
		int k = 0;
		int[][] iv = this.value();
		for (int i = 0; i < this.numRows(); i++) {
			for (int j = 0; j < this.numCols(); j++) {
				ds[k] = (double) iv[i][j];
				k++;
			}
		}
		return new Double1D(ds);
	}

	@Override
	public int length1D() {
		return this.numCols() * this.numRows();
	}

	public int getValue(int i, int j) {
		return this.value()[i][j];
	}

	public int[] getValue(int i) {
		int[] r = new int[this.numCols()];
		for (int j = 0; j < this.numCols(); j++) {
			r[j] = this.value()[i][j];
		}
		return (r);
	}

	public int set(int i, int j, int t) {
		if (i > this.size() || i < 0 || j > this.numCols() || j < 0) {
			throw new IndexOutOfBoundsException(
					"Tried to set index out of range"); // FIXME
		}
		if (j == this.numCols()) {
			int rows = this.numRows();
			if (i == this.numRows()) {
				rows += 1;
			}
			int[][] is = new int[rows][this.numCols() + 1];
			for (int m = 0; m < this.numRows(); m++) {
				for (int n = 0; n < this.numCols(); n++) {
					is[m][n] = this.getValue(m, n);
				}
			}
			is[i][j] = t;
			this.ints = is;
		} else {
			if (i == this.size()) {
				int[] row = new int[this.numCols()];
				row[j] = t;
				this.add(new Integer1D(row)); // FIXME - no dims checking??
			} else {
				this.value()[i][j] = t;
			}
		}
		return t;
	}

	public Integer1D mode() {
		int width = this.numCols();
		int[] v = new int[width];
		for (int i=0; i<this.numCols(); i++) {
			int[] col = new int[this.numRows()];
			for (int j=0; j<this.numRows(); j++) {
				col[j] = this.ints[j][i];
			}
			v[i] = new Integer1D(col).mode().value();
		}
		return new Integer1D(v);
	}

	// public double det() {
	// return AlgebraStatic.detStatic(this.getDm());
	// }
	//
	// public Double2D cov() {
	// Double2D cthis = this.minus(this.mean());
	// return cthis.transpose().mult(cthis).mult(1/((double) this.size()-1));
	// }
	//
	// public Double0D max() {
	// double max = -Double.MAX_VALUE;
	// for (int i=0; i<this.numRows(); i++) {
	// for (int j=0; j<this.numCols(); j++) {
	// if (this.getDm().get(i,j) > max) {
	// max = this.getDm().get(i,j);
	// }
	// }
	// }
	// return new Double0D(max);
	// }
	//
	// public Double1D getCol(int i) {
	// return new Double1D(this.getDm().viewColumn(i));
	// }
	//
	// public Double2D mult(Double0D k) {
	// return this.mult(k.value());
	// }
	//
	// public Double2D divide(Double0D k) {
	// return this.divide(k.value());
	// }
	//
	// public Double2D divide(double k) {
	// return new Double2D(AlgebraStatic.divideStatic(k, this.getDm()));
	// }
	//
	// public Double2D diag() {
	// // FIXME - fail if not square
	// double[][] ds = new double[this.numRows()][this.numCols()];
	// for (int i=0; i<this.numRows(); i++) {
	// ds[i][i] = this.value()[i][i];
	// }
	// Double2D d = new Double2D(ds);
	// return d;
	// }
	//	
	// public Double1D vec() {
	// double[] ds = new double[this.numRows()*this.numCols()];
	// int k=0;
	// for (int i=0; i<this.numCols(); i++) {
	// for (int j=0; j<this.numRows(); j++) {
	// ds[k] = this.value()[j][i];
	// k++;
	// }
	// }
	// return new Double1D(ds);
	// }
}
