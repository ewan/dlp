package org.jqgibbs.mathstat;

import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.CholeskyDecomposition;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import cern.jet.math.Functions;

public class Double2D extends AbstractSequence<Double2D, Double1D> {
	public static final double MAX_MATRIX_COND_NUM = 1e16;
	
	private DoubleMatrix2D dm;

	public Double2D(Double1D... d1Ds) {
		this.setDm(new double[0][0]);
		this.addAll(Arrays.asList(d1Ds));
	}

	public Double2D(double[][] ds) {
		this.setDm(ds);
	}

	public Double2D(DoubleMatrix2D dm) {
		this.setDm(dm);
	}

	public Double2D(Double[][] doubles) {
		double[][] ds = new double[doubles.length][];
		for (int i = 0; i < doubles.length; i++) {
			if (doubles[i] == null) {
				throw new NullPointerException();
			}
			ds[i] = new Double1D((Double[]) doubles[i]).value();
		}
		this.setDm(ds);
	}

	public Double2D(int m, int n, List<int[]> indices,
			AbstractSequence<?, Double0D> values) {
		if (!(indices.size() == values.size())) {
			throw new IllegalArgumentException(
					"Length of indices was different from "
							+ "length of values");
		}
		this.setDm(new SparseDoubleMatrix2D(m, n));
		// Could construct the double[][] matrix then make the SDM out of it -
		// FIXME
		for (int i = 0; i < indices.size(); i++) {
			try {
				this.getDm().set(indices.get(i)[0], indices.get(i)[1],
						values.get(i).value());
			} catch (ArrayIndexOutOfBoundsException e) {
				throw new IllegalArgumentException("Index " + String.valueOf(i)
						+ " " + "did not have two elements");
			}
		}
	}

	public Double2D(int m, int n, double... ds) {
		if (!(ds.length == m * n)) {
			if (ds.length == 0) {
				ds = new double[m * n];
			} else {
				throw new IllegalArgumentException(String.format(
						"Incorrect number of values provided"
								+ " (expected %d, got %d)", m * n, ds.length));
			}
		}
		double[][] dm = new double[m][n];
		int di = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				dm[i][j] = ds[di];
				di++;
			}
		}
		this.setDm(new DenseDoubleMatrix2D(dm));
	}

	private synchronized void setDm(DoubleMatrix2D dm) {
		this.dm = dm;
	}

	private void setDm(double[][] ds) {
		this.setDm(new DenseDoubleMatrix2D(ds));
	}

	synchronized DoubleMatrix2D getDm() {
		return this.dm;
	}

	@Override
	public Double3D sequence() {
		return new Double3D(this);
	}

	public boolean add(Double1D d1D) {
		double[][] ds = Arrays.copyOf(this.getDm().toArray(), this.size() + 1);
		ds[this.size()] = d1D.value();
		this.setDm(ds);
		return true;
	}

	public boolean addAll(Collection<? extends Double1D> c) {
		double[][] ds = Arrays.copyOf(this.getDm().toArray(), this.size()
				+ c.size());
		int i = this.size();
		for (Double1D d : c) {
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
		double[] d;
		if (o == null) {
			throw new NullPointerException();
		} else if (o instanceof Double1D) {
			d = ((Double1D) o).value();
		} else if (o instanceof Double[]) {
			d = new Double1D((Double[]) o).value();
		} else {
			return false;
		}
		for (int i = 0; i < this.size(); i++) {
			if (this.get(i).value() == d) {
				return true;
			}
		}
		return false;
	}

	public boolean isEmpty() {
		return this.size() == 0;
	}

	public Iterator<Double1D> iterator() {
		return new Iterator<Double1D>() {
			private int curr = 0;

			public boolean hasNext() {
				return (this.curr < Double2D.this.size());
			}

			public Double1D next() {
				Double1D d = Double2D.this.get(this.curr);
				this.curr++;
				return d;
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	public int size() {
		return this.numRows();
	}

	public double[][] value() {
		return this.getDm().toArray();
	}

	public Double1D get(int i) {
		// FIXME?
		// The Colt embedded in the returned Double1D will be backed by the same
		// matrix as the current Colt, but this isn't unsafe as far as I can
		// tell, unless those algebra methods are doing something sneaky,
		// because every destructive method copies the underlying array.
		return new Double1D(this.getDm().viewRow(i));
	}

	public Double1D mean() {
		double[] mean = new double[this.numCols()];
		double[] d;
		double sum;
		for (int c = 0; c < this.numCols(); c++) {
			sum = 0;
			d = this.getDm().viewColumn(c).toArray();
			for (int i = 0; i < d.length; i++) {
				sum = sum + d[i];
			}
			mean[c] = sum / d.length;
		}
		return new Double1D(mean);
	}

	public Double1D mode() {
		double[] mode = new double[this.numCols()];
		for (int c = 0; c < this.numCols(); c++) {
			mode[c] = this.getCol(c).mode().value();
		}
		return new Double1D(mode);
	}

	public boolean square() {
		return this.getDm().columns() == this.getDm().rows();
	}

	public int numCols() {
		return this.getDm().columns();
	}

	public int numRows() {
		return this.getDm().rows();
	}

	public Double2D inverse() {
		return new Double2D(AlgebraStatic.inverseStatic(this.dm));
	}

	public Double2D cholesky() {
		if (!this.square()) {
			throw new UnsupportedOperationException(
					"Cannot perform Cholesky decomposition for non-square matrix");
		}
		CholeskyDecomposition cd = new CholeskyDecomposition(this.getDm());
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
		return new Double2D(AlgebraStatic.multStatic(this.getDm(), other
				.getDm()));
	}

	public Double1D mult(Double1D other) {
		return new Double1D(AlgebraStatic.multStatic(this.getDm(), other
				.getDm()));
	}

	public Double2D transpose() {
		return new Double2D(AlgebraStatic.transposeStatic(this.getDm()).copy());
	}

	public DoubleMatrix2D toColt() {
		return this.getDm().copy();
	}

	@Override
	public Object clone() throws CloneNotSupportedException {
		return new Double2D(this.getDm().copy());
	}

	public Double2D mult(double k) {
		return new Double2D(AlgebraStatic.multStatic(k, this.getDm()));
	}

	public Double2D plus(Double2D o) {
		return new Double2D(AlgebraStatic.plusStatic(this.getDm(), o.getDm()));
	}

	public Double2D minus(Double1D v) {
		return new Double2D(AlgebraStatic.minusStatic(this.getDm(), v.getDm()));
	}

	public Double2D minus(Double2D o) {
		return new Double2D(AlgebraStatic.minusStatic(this.getDm(), o.getDm()));
	}

	@Override
	public Double2D getAll(int... dis) {
		double[][] d = this.value();
		double[][] all = new double[dis.length][this.numCols()];
		int ai = 0;
		for (int di : dis) {
			for (int j = 0; j < this.numCols(); j++) {
				all[ai][j] = d[di][j];
			}
			ai++;
		}
		DoubleMatrix2D dm = new DenseDoubleMatrix2D(dis.length, this.numCols());
		dm.assign(all);
		return new Double2D(dm);
	}

	public Double1D sum() {
		double[] sum = new double[this.numCols()];
		double[] d;
		for (int c = 0; c < this.numCols(); c++) {
			sum[c] = 0;
			d = this.getDm().viewColumn(c).toArray();
			for (int i = 0; i < d.length; i++) {
				sum[c] = sum[c] + d[i];
			}
		}
		return new Double1D(sum);
	}

	@Override
	public Double2D cloneFromVector(Double1D v) {
		double[][] ds = new double[this.numRows()][this.numCols()];
		int k = 0;
		for (int i = 0; i < this.numRows(); i++) {
			for (int j = 0; j < this.numCols(); j++) {
				if (k >= v.size()) {
					// FIXME - should warn or throw here (and of course in the
					// other case)...
					return new Double2D(ds);
				}
				ds[i][j] = v.get(k).value();
				k++;
			}
		}
		return new Double2D(ds);
	}

	@Override
	public String toString() {
		String s = "";
		String prefix = "";
		for (int i = 0; i < this.numRows(); i++) {
			for (int j = 0; j < this.numCols(); j++) {
				s = s + prefix + String.valueOf(this.getDm().get(i, j));
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
		return this.value() == ((Double2D) o).value();
	}

	@Override
	public int hashCode() {
		int h = 0;
		for (int i = 0; i < this.size(); i++) {
			h = h + this.get(i).hashCode();
		}
		return h;
	}

	@Override
	public Double1D set(int i, Double1D t) {
		if (i > this.size() || i < 0) {
			throw new IndexOutOfBoundsException(
					"Tried to set index out of range"); // FIXME
		}
		if (i == this.size()) {
			this.add(t); // FIXME - no dims checking??
		} else {
			for (int j = 0; j < this.numCols(); j++) {
				this.getDm().set(i, j, t.getDm().get(j));
			}
		}
		return t;
	}

	public double det() {
		return AlgebraStatic.detStatic(this.getDm());
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
				if (this.getDm().get(i, j) > max) {
					max = this.getDm().get(i, j);
				}
			}
		}
		return new Double0D(max);
	}

	public Double1D getCol(int i) {
		return new Double1D(this.getDm().viewColumn(i));
	}

	public Double2D mult(Double0D k) {
		return this.mult(k.value());
	}

	public Double2D divide(Double0D k) {
		return this.divide(k.value());
	}

	public Double2D divide(double k) {
		return new Double2D(AlgebraStatic.divideStatic(k, this.getDm()));
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

	public Double2D kron(Double2D r) {
		return new Double2D(AlgebraStatic.kronStatic(this.getDm(), r.getDm()));
	}

	public static Double2D ident(int h) {
		Double2D n = new Double2D(h, h);
		for (int i = 0; i < h; i++) {
			n.getDm().setQuick(i, i, 1);
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

	public Double2D getColAll(int... is) {
		double[][] d = new double[this.numRows()][is.length];
		for (int n = 0; n < is.length; n++) {
			for (int i = 0; i < this.numRows(); i++) {
				d[i][n] = this.getDm().get(i, is[n]);
			}
		}
		return new Double2D(d);
	}

	public Integer2D toInteger2D() {
		int[][] ints = new int[this.numRows()][this.numCols()];
		for (int i = 0; i < this.numCols(); i++) {
			for (int j = 0; j < this.numRows(); j++) {
				ints[j][i] = (int) this.getDm().get(j, i);
			}
		}
		return new Integer2D(ints);
	}

	@Override
	public int length1D() {
		return this.numCols() * this.numRows();
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

	public Double2D plus(Integer2D other) {
		double[][] td = this.value();
		double[][] sd = new double[other.numRows()][other.numCols()];
		int[][] oi = other.value();
		for (int i=0; i<other.numRows(); i++) {
			for (int j=0; j<other.numCols(); j++) {
				sd[i][j] = td[i][j] + oi[i][j];
			}
		}
		return new Double2D(sd);
	}
	
	public Double2D standardized() {
		Double1D mean = this.mean();
		Double2D centred = this.minus(mean);
		double[][] scale = centred.transpose().mult(centred).diag().inverse().value();
		for (int i=0; i<scale.length; i++) {
			scale[i][i] = Math.sqrt(scale[i][i]);
		}
		Double2D stdized = centred.mult(new Double2D(scale));
		return stdized;
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

}
