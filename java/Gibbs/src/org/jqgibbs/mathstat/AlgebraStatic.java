package org.jqgibbs.mathstat;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.jet.math.Functions;

public final class AlgebraStatic extends cern.colt.matrix.linalg.Algebra {

	private AlgebraStatic() {}
	
	private static cern.colt.matrix.linalg.Algebra a = new cern.colt.matrix.linalg.Algebra();

	private static final long serialVersionUID = 7354941851041130909L;

	public static DoubleMatrix2D inverseStatic(DoubleMatrix2D A) {
		// FIXME - This seems to crash unexpectedly sometimes
		return AlgebraStatic.a.inverse(A);
	}

	public static DoubleMatrix2D multStatic(DoubleMatrix2D A, DoubleMatrix2D B) {
		return AlgebraStatic.a.mult(A, B);
	}

	public static DoubleMatrix2D transposeStatic(DoubleMatrix2D A) {
		return AlgebraStatic.a.transpose(A);
	}

	public static DoubleMatrix1D multStatic(DoubleMatrix2D A, DoubleMatrix1D B) {
		return AlgebraStatic.a.mult(A, B);
	}

	public static DoubleMatrix2D multStatic(final double k, DoubleMatrix2D colt) {
		DoubleMatrix2D r = colt.copy().assign(Functions.mult(k));
		return (r);
	}

	public static DoubleMatrix2D plusStatic(DoubleMatrix2D colt,
			DoubleMatrix2D colt2) {
		return colt.copy().assign(colt2, Functions.plus);
	}

	public static DoubleMatrix1D plusStatic(DoubleMatrix1D colt,
			DoubleMatrix1D colt2) {
		return colt.copy().assign(colt2, Functions.plus);
	}

	public static DoubleMatrix2D multOuterStatic(DoubleMatrix1D dm,
			DoubleMatrix1D dm2) {
		DoubleMatrix2D r = new DenseDoubleMatrix2D(dm.size(), dm2.size());
		AlgebraStatic.a.multOuter(dm, dm2, r);
		return r;
	}

	public static DoubleMatrix2D minusStatic(DoubleMatrix2D colt,
			DoubleMatrix1D v) {
		DoubleMatrix2D copy = colt.copy();
		DoubleMatrix1D r;
		for (int ri = 0; ri < copy.rows(); ri++) {
			r = copy.viewRow(ri);
			r.assign(v, Functions.minus);
		}
		return copy;
	}

	public static DoubleMatrix1D multStatic(double k, DoubleMatrix1D colt) {
		DoubleMatrix1D r = colt.copy().assign(Functions.mult(k));
		return (r);
	}

	public static DoubleMatrix1D minusStatic(DoubleMatrix1D colt1,
			DoubleMatrix1D colt2) {
		return colt1.copy().assign(colt2, Functions.minus);
	}

	public static DoubleMatrix2D minusStatic(DoubleMatrix2D colt1,
			DoubleMatrix2D colt2) {
		return colt1.copy().assign(colt2, Functions.minus);
	}

	public static DoubleMatrix1D multStatic(DoubleMatrix1D B, DoubleMatrix2D A) {
		return AlgebraStatic.a.mult(AlgebraStatic.a.transpose(A), B);
	}

	public static double multStatic(DoubleMatrix1D colt1, DoubleMatrix1D colt2) {
		return AlgebraStatic.a.mult(colt1, colt2);
	}

	public static double detStatic(DoubleMatrix2D colt) {
		return AlgebraStatic.a.det(colt);
	}

	public static DoubleMatrix1D minusStatic(DoubleMatrix1D colt, double k) {
		DoubleMatrix1D r = colt.copy().assign(Functions.minus(k));
		return (r);
	}

	public static DoubleMatrix1D plusStatic(DoubleMatrix1D colt, double k) {
		DoubleMatrix1D r = colt.copy().assign(Functions.plus(k));
		return (r);
	}

	public static DoubleMatrix2D divideStatic(double k, DoubleMatrix2D colt) {
		DoubleMatrix2D r = colt.copy().assign(Functions.div(k));
		return (r);
	}

	public static DoubleMatrix1D divideStatic(double k, DoubleMatrix1D colt) {
		DoubleMatrix1D r = colt.copy().assign(Functions.div(k));
		return (r);
	}

	public static void assignBlockStatic(DoubleMatrix2D receiver,
			int top, int left, DoubleMatrix2D sender) {
		// FIXME - no dim checks		
		for (int i=top; i<(top+sender.rows()); i++) {
			for (int j=left; j<(left+sender.columns()); j++) {
				receiver.setQuick(i, j, sender.get(i-top, j-left));
			}
		}
	}

	public static DoubleMatrix2D kronStatic(DoubleMatrix2D left,
			DoubleMatrix2D right) {
		int k = left.rows();
		int l = left.columns();
		int m = right.rows();
		int n = right.columns();
		DoubleMatrix2D prod = new DenseDoubleMatrix2D(k * m, l * n);
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < l; j++) {
				DoubleMatrix2D lijr = AlgebraStatic.multStatic(left.get(i, j),
						right);
				AlgebraStatic.assignBlockStatic(prod, i * m, j * n, lijr);
			}
		}
		return prod;
	}

	public static DoubleMatrix2D unvec(DoubleMatrix1D dm, int h) {
		// FIXME - no dim check
		int rows = h;
		int cols = dm.size()/h;
		DoubleMatrix2D mat = new DenseDoubleMatrix2D(rows, cols);
		int k = 0;
		for (int j=0; j<cols; j++) {
			for (int i=0; i<rows; i++) {
				mat.setQuick(i, j, dm.get(k));
				k++;
			}
		}
		return mat;
	}
}
