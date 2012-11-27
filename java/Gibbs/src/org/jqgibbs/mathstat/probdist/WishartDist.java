package org.jqgibbs.mathstat.probdist;

import java.util.LinkedList;
import java.util.List;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Numeric;

public class WishartDist extends ProbDist<Double2D> {

	private static Double0D D0 = new Double0D(0);
	private static Double0D D1 = new Double0D(1);
	
	private Double2D Psi;
	private Double0D K;

	private GammaDist[] gammaDists;
	private NormalDist normalDist;
	
	protected void checkInitialized(Numeric... parms) {
		if(parms.length == 0) return;
		this.Psi = (Double2D)parms[0];
		this.K = (Double0D)parms[1];
		setUpFromParms();
	}
	
	public WishartDist(Double2D Psi, Double0D K, boolean checkParms) throws ProbDistParmException {
		this.Psi = Psi;
		this.K = K;
		if(checkParms) {
			checkParms();
		}
		this.setUpFromParms();
	}
	
	public WishartDist(Double2D Psi, Double0D K) throws ProbDistParmException {
		this(Psi, K, CHECK_PARMS);
	}
	
	public void initializeParms(Double2D Psi, Double0D K) {
		this.Psi = Psi;
		this.K = K;
	}
	
	public void checkParms() throws ProbDistParmException {
		if(!Psi.square()) {
			throw new ProbDistParmException("Expected square matrix for Psi");
		}
		if(K.value() <= ((double) getDims() - 1.0)) {
			String error = "K: Must be at > number of dimensions - 1 ("
				+ String.valueOf(WishartDist.this.getDims() - 1) + ")";
			throw new ProbDistParmException(error);
		}
	}
	
	private void setUpFromParms() {
		try {
			if (this.getGammaDists() == null) {
				this.gammaDists = new GammaDist[this.getDims()];
				for (int i=0; i<this.getDims(); i++) {
					this.gammaDists[i] = new GammaDist(new Double0D((this.getK().value()-i)/2), new Double0D(0.5));
				}
			} else {
				for (int i=0; i<this.getDims(); i++) {
					this.gammaDists[i] = new GammaDist(new Double0D((this.getK().value()-i)/2), new Double0D(0.5));
				}				
			}
			
			if (this.getNormalDist() == null) {
				this.normalDist = new NormalDist(WishartDist.D0, WishartDist.D1);
			} else {
				this.getNormalDist().initializeParms(WishartDist.D0, WishartDist.D1);
			}
		} catch (ProbDistParmException e) {
			throw new RuntimeException("Error initializing Wishart generator", e);
		}
	}

	private Double2D getPsi() {
		return Psi;
	}

	private Double0D getK() {
		return K;
	}

	private int getDims() {
		return this.getPsi().numCols();
	}

	private GammaDist[] getGammaDists() {
		return gammaDists;
	}

	private NormalDist getNormalDist() {
		return normalDist;
	}

	/*
	 * Modified Bartlett decomposition (Ku and Bloomfield).
	 * Originally derived from MCMCpack 1.0-6 for R. Copyright 2009
	 * Andrew D. Martin <admartin@wustl.edu>, Kevin M. Quinn
	 * <kevin_quinn@harvard.edu>, Jong Hee Park <jhp@uchicago.edu>.
	 */
	protected Double2D genVariate() throws ProbDistParmException {
		// Bartlett decomposition: step 1: Cholesky decomposition of scale
		// matrix
		Double2D L;
		try {
			L = this.getPsi().cholesky();
		} catch (UnsupportedOperationException e) {
			throw new ProbDistParmException("Invalid psi parameter", e);
		}
		Double2D A = null;
		boolean aOk = false;
		while (!aOk) {
			// Bartlett decomposition: step 2: construct Z, lower-triangular with
			// the diagonal sqrt of independent gamma
			// with dof v...(v-d+1), rate 0.5, the below-diagonals independent N(0,1)

			// Indices of a d x d matrix that will be sqrt gamma (diagonal), then indices
			// that will be normal
			// (below-diagonals)
			List<int[]> indices = new LinkedList<int[]>();
			for (int i = 0; i < this.getDims(); i++) {
				indices.add(new int[] { i, i });
			} // Sqrt gamma
			for (int i = 0; i < (this.getDims() - 1); i++) {
				for (int j = (i + 1); j < this.getDims(); j++) {
					indices.add(new int[] { i, j });
				} // Normal
			}
			// Sqrt gamma variates, then normal variates
			double[] v = new double[this.getDims()];
			for (int i = 0; i < this.getDims(); i++) {
				v[i] = this.getGammaDists()[i].variate().sqrt().value();
			}
			Double1D values = new Double1D(v);
			values.addAll(this.getNormalDist().variatesIID((this.getDims() - 1) * this.getDims() / 2));
			Double2D Z = new Double2D(this.getDims(), this.getDims(), indices, values);
			// Bartlett decomposition: step 3: return LZZ'L'
			Double2D LZ = L.mult(Z);
			A = LZ.mult(LZ.transpose());
			aOk = A.isWellConditioned();
			if (!aOk) {
				System.err.println("Warning: threw out ill-conditioned Wishart sample");
			}
		}
		return A;
	}
	
	@Override
	protected double getDensity(Double2D pt) throws ProbDistParmException {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}
}
