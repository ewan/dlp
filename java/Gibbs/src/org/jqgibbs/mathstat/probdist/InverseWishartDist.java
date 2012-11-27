package org.jqgibbs.mathstat.probdist;

import java.util.List;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Double3D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.stat.Gamma;

public class InverseWishartDist extends ProbDist<Double2D> {

	private WishartDist wishartDist;
	private Double2D Psi;
	private Double0D K;
	
	public static Double3D variates(List<Double2D> postPsi, List<Double0D> postKappa) throws ProbDistParmException {
		Double3D sequence = new Double3D();
		sequence.add(new InverseWishartDist(postPsi.get(0), postKappa.get(0)).variate());
		
		for(int i = 1; i < postPsi.size(); i++) {
			sequence.add(new InverseWishartDist(postPsi.get(i), postKappa.get(i), false).variate());
		}
		return sequence;
	}

	public InverseWishartDist(Double2D Psi, Double0D K, boolean checkParms) throws ProbDistParmException {
		this.Psi = Psi;
		this.K = K;
		setUpFromParms(checkParms);
	}
	
	public InverseWishartDist(Double2D Psi, Double0D K) throws ProbDistParmException {
		this(Psi, K, CHECK_PARMS);
	}

	private WishartDist getWishartDist() {
		return this.wishartDist;
	}
	
	protected void checkInitialized(Numeric... parms) {
		if(parms.length == 0) return;
		Double2D Psi = (Double2D)parms[0];
		Double0D K = (Double0D)parms[1];
		this.Psi = Psi;
		this.K = K;
		try {
			setUpFromParms(false);
		} catch (ProbDistParmException e) {
			// Do nothing (should never occur)
		}
	}

	@Override
	protected Double2D genVariate() throws ProbDistParmException {
		assert this.initialized;
		Double2D W = this.getWishartDist().variate();
		return W.inverse();
	}

	private void setUpFromParms(boolean checkParms) throws ProbDistParmException {
		if (this.getWishartDist() == null) {
			Double2D psiInv = this.getPsi().inverse();
			Double0D K = this.getK();
			this.wishartDist = new WishartDist(psiInv, K, checkParms);
		} else {
			this.getWishartDist().initializeParms(this.getPsi().inverse(),
					this.getK());
		}
	}

	private Double2D getPsi() {
		return this.Psi;
	}

	private Double0D getK() {
		return this.K;
	}

	@Override
	protected double getDensity(Double2D pt) {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}	
	
	public double getLogDensity(Double2D pt) {
		int d = this.getPsi().numCols();
		double k = this.getK().value();
		double ld = -k*d*Math.log(2)/2;
 		double detPt = Math.pow(pt.det(), (k+d+1)/2);
		double detPsi = Math.pow(this.getPsi().det(), k/2);
		double gamK2 = 0; 
		double a = k/2;
		for (int i=d; i>1; i--) {
			gamK2 += Math.log(Math.pow(Math.PI, ((double)i-1)/2));
			gamK2 += Math.log(Gamma.gamma(a));
			a -= 0.5;
		}
		gamK2 += Math.log(Gamma.gamma(a));
		ld -= Math.log(detPt);
		ld += Math.log(detPsi);
		ld -= gamK2;
		Double2D PsiPtInv = this.getPsi().mult(pt.inverse());
		ld -= 0.5*PsiPtInv.trace().value();
		return ld;
	}
}