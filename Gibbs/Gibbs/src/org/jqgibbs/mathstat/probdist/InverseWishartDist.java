package org.jqgibbs.mathstat.probdist;

import java.util.ArrayList;
import java.util.List;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Numeric;

import cern.jet.stat.Gamma;

public class InverseWishartDist extends ProbDistInitializeDirectly<Double2D> {

	private WishartDist wishartDist;
	private List<String> parmNames;
	private List<ProbDistParmCheck[]> parmCheck;
	private List<Class<? extends Numeric<?>>> parmClasses;

	public InverseWishartDist(Numeric<?>... parms) throws ProbDistParmException {
		super(parms);
	}

	private WishartDist getWishartDist() {
		return this.wishartDist;
	}

	@Override
	protected Double2D genVariate() throws ProbDistParmException {
		assert this.initialized;
		Double2D W = this.getWishartDist().variate();
		return W.inverse();
	}

	@Override
	protected void installParmChecks() {
		// Names
		this.parmNames = new ArrayList<String>(2);
		this.parmNames.add("Psi");
		this.parmNames.add("K");
		// Checks
		this.parmCheck = new ArrayList<ProbDistParmCheck[]>(2);
		this.parmCheck.add(null); // Don't bother testing; Wishart will do it
		this.parmCheck.add(null); // for us
		// Classes
		this.parmClasses = new ArrayList<Class<? extends Numeric<?>>>(2);
		this.parmClasses.add(Double2D.class);
		this.parmClasses.add(Double0D.class);
	}

	@Override
	protected void setUpFromParms() throws ProbDistParmException {
		if (this.getWishartDist() == null) {
			Double2D psiInv = this.getPsi().inverse();
			Double0D K = this.getK();
			this.wishartDist = new WishartDist(psiInv, K);
		} else {
			this.getWishartDist().initializeParms(this.getPsi().inverse(),
					this.getK());
		}
	}

	private Double2D getPsi() {
		return (Double2D) this.parms[0];
	}

	private Double0D getK() {
		return (Double0D) this.parms[1];
	}

	@Override
	protected List<ProbDistParmCheck[]> getParmCheck() {
		return this.parmCheck;
	}

	@Override
	protected List<Class<? extends Numeric<?>>> getParmClasses() {
		return this.parmClasses;
	}

	@Override
	protected List<String> getParmNames() {
		return this.parmNames;
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