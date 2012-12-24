package org.jqgibbs.models;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.jqgibbs.ChainLink;
import org.jqgibbs.Model;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Double3D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Integer1D;
import org.jqgibbs.mathstat.Integer2D;
import org.jqgibbs.mathstat.ListSequence;
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.mathstat.RandomVar;
import org.jqgibbs.mathstat.probdist.BetaDist;
import org.jqgibbs.mathstat.probdist.CategoricalDist;
import org.jqgibbs.mathstat.probdist.GammaDist;
import org.jqgibbs.mathstat.probdist.InverseWishartDist;
import org.jqgibbs.mathstat.probdist.MVNormalDist;
import org.jqgibbs.mathstat.probdist.MatrixNormalDist;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistInitializeByChain;
import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;
import org.jqgibbs.mathstat.probdist.SeqInverseWishartDist;
import org.jqgibbs.mathstat.probdist.SeqMVNormalDist;
import org.jqgibbs.mathstat.probdist.SeqMatrixNormalDist;
import cern.colt.matrix.linalg.SingularValueDecomposition;

import cern.jet.random.Uniform;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.stat.Gamma;

public class FLGFWModel extends Model {
	
	private static double gammaD(double x, int d) { // FIXME
		if (d == 1) {
			return Gamma.gamma(x);
		} else { // FIXME check d0
			return Math.pow(Math.PI, (d - 1) / 2) * FLGFWModel.gammaD(x, d - 1)
					* Gamma.gamma(x + (1 - d) / 2);
		}
	}

	private static Uniform unifGen = new Uniform(
			new MersenneTwister(new Date()));

	class PostGammaDist extends ProbDistInitializeByChain<Integer1D> {
		public PostGammaDist(Numeric<?>... fixed)
				throws ProbDistParmException {
			super(fixed);
		}
		
		private CategoricalDist catDist;
		private Double1D postProb;
		
		protected Integer0D catVariate() throws ProbDistParmException {
			if (this.catDist == null) {
				this.catDist = new CategoricalDist(this
						.getPostProb());
				return this.catDist.variate();
			} else {
				return this.catDist.variate(this.getPostProb());
			}
		}			
		
		private Double1D getPostProb() {
			return this.postProb;
		}

		protected void setPostProb(Double1D prob) {
			this.postProb = prob;
		}

		protected Integer1D getGamma() {
			return (Integer1D) this.chainParms[1];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[0];
		}

		protected Double3D getA() {
			return (Double3D) this.chainParms[2];
		}

		protected Integer2D getB() {
			return (Integer2D) this.chainParms[3];
		}

		protected Double2D getA0() {
			return (Double2D) this.chainParms[4];
		}
		
		protected Double3D getSg() {
			return (Double3D) this.chainParms[5];
		}

		protected Double2D getM() {
			return (Double2D) this.chainParms[6];
		}

		protected Double1D getM0() {
			return (Double1D) this.chainParms[7];
		}
		
		protected Double0D getRho() {
			return (Double0D) this.chainParms[8];
		}
		
		protected Double0D getRho0() {
			return (Double0D) this.chainParms[9];
		}
		
		protected Double0D getP() {
			return (Double0D) this.fixedParms[0];
		}

		protected Double0D getLambda() {
			return (Double0D) this.fixedParms[1];
		}

		protected Double0D getBeta() {
			return (Double0D) this.fixedParms[2];
		}

		protected Double2D getS() {
			return (Double2D) this.fixedParms[3];
		}

		protected Double0D getLambda0() {
			return (Double0D) this.fixedParms[4];
		}

		protected Double0D getBeta0() {
			return (Double0D) this.fixedParms[5];
		}
		
		protected Double2D getPsi() {
			return (Double2D) this.fixedParms[6];
		}
		
		protected Double0D getKappa() {
			return (Double0D) this.fixedParms[7];
		}

		protected Double2D getS0() {
			return (Double2D) this.fixedParms[8];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(6);
			this.chainParmNames.add("Z");
			this.chainParmNames.add("Gamma");
			this.chainParmNames.add("A");
			this.chainParmNames.add("B");
			this.chainParmNames.add("A0");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("M");
			this.chainParmNames.add("M0");
			this.chainParmNames.add("rho");
			this.chainParmNames.add("rho0");
			this.fixedParmNames = new ArrayList<String>(3);
			this.fixedParmNames.add("p");
			this.fixedParmNames.add("lambda");
			this.fixedParmNames.add("beta");
			this.fixedParmNames.add("S");
			this.fixedParmNames.add("lambda0");
			this.fixedParmNames.add("beta0");
			this.fixedParmNames.add("Psi");
			this.fixedParmNames.add("kappa");
			this.fixedParmNames.add("S0");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(6);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(3);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					6);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(3);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double2D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			return; // FIXME??
		}
		
		
		protected void resampleA() throws ProbDistParmException {
			int[] gamma =  this.getGamma().value();
			int h = gamma.length;
			Integer1D active = this.getZ().items();
			int d = FLGFWModel.this.dims;
			double[] z = new double[d];
			for (int i=0; i<d; i++) {
				z[0] = 0;
			}
			Double1D zero = new Double1D(z);
			Double2D[] aNewA = (Double2D[]) this.getA().toArray();
			for (int i=0; i<h; i++) {
				if (gamma[i] == 0) {
					for (Integer0D k : active) {
						aNewA[k.value()].set(i, zero);
					}
				}
			}
			int[] gamma1 = this.getGamma().which(new Integer0D(1)).value();
			int h1 = gamma1.length;
			MatrixNormalDist mnDist = new MatrixNormalDist();
			for (Integer0D k : active) {
				Integer1D zK = this.getZ().which(k);
				Double2D xK = this.getSamplerData().getAll(zK.value());
				Integer2D bPlusT = this.getB().getAll(zK.value()).transpose().getAll(gamma1);
				int[] ones = new int[xK.size()];
				for (int j=0; j<ones.length; j++) {
					ones[j] = 1;
				}
				bPlusT.add(new Integer1D(ones)); 
				double[][] poop = new double[h1+1][h1+1]; 
				if (h1 > 0) {
					for (int i=0; i<h1; i++) {
						poop[i][i] = this.getRho().value();
					}
				}
				poop[h1][h1] = this.getRho0().value();
				Double2D poop2D = new Double2D(poop);	
				Double2D mPlus = (Double2D) this.getM().transpose().transpose().getAll(gamma1);
				mPlus.add(this.getM0());
				Double2D s = this.getSg().get(k.value());
				Double2D o = poop2D.plus(bPlusT.mult(bPlusT.transpose())).inverse();
				Double2D m = o.mult(bPlusT.mult(xK).plus(poop2D.mult(mPlus)));
				Double2D aRed = mnDist.variate(m,s,o);
				for (int i=0; i<h1; i++) {
					aNewA[k.value()].set(gamma1[i], aRed.get(i));
				}
				this.getA0().set(k.value(), aRed.get(h1)); // !!! FIXME 
			}
			for (Integer0D k : active) {
				this.getA().set(k.value(), aNewA[k.value()]);
			}
		}
			
		protected void resampleM() throws ProbDistParmException {
			Integer1D active = this.getZ().items();
			int[] gamma = this.getGamma().value();
			int h = gamma.length;
			Double1D zero = new Double1D(new double[FLGFWModel.this.dims]);
			for (int i=0; i<h; i++) {
				if (gamma[i] == 0) {
					this.getM().set(i, zero);
				}
			}			
			MVNormalDist mvnDist = new MVNormalDist();
			Double2D sInv = this.getS0().inverse();
			Double2D sgA = null;
			for (Integer0D k : active) {
				Double2D sgRho = this.getSg().get(k.value()).inverse().mult(this.getRho());
				sInv = sInv.plus(sgRho);
				if (sgA == null) {
					sgA = sgRho.mult(this.getA().get(k.value()).transpose());
				} else {
					sgA = sgA.plus(sgRho.mult(this.getA().get(k.value()).transpose()));
				}
			}	
			sgA = sgA.transpose();
			for (int i=0; i<h; i++) {
				if (gamma[i] == 1) {
					this.getM().set(i, mvnDist.variate(sgA.get(i), sInv.inverse()));
				}
			}			
			
		}
		
		protected void resampleM0() throws ProbDistParmException {
			Integer1D active = this.getZ().items();
			MVNormalDist mvnDist = new MVNormalDist();
			Double2D sInv = this.getS().inverse();
			Double1D sgA0 = null;
			for (Integer0D k : active) {
				Double2D sgRho = this.getSg().get(k.value()).inverse().mult(this.getRho0());
				sInv = sInv.plus(sgRho);
				if (sgA0 == null) {
					sgA0 = sgRho.mult(this.getA0().get(k.value()));
				} else {
					sgA0 = sgA0.plus(sgRho.mult(this.getA0().get(k.value())));
				}
			}
			Double1D m0New = mvnDist.variate(sgA0, sInv.inverse());	
			for (int i=0; i<FLGFWModel.this.dims; i++) {
				this.getM0().set(i, m0New.get(i));		
			}
		}
		
		protected void resampleSg() throws ProbDistParmException {
			Integer1D active = this.getZ().items();
			int h1 = this.getGamma().sum().value();
			int[] gamma1 = this.getGamma().which(new Integer0D(1)).value();
			Double2D psi;
			Double0D kappa;
			InverseWishartDist iwDist = null;
			for (Integer0D k : active) {
				// set df parameter
				Integer1D zK = this.getZ().which(k);
				Double2D xK = this.getSamplerData().getAll(zK.value());
				kappa = this.getKappa().plus(xK.size());
				// set scale matrix 
				psi = this.getPsi().plus(xK.transpose().mult(xK));
				Double2D aPlus = (Double2D) this.getA().get(k.value()).transpose().transpose();
				aPlus.add(this.getA0().get(k.value()));
				Integer2D bPlusT = this.getB().getAll(zK.value()).transpose().getAll(gamma1);
				int[] ones = new int[xK.size()];
				for (int j=0; j<ones.length; j++) {
					ones[j] = 1;
				}
				bPlusT.add(new Integer1D(ones)); 
				double[][] poop = new double[h1+1][h1+1]; 
				if (h1 > 0) {
					for (int i=0; i<h1; i++) {
						poop[i][i] = this.getRho().value();
					}
				}
				poop[h1][h1] = this.getRho0().value();
				Double2D poop2D = new Double2D(poop);
				Double2D mPlus = (Double2D) this.getM().transpose().transpose().getAll(gamma1);
				mPlus.add(this.getM0());
				psi = psi.plus(mPlus.transpose().mult(poop2D).mult(mPlus));
				Double2D mX = poop2D.mult(mPlus).plus(bPlusT.mult(xK));
				psi = psi.minus(mX.transpose().mult(poop2D.plus(bPlusT.mult(bPlusT.transpose())).inverse()).mult(mX));
				if (iwDist == null) {
					iwDist = new InverseWishartDist(psi, kappa);
					this.getSg().set(k.value(), iwDist.variateFast());
				} else {
					this.getSg().set(k.value(), iwDist.variateFast(psi, kappa));
				}
			}
		}
		
		protected void resampleRho() throws ProbDistParmException {
			int d = this.getM().numCols();
			Integer1D active = this.getZ().items();
			int K = active.size();
			int[] gamma1 = this.getGamma().which(new Integer0D(1)).value();
			int h1 = gamma1.length;
			Double0D l = this.getLambda().plus(d*h1*K/2);
			Double0D b = this.getBeta();
			for (Integer0D k : active) {
				for (int i=0; i<h1; i++) {
					Double1D ac = this.getA().get(k.value()).get(gamma1[i]).minus(this.getM().get(gamma1[i]));
					b = b.plus(0.5*ac.mult(this.getSg().get(k.value()).inverse()).mult(ac));
				}
			}
			GammaDist gammaDist = new GammaDist(l, b);
			Double0D v = gammaDist.variateFast();
			this.getRho().set(v.value());
		}
		
		protected void resampleRho0() throws ProbDistParmException {
			int d = this.getM0().size();
			Integer1D active = this.getZ().items();
			int K = active.size();
			Double0D l = this.getLambda0().plus(d*K/2);
			Double0D b = this.getBeta0();
			for (Integer0D k : active) {
				Double1D ac = this.getA0().get(k.value()).minus(this.getM0());
				b = b.plus(0.5*ac.mult(this.getSg().get(k.value()).inverse()).mult(ac));
			}
			GammaDist gammaDist = new GammaDist(l, b);
			Double0D v = gammaDist.variateFast();
			this.getRho0().set(v.value());
		}

		@Override
		protected Integer1D genVariate() throws ProbDistParmException {
			// get current gamma
			int[] gamma = null;
			try {
				gamma = ((Integer1D) this.getGamma().clone()).value();
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			int h = gamma.length;
			// split up data assigned to each point in the DP mixture
			int K = this.getA().size();
			Integer1D activeZ = this.getZ().items();
			int Ka = activeZ.size();
			Integer1D[] zK = new Integer1D[K];
			Double2D[] xK = new Double2D[K];
			Integer2D[] bK = new Integer2D[K];
			for (Integer0D k : activeZ) {
				zK[k.value()] = this.getZ().which(k);
				xK[k.value()] = this.getSamplerData().getAll(zK[k.value()].value());
				bK[k.value()] = this.getB().getAll(zK[k.value()].value());
			}
			//
			int d = FLGFWModel.this.dims;
			Double1D[] be = new Double1D[K];
			double[] bsum = new double[K];
			double[] bro = new double[K];
			LinkedList<Integer> samplingOrder = new LinkedList<Integer>();
			LinkedList<Integer> remainingIndices = new LinkedList<Integer>();
			for (int i=0; i<h; i++) {
				remainingIndices.addLast(i);
			}
			int whichNext;
			for (int i=1; i<h; i++) {
				whichNext = unifGen.nextIntFromTo(0, i);
				samplingOrder.addLast(remainingIndices.get(whichNext));
				remainingIndices.remove(whichNext);
			}
			samplingOrder.addLast(remainingIndices.get(0));
			for (Integer index : samplingOrder) {
				int i = index.intValue();
				// Get non-zero gammas
				int rowsRed = 0;
				int[] indicesRed = new int[h];
				for (int r = 0; r < h; r++) {
					if (!(r==i || gamma[r] == 0)) {
						indicesRed[rowsRed] = r;
						++rowsRed;
					}
				}
				// be[k] := \sum_j (x_j - a_{0,k} - \sum_{h-i} b_{h,j} * a_{h,k})
				for (Integer0D k : activeZ) {
					// \sum_j x_j
					be[k.value()] = xK[k.value()].sum();
					// \sum_j a_{0,k}
					int nK = xK[k.value()].numRows();
					be[k.value()] = be[k.value()].minus(this.getA0().get(k.value()).mult(nK));	
					// \sum_j \sum_{h-i} b_{h,j} * a_{h,k}
					if (rowsRed > 0) {
						double[][] A = this.getA().get(k.value()).value();
						int[] B = bK[k.value()].sum().value();
						double[][] redA = new double[rowsRed][d];
						double[] redB = new double[rowsRed];
						for (int s = 0; s < rowsRed; s++) {
							int r = indicesRed[s];
							for (int c=0; c<d; c++) {
								redA[s][c] = A[r][c];
							}
							redB[s] = B[r];
						}
						Double1D ba = (new Double1D(redB)).mult(new Double2D(redA));
						be[k.value()] = be[k.value()].minus(ba);
					}
				}
				// bsum[k] := \sum_j (b_{i,j})
				for (Integer0D k : activeZ) {
					double bs = 0;
					int[][] bKval = bK[k.value()].value();
					for (int j=0; j<bKval.length; j++) {
						bs += bKval[j][i];
					}
					bsum[k.value()] = bs;
				}
				// bro[k] := bsum[k]^2 + rho
				for (Integer0D k : activeZ) {
					bro[k.value()] = Math.pow(bsum[k.value()], 2) + this.getRho().value();
				}	
				// draw gamma
				// p(g=1|X)
				double logP = Math.log(this.getP().value());
				logP += ((d*Ka)/2)*Math.log(this.getRho().value());
				logP -= 0.5*Math.log(this.getS().det());
//				logP -= (d/2)*Math.log(2*Math.PI);
				logP -= 0.5*this.getM().get(i).mult(this.getS().inverse()).mult(this.getM().get(i));
				Double2D ms = this.getS().inverse();
				for (Integer0D k : activeZ) {
					double rc = this.getRho().value() - Math.pow(this.getRho().value()/bro[k.value()], 2);
					Double2D rSgInv = this.getSg().get(k.value()).inverse().mult(rc);
					ms = ms.plus(rSgInv);
				}
				ms = ms.inverse();
				logP += 0.5*Math.log(ms.det());
				for (Integer0D k : activeZ) {
					logP += 0.5*Math.log(this.getSg().get(k.value()).det());
					Double2D sgInv = this.getSg().get(k.value()).inverse();
					logP -= (d/2)*Math.log(bro[k.value()]);
//					logP -= 0.5*(1-Math.pow(bsum[k.value()]/bro[k.value()], 2))*be[k.value()].mult(sgInv).mult(be[k.value()]);
//					logP += 0.5*Math.pow(bsum[k.value()]/bro[k.value()], 2)*be[k.value()].mult(sgInv).mult(be[k.value()]);
//					logP -= 0.5*(this.getRho().value()-Math.pow(this.getRho().value()/bro[k.value()], 2))*this.getM().get(i).mult(sgInv).mult(this.getM().get(i));
//					logP -= (this.getRho().value()*bsum[k.value()]/bro[k.value()])*be[k.value()].mult(sgInv).mult(this.getM().get(i));
					logP += 0.5*Math.pow(bsum[k.value()]/bro[k.value()], 2)*be[k.value()].mult(sgInv).mult(be[k.value()]);
					Double1D bp = be[k.value()].mult(sgInv).mult(this.getRho().value() * bsum[k.value()] / Math.pow(bro[k.value()], 2));
					logP += 0.5*bp.mult(ms).mult(bp);
					
				}
				// p(g=0|X) \propto 1-p
				double logQ = Math.log(1-this.getP().value());
				// 2. exp{\sum_k be[k]^T \Sigma_k^{-1} be[k]}
//				for (Integer0D k : activeZ) {
//					Double2D sgK = this.getSg().get(k.value());	
//					int nK = xK[k.value()].size();
//					if (nK > 0) {
//						logQ -= 0.5*be[k.value()].mult(sgK.inverse()).mult(be[k.value()]);
//					}
//				}
				// Normalizing constant
				double maxP;
				if (logP > logQ) {
					maxP = logP;
				} else {
					maxP = logQ;
				}
				Double1D p = (new Double1D(logQ, logP)).minus(maxP).exp();
				this.setPostProb(p);
				int gammaNew = this.catVariate().value();
				if (gammaNew != gamma[i]) {
					gamma[i] = gammaNew;
//					this.getGamma().set(i, new Integer0D(gamma[i]));
//					// Resample A
//					this.resampleA();
//					// Resample M
//					this.resampleM();
//					// Resample M0
//					this.resampleM0();
//					// Resample Sg
//					this.resampleSg();
//					// Resample rho
//					this.resampleRho();
//					// Resample rho0
//					this.resampleRho0();
				}
			}
			// return new values of gamma
			return new Integer1D(gamma);
		}

		@Override
		protected double getDensity(Integer1D pt) throws ProbDistParmException {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}
	
	class PriorSgDist extends ProbDistInitializeByChain<Double3D> {
		private SeqInverseWishartDist iwDists; // FIXME - Shouldn't this use the
		// IID?
		private ListSequence<Double2D> repPsi;
		private ListSequence<Double0D> repKappa;

		public PriorSgDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private ListSequence<Double2D> getRepPsi() {
			if (this.repPsi == null) {
				this.repPsi = new ListSequence<Double2D>();
				Integer1D active = this.getZ().items();
				for (int i = 0; i < active.size(); i++) {
					this.repPsi.add(this.getPsi());
				}
			}
			return this.repPsi;
		}

		private ListSequence<Double0D> getRepKappa() {
			if (this.repKappa == null) {
				this.repKappa = new ListSequence<Double0D>();
				Integer1D active = this.getZ().items();
				for (int i = 0; i < active.size(); i++) {
					this.repKappa.add(this.getKappa());
				}
			}
			return this.repKappa;
		}

		private Double2D getPsi() {
			return (Double2D) this.fixedParms[0];
		}

		private Double0D getKappa() {
			return (Double0D) this.fixedParms[1];
		}

		private Integer1D getZ() {
			return (Integer1D) this.chainParms[0];
		}

		private Double3D iwVariates() throws ProbDistParmException {
			if (this.iwDists == null) {
				this.iwDists = new SeqInverseWishartDist(this.getRepPsi(),
						this.getRepKappa());
				return this.iwDists.variate();
			} else {
				return this.iwDists.variate(this.getRepKappa(),
						this.getRepKappa());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(1);
			this.chainParmNames.add("Z");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("Psi");
			this.fixedParmNames.add("kappa");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
			this.chainParmCheck.add(null); // FIXME
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					1);
			this.chainParmClasses.add(Integer1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		@Override
		protected Double3D genVariate() throws ProbDistParmException {
			return this.iwVariates();
		}

		@Override
		protected double getDensity(Double3D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PostSgDist extends ProbDistInitializeByChain<Double3D> {
		private ListSequence<Double2D> postPsi;
		private ListSequence<Double0D> postKappa;
		private SeqInverseWishartDist iwDists;

		public PostSgDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected ListSequence<Double2D> getPostPsi() {
			return this.postPsi;
		}

		protected ListSequence<Double0D> getPostKappa() {
			return this.postKappa;
		}

		protected void setPostPsi(ListSequence<Double2D> postPsi) {
			this.postPsi = postPsi;
		}

		protected void setPostKappa(ListSequence<Double0D> postKappa) {
			this.postKappa = postKappa;
		}

		protected Double2D getPsi() {
			return (Double2D) this.fixedParms[0];
		}

		protected Double0D getKappa() {
			return (Double0D) this.fixedParms[1];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[0];
		}

		protected Integer2D getB() {
			return (Integer2D) this.chainParms[1];
		}
		
		protected Double3D getA() {
			return (Double3D) this.chainParms[2];
		}
		
		protected Integer1D getGamma() {
			return (Integer1D) this.chainParms[3];
		}
		
		protected Double0D getRho() {
			return (Double0D) this.chainParms[4];
		}

		protected Double0D getRho0() {
			return (Double0D) this.chainParms[5];
		}
		
		protected Double2D getA0() {
			return (Double2D) this.chainParms[6];
		}
	
		protected Double2D getM() {
			return (Double2D) this.chainParms[7];
		}
		
		protected Double1D getM0() {
			return (Double1D) this.chainParms[8];
		}
		
		private Double3D iwVariates() throws ProbDistParmException {
			if (this.iwDists == null) {
				this.iwDists = new SeqInverseWishartDist(this.getPostPsi(),
						this.getPostKappa());
				return this.iwDists.variate();
			} else {
				return this.iwDists.variate(this.getPostPsi(), this.getPostKappa());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(4);
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.chainParmNames.add("A");
			this.chainParmNames.add("Gamma");
			this.chainParmNames.add("rho");
			this.chainParmNames.add("rho0");
			this.chainParmNames.add("A0");
			this.chainParmNames.add("M");
			this.chainParmNames.add("M0");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("Psi");
			this.fixedParmNames.add("kappa");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					5);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			Integer1D active = this.getZ().items();
			this.setPostPsi(new ListSequence<Double2D>());
			this.setPostKappa(new ListSequence<Double0D>());
			int h1 = this.getGamma().sum().value();
			int[] gamma1 = this.getGamma().which(new Integer0D(1)).value();
			Double2D psi = this.getPsi();
			for (Integer0D k : active) {
				// set df parameter
				Integer1D zK = this.getZ().which(k);
				Double2D xK = this.getSamplerData().getAll(zK.value());
				this.getPostKappa().add(this.getKappa().plus(xK.size()));
				// set scale matrix 
				psi = psi.plus(xK.transpose().mult(xK));
				Double2D aPlus = (Double2D) this.getA().get(k.value()).transpose().transpose();
				aPlus.add(this.getA0().get(k.value()));
				Integer2D bPlusT = this.getB().getAll(zK.value()).transpose().getAll(gamma1);
				int[] ones = new int[xK.size()];
				for (int j=0; j<ones.length; j++) {
					ones[j] = 1;
				}
				bPlusT.add(new Integer1D(ones)); 
				double[][] poop = new double[h1+1][h1+1]; 
				if (h1 > 0) {
					for (int i=0; i<h1; i++) {
						poop[i][i] = this.getRho().value();
					}
				}
				poop[h1][h1] = this.getRho0().value();
				Double2D poop2D = new Double2D(poop);
				Double2D mPlus = (Double2D) this.getM().transpose().transpose().getAll(gamma1);
				mPlus.add(this.getM0());
				psi = psi.plus(mPlus.transpose().mult(poop2D).mult(mPlus));
				Double2D mX = poop2D.mult(mPlus).plus(bPlusT.mult(xK));
				psi = psi.minus(mX.transpose().mult(poop2D.plus(bPlusT.mult(bPlusT.transpose())).inverse()).mult(mX));
				this.getPostPsi().add(psi);
			}
		}

		@Override
		protected Double3D genVariate() throws ProbDistParmException {
			int d = FLGFWModel.this.dims;
			double[][][] v = new double[this.getA().size()][d][d];
			Double3D activeUpdates = this.iwVariates();
			int i = 0;
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], d);
				i++;
			}
			return new Double3D(v); // I have no idea why I copy this array
		}

		@Override
		protected double getDensity(Double3D pt) throws ProbDistParmException {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PriorADist extends ProbDistInitializeByChain<Double3D> {
		private SeqMatrixNormalDist matnDists;

		public PriorADist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private Double2D getOmega() {
			return (Double2D) this.chainParms[0];
		}

		private Double3D getActiveSg() {
			Double3D sg = (Double3D) this.chainParms[1];
			Integer1D active = this.getZ().items();
			return sg.getAll(active.value());
		}

		private Double3D getRepM() {
			Double2D M = (Double2D) this.chainParms[2];
			Integer1D active = this.getZ().items();
			Double3D repM = new Double3D(active.size(), this.getB().numCols(),
					FLGFWModel.this.dims);
			for (int i = 0; i < active.size(); i++) {
				repM.add(M);
			}
			return repM;
		}

		private Integer1D getZ() {
			return (Integer1D) this.chainParms[3];
		}

		private Integer2D getB() {
			return (Integer2D) this.chainParms[4];
		}

		private Double3D matnVariates() throws ProbDistParmException {
			if (this.matnDists == null) {
				this.matnDists = new SeqMatrixNormalDist(this.getRepM(),
						this.getActiveSg(), this.getOmega());
				return this.matnDists.variate();
			} else {
				return this.matnDists.variate(this.getRepM(),
						this.getActiveSg(), this.getOmega());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(4);
			this.chainParmNames.add("Omega");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("M");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.fixedParmNames = new ArrayList<String>(0);
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					4);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					0);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		@Override
		protected Double3D genVariate() throws ProbDistParmException {
			return this.matnVariates();
		}

		@Override
		protected double getDensity(Double3D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PostADist extends ProbDistInitializeByChain<Double3D> {
		private ListSequence<Double2D> postMs;
		private ListSequence<Double2D> postOmegas;

		public PostADist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected ListSequence<Double2D> getPostM() {
			return this.postMs;
		}

		protected ListSequence<Double2D> getPostOmegas() {
			return this.postOmegas;
		}

		protected void setPostM(ListSequence<Double2D> postMs) {
			this.postMs = postMs;
		}

		protected void setPostOmegas(ListSequence<Double2D> postOmegas) {
			this.postOmegas = postOmegas;
		}

		protected Double3D getA() {
			return (Double3D) this.chainParms[0];
		}

		protected Double3D getSg() {
			return (Double3D) this.chainParms[1];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[2];
		}

		protected Integer2D getB() {
			return (Integer2D) this.chainParms[3];
		}
		
		protected Double2D getM() {
			return (Double2D) this.chainParms[4];
		}
		
		protected Integer1D getGamma() {
			return (Integer1D) this.chainParms[5];
		}
			
		protected Double2D getA0() {
			return (Double2D) this.chainParms[6];
		}
		
		protected Double0D getRho() {
			return (Double0D) this.chainParms[7];
		}
		
		protected Double0D getRho0() {
			return (Double0D) this.chainParms[8];
		}
		
		protected Double1D getM0() {
			return (Double1D) this.chainParms[9];
		}
			
		protected Double0D getLambda() {
			return (Double0D) this.fixedParms[0];
		}

		protected Double0D getBeta() {
			return (Double0D) this.fixedParms[1];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(6);
			this.chainParmNames.add("A");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.chainParmNames.add("M");
			this.chainParmNames.add("Gamma");
			this.chainParmNames.add("A0");
			this.chainParmNames.add("rho");
			this.chainParmNames.add("rho0");
			this.chainParmNames.add("M0");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("lambda");
			this.fixedParmNames.add("beta");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(6);
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					6);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		@Override
		protected Double3D genVariate() throws ProbDistParmException {
			int[] gamma =  this.getGamma().value();
			int h = gamma.length;
			Integer1D active = this.getZ().items();
			int d = FLGFWModel.this.dims;
			double[] z = new double[d];
			for (int i=0; i<d; i++) {
				z[0] = 0;
			}
			Double1D zero = new Double1D(z);
			Double3D aNew = null;
			try {
				aNew = (Double3D) this.getA().clone();
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			Double2D[] aNewA = (Double2D[]) aNew.toArray();
			for (int i=0; i<h; i++) {
				if (gamma[i] == 0) {
					for (Integer0D k : active) {
						aNewA[k.value()].set(i, zero);
					}
				}
			}
			int[] gamma1 = this.getGamma().which(new Integer0D(1)).value();
			int h1 = gamma1.length;
			MatrixNormalDist mnDist = new MatrixNormalDist();
			for (Integer0D k : active) {
				Integer1D zK = this.getZ().which(k);
				Double2D xK = this.getSamplerData().getAll(zK.value());
				Integer2D bPlusT = this.getB().getAll(zK.value()).transpose().getAll(gamma1);
				int[] ones = new int[xK.size()];
				for (int j=0; j<ones.length; j++) {
					ones[j] = 1;
				}
				bPlusT.add(new Integer1D(ones)); 
				double[][] poop = new double[h1+1][h1+1]; 
				if (h1 > 0) {
					for (int i=0; i<h1; i++) {
						poop[i][i] = this.getRho().value();
					}
				}
				poop[h1][h1] = this.getRho0().value();
				Double2D poop2D = new Double2D(poop);	
				Double2D mPlus = (Double2D) this.getM().transpose().transpose().getAll(gamma1);
				mPlus.add(this.getM0());
				Double2D s = this.getSg().get(k.value());
				Double2D o = poop2D.plus(bPlusT.mult(bPlusT.transpose())).inverse();
				Double2D m = o.mult(bPlusT.mult(xK).plus(poop2D.mult(mPlus)));
				Double2D aRed = mnDist.variate(m,s,o);
				for (int i=0; i<h1; i++) {
					aNewA[k.value()].set(gamma1[i], aRed.get(i));
				}
				this.getA0().set(k.value(), aRed.get(h1)); // !!! FIXME 
			}
			for (Integer0D k : active) {
				aNew.set(k.value(), aNewA[k.value()]);
			}
			return aNew;
		}
			
		@Override
		protected double getDensity(Double3D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}
	
	
	class PriorA0Dist extends ProbDistInitializeByChain<Double2D> {
		public PriorA0Dist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}
		
		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(4);
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("M");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.fixedParmNames = new ArrayList<String>(0);
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					4);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					0);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		@Override
		protected Double2D genVariate() throws ProbDistParmException {
			return new Double2D(); //FIXME
		}

		@Override
		protected double getDensity(Double2D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}
	

	class CRPDist extends ProbDistInitializeByChain<Integer1D> {
		public CRPDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double1D postProb;
		private CategoricalDist catDist;
		private HashMap<Integer0D, MVNormalDist> mvnDists; // Shouldn't you be
		// using (and
		// hacking)
		// SeqProbDist for
		// this? FIXME
		private InverseWishartDist sgDist;
		private MatrixNormalDist aDist;

		private Double2D poop;
		private Double2D mPlus;
		private Double2D poopMPlus;
		private Double2D mPlusPoopMplus;
		private double fxConstant;
		
//		private ProbDist<Double2D> baseSg;
//		private ProbDist<Double2D> baseA;
//
//		public ProbDist<Double2D> getBaseSg() {
//			return this.baseSg;
//		}
//
//		public void setBaseSg(ProbDist<Double2D> baseSg) {
//			this.baseSg = baseSg;
//		}
//
//		public ProbDist<Double2D> getBaseA() {
//			return this.baseA;
//		}
//
//		public void setBaseA(ProbDist<Double2D> baseA) {
//			this.baseA = baseA;
//		}
//
		protected Integer0D catVariate() throws ProbDistParmException {
			if (this.catDist == null) {
				this.catDist = new CategoricalDist(this
						.getPostProb());
				return this.catDist.variate();
			} else {
				return this.catDist.variate(this.getPostProb());
			}
		}

		protected double mvnLogDensity(Integer0D k, Double1D pt)
				throws ProbDistParmException {
			return this.mvnDists.get(k).logDensity(pt);
		}

		protected void setUpPostMvn(Integer0D k) throws ProbDistParmException {
			if (this.mvnDists == null) {
				this.mvnDists = new HashMap<Integer0D, MVNormalDist>();
			}
			if (!this.mvnDists.containsKey(k)) {
				double[] zeroD = new double[FLGFWModel.this.dims];
				Double1D zero = new Double1D(zeroD);
				this.mvnDists.put(k, new MVNormalDist(zero, this.getSg().get(
						k.value())));
			}
			// else do nothing

		}

		protected void resetPostMvn(Integer0D k) {
			if (this.mvnDists != null && this.mvnDists.containsKey(k)) {
				this.mvnDists.remove(k);
			}
		}

		private Double1D getPostProb() {
			return this.postProb;
		}

		protected void setPostProb(Double1D prob) {
			this.postProb = prob;
		}

		protected double marginalNew(Double1D pt, Integer1D b) {
			Double2D omegaInv = this.poop.plus(b.outer(b)); // FIXME b is reduced and augmented
			Double2D omegaM = pt.outer(b).transpose().plus(this.poopMPlus);
			Double2D psi = this.getPsi().plus(pt.outer(pt).plus(this.mPlusPoopMplus).minus(omegaM.transpose().mult(omegaInv.inverse()).mult(omegaM)));
			int d = pt.size();
			double lmd = this.fxConstant; 
			lmd -= (d/2) * Math.log(omegaInv.det());
			lmd -= ((this.getKappa().value() + 1)/2)*Math.log(psi.det()); 
			return lmd;
		}

		protected Double2D sgVariate(Double2D psi, double k)
				throws ProbDistParmException {
			if (this.sgDist == null) {
				this.sgDist = new InverseWishartDist(psi, new Double0D(k));
				return this.sgDist.variate();
			} else {
				return this.sgDist.variate(psi, new Double0D(k));
			}
		}

		protected Double2D aVariate(Double2D M, Double2D sg, Double2D omega)
				throws ProbDistParmException {
			if (this.aDist == null) {
				this.aDist = new MatrixNormalDist(M, sg, omega);
				return this.aDist.variate();
			} else {
				return this.aDist.variate(M, sg, omega);
			}
		}

		protected ListSequence<? extends Numeric<?>> getVariateDPPost(Double1D pt,
				Integer1D b) throws ProbDistParmException {
			int d = FLGFWModel.this.dims;
			int[] gamma = this.getGamma().value();
			int h = gamma.length;
			int[] gamma1 = this.getGamma().which(new Integer0D(1)).value();
			int h1 = gamma1.length;
			
			Double2D omegaInv = this.poop.plus(b.outer(b)); // FIXME b is reduced and augmented
			Double2D omegaM = pt.outer(b).transpose().plus(this.poopMPlus);
			Double2D omega = omegaInv.inverse();
			Double2D psi = this.getPsi().plus(pt.outer(pt).plus(this.mPlusPoopMplus).minus(omegaM.transpose().mult(omegaInv.inverse()).mult(omegaM)));

			Double2D sgNew = null;
			Double2D aPlusRed = null;
			boolean goodSample = false;
			while (!goodSample) {
				goodSample = true;
				sgNew = this.sgVariate(psi, this.getKappa().value() + 1);
				aPlusRed = this.aVariate(omega.mult(omegaM), sgNew, omega);
				double[][] a = aPlusRed.value();
				for (int m = 0; m < aPlusRed.numRows(); m++) {
					for (int n = 0; n < aPlusRed.numCols(); n++) {
						if (Double.isInfinite(a[m][n]) || Double.isNaN(a[m][n])) {
							System.err.println("Warning: bad sample for new A+: " + aPlusRed);
							System.err.println("New Sg: " + sgNew);
							System.err.println("Omega: " + omega);
							System.err.println("mean: " + omega.mult(omegaM));
							goodSample = false;
							break;
						}
					}
					if (!goodSample) {
						break;
					}
				}
			}
			
			Double1D zero = new Double1D(new double[d]);
			Double2D aPlusNew = new Double2D(h+1, d);
			for (int i=0; i<h; i++) {
				if (gamma[i] == 0) {
					aPlusNew.set(i, zero);
				}
			}
			for (int i=0; i<h1; i++) {
				aPlusNew.set(gamma1[i], aPlusRed.get(i));
			}
			aPlusNew.set(h1, aPlusRed.get(h1));
			return new ListSequence<Double2D>(sgNew, aPlusNew);
		}

		protected Double0D getRho() {
			return (Double0D) this.chainParms[0];
		}

		protected Double3D getSg() {
			return (Double3D) this.chainParms[1];
		}

		protected Double3D getA() {
			return (Double3D) this.chainParms[2];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[3];
		}

		protected Integer2D getB() {
			return (Integer2D) this.chainParms[4];
		}

		protected Double2D getM() {
			return (Double2D) this.chainParms[5];
		}

		protected Double0D getAl() {
			return (Double0D) this.chainParms[6];
		}
		
		protected Double0D getRho0() {
			return (Double0D) this.chainParms[7];
		}
		
		protected Double2D getA0() {
			return (Double2D) this.chainParms[8];
		}
		
		protected Integer1D getGamma() {
			return (Integer1D) this.chainParms[9];
		}

		protected Double1D getM0() {
			return (Double1D) this.chainParms[10];
		}
		
		protected Double2D getPsi() {
			return (Double2D) this.fixedParms[0];
		}

		protected Double0D getKappa() {
			return (Double0D) this.fixedParms[1];
		}

		// protected Double2D getPhi() {
		public Double2D getPhi() { // FIXME
			return (Double2D) this.fixedParms[2];
		}

		protected Double0D getLambda() {
			return (Double0D) this.fixedParms[3];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(7);
			this.chainParmNames.add("rho");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("A");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.chainParmNames.add("M");
			this.chainParmNames.add("al");
			this.chainParmNames.add("rho0");
			this.chainParmNames.add("A0");
			this.chainParmNames.add("Gamma");
			this.chainParmNames.add("M0");
			this.fixedParmNames = new ArrayList<String>(4);
			this.fixedParmNames.add("Psi");
			this.fixedParmNames.add("kappa");
			this.fixedParmNames.add("beta");
			this.fixedParmNames.add("lambda");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(7);
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					7);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					4);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			this.mvnDists = null;
			this.aDist = null;
			int h1 = this.getGamma().sum().value();
			int[] gamma1 = this.getGamma().which(new Integer0D(1)).value();
			int d = FLGFWModel.this.dims;
			double[][] poopMat = new double[h1+1][h1+1]; 
			if (h1 > 0) {
				for (int i=0; i<h1; i++) {
					poopMat[i][i] = this.getRho().value();
				}
			}
			poopMat[h1][h1] = this.getRho0().value();
			this.poop = new Double2D(poopMat);
			this.mPlus = (Double2D) this.getM().transpose().transpose().getAll(gamma1);
			this.mPlus.add(this.getM0());
			this.poopMPlus = this.poop.mult(this.mPlus);
			this.mPlusPoopMplus = this.mPlus.transpose().mult(this.poopMPlus);
			this.fxConstant = (-d/2)*Math.log(2*Math.PI);
			this.fxConstant -= (d/2)*Math.log(2);
			this.fxConstant += (d/2)*Math.log(this.getRho0().value()) + (h1*d/2)*Math.log(this.getRho().value());
			this.fxConstant += (this.getKappa().value()/2)*Math.log(this.getPsi().det());
			this.fxConstant += Math.log(FLGFWModel.gammaD((this.getKappa().value()+1)/2, d));
			this.fxConstant -= Math.log(FLGFWModel.gammaD(this.getKappa().value()/2, d));
		}

		@Override
		protected Integer1D genVariate() throws ProbDistParmException {
			// FIXME - This sampling scheme destructively modifies
			// A, Sg - HACK
			int N = this.getSamplerData().size();
			Integer1D Z = null;
			try {
				Z = (Integer1D) this.getZ().clone();
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace(); // FIXME
			}
			Integer0D unused = new Integer0D(-1);
			int h = this.getGamma().size();
			int[] effectRows = new int[h];
			for (int i=0; i<h; i++) {
				effectRows[i] = i;
			}	
			
			int[] gamma1 = this.getGamma().which(new Integer0D(1)).value();
			Integer2D bPlusT = this.getB().transpose().getAll(gamma1);
			int[] ones = new int[this.getSamplerData().size()];
			for (int j=0; j<ones.length; j++) {
				ones[j] = 1;
			}
			bPlusT.add(new Integer1D(ones)); 	
			Integer2D bPlus = bPlusT.transpose();
			
			for (int i = 0; i < N; i++) {
				Double1D Xi = this.getSamplerData().get(i);
				Integer1D active = Z.items();
				Z.set(i, unused);
				double[] logP = new double[active.size() + 1];
				double maxLogP = -Double.MAX_VALUE;
				int ai = 0;
				// Existing categories
				Integer1D b = this.getB().get(i);
				for (Integer0D k : active) {
					Integer1D zK = Z.which(k);
					int nZK = zK.size();
					if (nZK > 0) {
						Double1D aTB = b.mult(this.getA().get(k.value()));
						Double1D xIA = this.getSamplerData().get(i).minus(aTB).minus(this.getA0().get(k.value()));
						this.setUpPostMvn(k);
						logP[ai] = this.mvnLogDensity(k, xIA) + Math.log(nZK);
						if (logP[ai] == Double.POSITIVE_INFINITY) {
							logP[ai] = Double.MAX_VALUE;
						}
						if (logP[ai] > maxLogP) {
							maxLogP = logP[ai];
						}
					} else {
						logP[ai] = -Double.MAX_VALUE;
					}
					ai++;
				}
				// New category
				int newc = ai;
				Integer1D bp = bPlus.get(i);
				logP[newc] = this.marginalNew(Xi, bp)
						+ Math.log(this.getAl().value());
				if (logP[ai] == Double.POSITIVE_INFINITY) {
					logP[ai] = Double.MAX_VALUE;
				}
				if (logP[newc] > maxLogP) {
					maxLogP = logP[newc];
				}
				// Select a category for this point
				Double1D p = (new Double1D(logP)).minus(maxLogP).exp();
				this.setPostProb(p);
				int az = this.catVariate().value();
				// Sample new category if necessary
				Integer0D zedNew = null;
				if (az == newc) {
					int maxActive = active.get(active.size() - 1).value();
					zedNew = Z.minNotIn(0, maxActive);
					this.resetPostMvn(zedNew);
					ListSequence<? extends Numeric<?>> newcParms = this.getVariateDPPost(
							Xi, bp);
					Double2D newSg = (Double2D) newcParms.get(0);
					Double2D newAPlus = (Double2D) newcParms.get(1);

					Double2D newA = newAPlus.getAll(effectRows);
					Double1D newA0 = newAPlus.get(h);
					this.getA().set(zedNew.value(), newA);
					this.getA0().set(zedNew.value(), newA0);
					this.getSg().set(zedNew.value(), newSg);
				} else {
					zedNew = active.get(az);
				}
				Z.set(i, zedNew);
			}
			return Z;
		}

		@Override
		protected double getDensity(Integer1D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PostXDist extends ProbDistInitializeByChain<Double0D> {
		private BetaDist betaDist;
		private Double0D postXA;
		private Double0D postXB;

		public PostXDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(1);
			this.chainParmNames.add("al");
			this.fixedParmNames = new ArrayList<String>(0);
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					1);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					0);
		}

		private Double0D betaVariate() throws ProbDistParmException {
			if (this.betaDist == null) {
				this.betaDist = new BetaDist(this.getPostXA(), this.getPostXB());
				return this.betaDist.variate();
			} else {
				return this.betaDist
						.variate(this.getPostXA(), this.getPostXB());
			}
		}

		protected Double0D getAl() {
			return (Double0D) this.chainParms[0];
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPostXA(this.getAl().plus(1));
			this.setPostXB(new Double0D((double) this.getSamplerData().size()));
		}

		private Double0D getPostXA() {
			return this.postXA;
		}

		private void setPostXA(Double0D xA) {
			this.postXA = xA;
		}

		private Double0D getPostXB() {
			return this.postXB;
		}

		private void setPostXB(Double0D xB) {
			this.postXB = xB;
		}

		@Override
		protected Double0D genVariate() throws ProbDistParmException {
			return this.betaVariate();
		}

		@Override
		protected double getDensity(Double0D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PostAlDist extends ProbDistInitializeByChain<Double0D> {
		private GammaDist gammaDist;
		private CategoricalDist catDist;
		private Double0D postAlA;
		private Double0D postAlB;
		private Double1D postMix;

		public PostAlDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(2);
			this.chainParmNames.add("x");
			this.chainParmNames.add("Z");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("ala");
			this.fixedParmNames.add("alb");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		private Double0D gammaVariate() throws ProbDistParmException {
			if (this.gammaDist == null) {
				this.gammaDist = new GammaDist(this.getPostAlA(),
						this.getPostAlB());
				return this.gammaDist.variate();
			} else {
				return this.gammaDist.variate(this.getPostAlA(),
						this.getPostAlB());
			}
		}

		private Integer0D catVariate() throws ProbDistParmException {
			if (this.catDist == null) {
				this.catDist = new CategoricalDist(
						this.getPostMix());
				return this.catDist.variate();
			} else {
				return this.catDist.variate(this.getPostMix());
			}
		}

		protected Double0D getX() {
			return (Double0D) this.chainParms[0];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[1];
		}

		protected Double0D getAlA() {
			return (Double0D) this.fixedParms[0];
		}

		protected Double0D getAlB() {
			return (Double0D) this.fixedParms[1];
		}

		@Override
		protected void setUpFromChainParms() {
			int K = this.getZ().items().size();
			int N = this.getSamplerData().size();
			double logX = Math.log(this.getX().value());
			double p1 = this.getAlA().value() + K - 1;
			double p2 = N * (this.getAlB().value() - logX);
			this.setPostMix(new Double1D(p1, p2));
			int i;
			try {
				i = this.catVariate().value();
			} catch (ProbDistParmException e) {
				// FIXME
				throw new RuntimeException(e);
			}
			if (i == 0) {
				this.setPostAlA(this.getAlA().plus(K));
			} else {
				this.setPostAlA(this.getAlA().plus(K - 1));
			}
			this.setPostAlB(this.getAlB().plus(-logX)); // FIXME
		}

		private Double0D getPostAlA() {
			return this.postAlA;
		}

		private void setPostAlA(Double0D alA) {
			this.postAlA = alA;
		}

		private Double0D getPostAlB() {
			return this.postAlB;
		}

		private void setPostAlB(Double0D alB) {
			this.postAlB = alB;
		}

		private Double1D getPostMix() {
			return this.postMix;
		}

		private void setPostMix(Double1D mix) {
			this.postMix = mix;
		}

		@Override
		protected Double0D genVariate() throws ProbDistParmException {
			// FIXME - create
			return this.gammaVariate();
		}

		@Override
		protected double getDensity(Double0D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PostBDist extends ProbDistInitializeByChain<Integer2D> {
		public PostBDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private Integer2D getB() {
			return (Integer2D) this.chainParms[0];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(1);
			this.chainParmNames.add("B");
			this.fixedParmNames = new ArrayList<String>(0);
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
			this.chainParmCheck.add(null); // FIXME
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					1);
			this.chainParmClasses.add(Integer2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					0);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		@Override
		protected Integer2D genVariate() throws ProbDistParmException {
			return this.getB();
		}

		@Override
		protected double getDensity(Integer2D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}
	

	class PostA0Dist extends ProbDistInitializeByChain<Double2D> { // FIXME cause its off somewhere else
		public PostA0Dist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private Double2D getA0() {
			return (Double2D) this.chainParms[0];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(1);
			this.chainParmNames.add("A0");
			this.fixedParmNames = new ArrayList<String>(0);
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
			this.chainParmCheck.add(null); // FIXME
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					1);
			this.chainParmClasses.add(Double2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					0);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		@Override
		protected Double2D genVariate() throws ProbDistParmException {
			return this.getA0();
		}

		@Override
		protected double getDensity(Double2D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}
		
	
	class PostMDist extends ProbDistInitializeByChain<Double2D> {
		public PostMDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double3D getA() {
			return (Double3D) this.chainParms[0];
		}

		protected Double3D getSg() {
			return (Double3D) this.chainParms[1];
		}

		protected Double2D getM() {
			return (Double2D) this.chainParms[2];
		}
		
		protected Integer1D getGamma() {
			return (Integer1D) this.chainParms[3];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[4];
		}

		protected Double0D getRho() {
			return (Double0D) this.chainParms[5];
		}
		
		protected Double2D getS() {
			return (Double2D) this.fixedParms[0];
		}
		
		protected Double0D getBeta() {
			return (Double0D) this.fixedParms[1];
		}

		protected Double0D getLambda() {
			return (Double0D) this.fixedParms[2];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(5);
			this.chainParmNames.add("A");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("M");
			this.chainParmNames.add("Gamma");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("rho");
			this.fixedParmNames = new ArrayList<String>(0);
			this.fixedParmNames.add("S");
			this.fixedParmNames.add("beta");
			this.fixedParmNames.add("lambda");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					4);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					0);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		@Override
		protected Double2D genVariate() throws ProbDistParmException {
			Integer1D active = this.getZ().items();
			int[] gamma = this.getGamma().value();
			int h = gamma.length;
			Double1D zero = new Double1D(new double[FLGFWModel.this.dims]);
			Double2D mNew = null;
			try {
				mNew = (Double2D) this.getM().clone();
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			for (int i=0; i<h; i++) {
				if (gamma[i] == 0) {
					mNew.set(i, zero);
				}
			}			
			MVNormalDist mvnDist = new MVNormalDist();
			Double2D sInv = this.getS().inverse();
			Double2D sgA = null;
			for (Integer0D k : active) {
				Double2D sgRho = this.getSg().get(k.value()).inverse().mult(this.getRho());
				sInv = sInv.plus(sgRho);
				if (sgA == null) {
					sgA = sgRho.mult(this.getA().get(k.value()).transpose());
				} else {
					sgA = sgA.plus(sgRho.mult(this.getA().get(k.value()).transpose()));
				}
			}	
			sgA = sgA.transpose();
			for (int i=0; i<h; i++) {
				if (gamma[i] == 1) {
					mNew.set(i, mvnDist.variate(sgA.get(i), sInv.inverse()));
				}
			}			
			return mNew;
		}

		@Override
		protected double getDensity(Double2D pt) throws ProbDistParmException {
			// TODO Auto-generated method stub
			return 0;
		}
	}
	
	class PostM0Dist extends ProbDistInitializeByChain<Double1D> {
		public PostM0Dist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double2D getA0() {
			return (Double2D) this.chainParms[0];
		}

		protected Double3D getSg() {
			return (Double3D) this.chainParms[1];
		}

		protected Double2D getM() {
			return (Double2D) this.chainParms[2];
		}
		
		protected Integer1D getGamma() {
			return (Integer1D) this.chainParms[3];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[4];
		}

		protected Double1D getM0() {
			return (Double1D) this.chainParms[5];
		}
		
		protected Double0D getRho0() {
			return (Double0D) this.chainParms[6];
		}

		protected Double2D getS0() {
			return (Double2D) this.fixedParms[0];
		}
		
		protected Double0D getBeta() {
			return (Double0D) this.fixedParms[1];
		}

		protected Double0D getLambda() {
			return (Double0D) this.fixedParms[2];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(5);
			this.chainParmNames.add("A0");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("M");
			this.chainParmNames.add("Gamma");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("M0");
			this.chainParmNames.add("rho0");
			this.fixedParmNames = new ArrayList<String>(0);
			this.fixedParmNames.add("S0");
			this.fixedParmNames.add("beta");
			this.fixedParmNames.add("lambda");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					4);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					0);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		@Override
		protected Double1D genVariate() throws ProbDistParmException {
			Integer1D active = this.getZ().items();
			MVNormalDist mvnDist = new MVNormalDist();
			Double2D sInv = this.getS0().inverse();
			Double1D sgA0 = null;
			for (Integer0D k : active) {
				Double2D sgRho = this.getSg().get(k.value()).inverse().mult(this.getRho0());
				sInv = sInv.plus(sgRho);
				if (sgA0 == null) {
					sgA0 = sgRho.mult(this.getA0().get(k.value()));
				} else {
					sgA0 = sgA0.plus(sgRho.mult(this.getA0().get(k.value())));
				}
			}	
			return mvnDist.variate(sgA0, sInv.inverse());
		}	
		

		@Override
		protected double getDensity(Double1D pt) throws ProbDistParmException {
			// TODO Auto-generated method stub
			return 0;
		}
	}
	
	class PostRhoDist extends ProbDistInitializeByChain<Double0D> {
		public PostRhoDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double3D getA() {
			return (Double3D) this.chainParms[0];
		}

		protected Double3D getSg() {
			return (Double3D) this.chainParms[1];
		}

		protected Double2D getM() {
			return (Double2D) this.chainParms[2];
		}
		
		protected Integer1D getGamma() {
			return (Integer1D) this.chainParms[3];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[4];
		}

		protected Double0D getBeta() {
			return (Double0D) this.fixedParms[1];
		}

		protected Double0D getLambda() {
			return (Double0D) this.fixedParms[0];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(5);
			this.chainParmNames.add("A");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("M");
			this.chainParmNames.add("Gamma");
			this.chainParmNames.add("Z");
			this.fixedParmNames = new ArrayList<String>(0);
			this.fixedParmNames.add("lambda");
			this.fixedParmNames.add("beta");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					4);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					0);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		@Override
		protected Double0D genVariate() throws ProbDistParmException {
			int d = this.getM().numCols();
			Integer1D active = this.getZ().items();
			int K = active.size();
			int[] gamma1 = this.getGamma().which(new Integer0D(1)).value();
			int h1 = gamma1.length;
			Double0D l = this.getLambda().plus(d*h1*K/2);
			Double0D b = this.getBeta();
			for (Integer0D k : active) {
				for (int i=0; i<h1; i++) {
					Double1D ac = this.getA().get(k.value()).get(gamma1[i]).minus(this.getM().get(gamma1[i]));
					b = b.plus(0.5*ac.mult(this.getSg().get(k.value()).inverse()).mult(ac));
				}
			}
			GammaDist gammaDist = new GammaDist(l, b);
			Double0D v = gammaDist.variateFast();
			return v;
		}

		@Override
		protected double getDensity(Double0D pt) throws ProbDistParmException {
			// TODO Auto-generated method stub
			return 0;
		}
	}
	
	class PostRho0Dist extends ProbDistInitializeByChain<Double0D> {
		public PostRho0Dist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double2D getA0() {
			return (Double2D) this.chainParms[0];
		}

		protected Double3D getSg() {
			return (Double3D) this.chainParms[1];
		}

		protected Double1D getM0() {
			return (Double1D) this.chainParms[2];
		}
		
		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[3];
		}

		protected Double0D getBeta0() {
			return (Double0D) this.fixedParms[1];
		}

		protected Double0D getLambda0() {
			return (Double0D) this.fixedParms[0];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(5);
			this.chainParmNames.add("A0");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("M0");
			this.chainParmNames.add("Z");
			this.fixedParmNames = new ArrayList<String>(0);
			this.fixedParmNames.add("lambda0");
			this.fixedParmNames.add("beta0");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					4);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					0);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		@Override
		protected Double0D genVariate() throws ProbDistParmException {
			int d = this.getM0().size();
			Integer1D active = this.getZ().items();
			int K = active.size();
			Double0D l = this.getLambda0().plus(d*K/2);
			Double0D b = this.getBeta0();
			for (Integer0D k : active) {
				Double1D ac = this.getA0().get(k.value()).minus(this.getM0());
				b = b.plus(0.5*ac.mult(this.getSg().get(k.value()).inverse()).mult(ac));
			}
			GammaDist gammaDist = new GammaDist(l, b);
			return gammaDist.variateFast();
		}

		@Override
		protected double getDensity(Double0D pt) throws ProbDistParmException {
			// TODO Auto-generated method stub
			return 0;
		}
	}
		
	@SuppressWarnings("unchecked")
	public static ChainLink pointEstimate(List<ChainLink> c, Double2D d)
			throws ProbDistParmException {
		// MAP
		double max = Double.NEGATIVE_INFINITY;
		int argmax = c.size() - 1;
		MVNormalDist mvnDist = new MVNormalDist();
		for (int i = 0; i < c.size(); i++) {
			double apost = 0;
			ChainLink cli = c.get(i);
			// Hyperparameters
			RandomVar<Double0D> rvAl = (RandomVar<Double0D>) cli.get("al");
			Double0D al = rvAl.getNumericValue();
			apost += rvAl.getPrior().logDensity(al);
			// Model parameters and likelihood
			RandomVar<Double3D> rvSg = (RandomVar<Double3D>) cli.get("Sg");
			RandomVar<Double3D> rvA = (RandomVar<Double3D>) cli.get("A");
			// Gamma
			RandomVar<Integer1D> rvGamma = (RandomVar<Integer1D>) cli.get("Gamma");
			Integer1D gamma = rvGamma.getNumericValue();
			// B
			RandomVar<Integer2D> rvB = (RandomVar<Integer2D>) cli.get("B");
			Integer2D B = rvB.getNumericValue();
			// skip
			// Z
			RandomVar<Integer1D> rvZ = (RandomVar<Integer1D>) cli.get("Z");
			Integer1D Z = rvZ.getNumericValue();
			for (int j = 0; j < Z.size(); j++) {
				int[] zIndices = new int[j];
				for (int k = 0; k < j; k++) {
					zIndices[k] = k;
				}
				Integer1D ZPreJ = Z.getAll(zIndices);
				int nZ = ZPreJ.which(new Integer0D(Z.getValue(j))).size();
				if (nZ != 0) {
					apost += Math.log(nZ) - Math.log(ZPreJ.size() + al.value());
				} else {
					apost += Math.log(al.value())
							- Math.log(ZPreJ.size() + al.value());
				}
			}
			Integer1D active = Z.items();
			for (Integer0D z : active) {
				try {
					// Model parameters
					Double2D Sg = ((Double3D) rvSg.getNumericValue()).get(z
							.value());
					Double2D A = ((Double3D) rvA.getNumericValue()).get(z
							.value());
//					CRPDist postZ = (CRPDist) rvZ.getPosterior();
//					apost += postZ.getBaseSg().logDensity(Sg);
//					apost += postZ.getBaseA().logDensity(A, M, Sg, // EEK FIXME
//							postZ.getPhi());
					// Likelihood
					Integer1D whichz = Z.which(z);
					Integer2D db = B.getAll(whichz.value());
					Double2D dz = d.getAll(whichz.value());
					for (int j = 0; j < dz.size(); j++) {
						Double1D pt = dz.get(j);
						Integer1D b = db.get(j);
						Double1D mean = b.mult(A);
						Double2D cov = Sg;
						apost += mvnDist.logDensity(pt, mean, cov);
					}
				} catch (Exception e) {
					System.err.println(e.getMessage());
					StackTraceElement[] st = e.getStackTrace();
					for (int j = 0; j < st.length; j++) {
						System.err.println(st[j].toString());
					}
					throw new ProbDistParmException(e);
				}
			}
			// Set max
			if (apost > max) {
				argmax = i;
				max = apost;
			}
		}
		return c.get(argmax);
	}
	
	public FLGFWModel() {
		super();
	}

	public FLGFWModel(Map<String, Numeric<? extends Numeric<?>>> hypers,
			Map<String, Numeric<? extends Numeric<?>>> init, int dims)
			throws ProbDistParmException {
		super(hypers, dims);
		// Set up parameters
		this.params = new HashMap<String, RandomVar<? extends Numeric<?>>>();

		// Sg
		PriorSgDist priorSg = new PriorSgDist((Double2D) this.getHyper("Psi"),
				(Double0D) this.getHyper("kappa"));
		PostSgDist postSg = new PostSgDist((Double2D) this.getHyper("Psi"),
				(Double0D) this.getHyper("kappa"));
		Double3D defaultSg = (Double3D) init.get("Sg"); // FIXME
		RandomVar<Double3D> rvSg = new RandomVar<Double3D>("Sg", priorSg,
				postSg, defaultSg);
		this.params.put("Sg", rvSg);
		
		// M0
		MVNormalDist priorM0 = new MVNormalDist();
		PostM0Dist postM0 = new PostM0Dist((Double2D) this.getHyper("S"), (Double0D) this.getHyper("beta"), (Double0D) this.getHyper("lambda"));
		Double1D defaultM0 = (Double1D) init.get("M0"); // FIXME
		RandomVar<Double1D> rvM0 = new RandomVar<Double1D>("M0", priorM0, postM0,
				defaultM0);
		this.params.put("M0", rvM0);
		
		// M
		MatrixNormalDist priorM = new MatrixNormalDist();
		PostMDist postM = new PostMDist((Double2D) this.getHyper("S0"), (Double0D) this.getHyper("beta"), (Double0D) this.getHyper("lambda"));
		Double2D defaultM = (Double2D) init.get("M"); // FIXME
		RandomVar<Double2D> rvM = new RandomVar<Double2D>("M", priorM, postM,
				defaultM);
		this.params.put("M", rvM);

		// A0
		PriorA0Dist priorA0 = new PriorA0Dist();
		PostA0Dist postA0 = new PostA0Dist();
		Double2D defaultA0 = (Double2D) init.get("A0"); // FIXME
		RandomVar<Double2D> rvA0 = new RandomVar<Double2D>("A0", priorA0, postA0, defaultA0);
		this.params.put("A0", rvA0);
		
		// A
		PriorADist priorA = new PriorADist();
		PostADist postA = new PostADist((Double0D) this.getHyper("lambda"), (Double0D) this.getHyper("beta"));
		Double3D defaultA = (Double3D) init.get("A"); // FIXME
		RandomVar<Double3D> rvA = new RandomVar<Double3D>("A", priorA, postA, defaultA);
		this.params.put("A", rvA);

		// rho
		GammaDist priorRho = new GammaDist((Double0D) this.getHyper("lambda"), (Double0D) this.getHyper("beta"));
		PostRhoDist postRho = new PostRhoDist((Double0D) this.getHyper("lambda"), (Double0D) this.getHyper("beta"));
		Double0D defaultRho = (Double0D) init.get("rho"); // FIXME
		RandomVar<Double0D> rvRho = new RandomVar<Double0D>("rho", priorRho, postRho, defaultRho);
		this.params.put("rho", rvRho);
		
		// rho0
		GammaDist priorRho0 = new GammaDist((Double0D) this.getHyper("lambda0"), (Double0D) this.getHyper("beta0"));
		PostRho0Dist postRho0 = new PostRho0Dist((Double0D) this.getHyper("lambda0"), (Double0D) this.getHyper("beta0"));
		Double0D defaultRho0 = (Double0D) init.get("rho0"); // FIXME
		RandomVar<Double0D> rvRho0 = new RandomVar<Double0D>("rho0", priorRho0, postRho0, defaultRho0);
		this.params.put("rho0", rvRho0);		
		
		// Z
		ProbDist<Integer1D> priorZ = new Deterministic<Integer1D>(init.get("Z"));
		CRPDist postZ = new CRPDist((Double2D) this.getHyper("Psi"),
				(Double0D) this.getHyper("kappa"),
				(Double0D) this.getHyper("beta"), (Double0D) this.getHyper("lambda")); // FIXME
		// --
		// Hack!
//		postZ.setBaseSg(new InverseWishartDist(this.getHyper("Psi"), this.getHyper("kappa"))); // FIXME -- Eek! Hack!
//		postZ.setBaseA(new MatrixNormalDist()); // FIXME -- Eek! Hack!
//		postZ.setBaseA0(new MVNormalDist()); // FIXME -- Eek! Hack!
		Integer1D defaultZ = (Integer1D) init.get("Z");
		RandomVar<Integer1D> rvZ = new RandomVar<Integer1D>("Z", priorZ, postZ, defaultZ);
		this.params.put("Z", rvZ); // FIXME
		
		// Gamma
		ProbDist<Integer1D> priorGamma = new Deterministic<Integer1D>(init.get("Gamma"));
		PostGammaDist postGamma = new PostGammaDist(hypers.get("p"), hypers.get("lambda"), hypers.get("beta"), hypers.get("S"), hypers.get("lambda0"), hypers.get("beta0"), hypers.get("Psi"), hypers.get("kappa"), hypers.get("S0"));
		Integer1D defaultGamma = (Integer1D) init.get("Gamma"); // FIXME
		RandomVar<Integer1D> rvGamma = new RandomVar<Integer1D>("Gamma", priorGamma, postGamma, defaultGamma);
		this.params.put("Gamma", rvGamma); // FIXME
		

		// B
		ProbDist<Integer2D> priorB = new Deterministic<Integer2D>(init.get("B"));
		PostBDist postB = new PostBDist();
		Integer2D defaultB = (Integer2D) init.get("B"); // FIXME
		RandomVar<Integer2D> rvB = new RandomVar<Integer2D>("B", priorB, postB,
				defaultB);
		this.params.put("B", rvB); // FIXME

		// x
		ProbDist<Double0D> priorX = new BetaDist(
				(Double0D) this.getHyper("xa"), (Double0D) this.getHyper("xb"));
		PostXDist postX = new PostXDist();
		Double0D defaultX = (Double0D) init.get("x"); // FIXME
		RandomVar<Double0D> rvX = new RandomVar<Double0D>("x", priorX, postX,
				defaultX);
		this.params.put("x", rvX); // FIXME

		// al
		ProbDist<Double0D> priorAl = new GammaDist(
				(Double0D) this.getHyper("ala"),
				(Double0D) this.getHyper("alb"));
		ProbDistInitializeByChain<Double0D> postAl = new PostAlDist(
				(Double0D) this.getHyper("ala"),
				(Double0D) this.getHyper("alb"));
		Double0D defaultAl = (Double0D) init.get("al"); // FIXME
		RandomVar<Double0D> rvAl = new RandomVar<Double0D>("al", priorAl,
				postAl, defaultAl);
		this.params.put("al", rvAl); // FIXME
	}

	@Override
	public ChainLink getInitialLink() {
		return new ChainLink(this.getParam("Z"), this.getParam("B"), 
				this.getParam("Gamma"), this.getParam("A"), this.getParam("rho0"), this.getParam("rho"), this.getParam("Sg"),  this.getParam("x"), this.getParam("A0"),
				this.getParam("al"), this.getParam("M"), this.getParam("M0"));
	}
}