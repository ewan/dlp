package org.jqgibbs.models;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
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
import org.jqgibbs.mathstat.probdist.ProbDistMC;
import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;
import org.jqgibbs.mathstat.probdist.SeqInverseWishartDist;
import org.jqgibbs.mathstat.probdist.SeqMVNormalDist;
import org.jqgibbs.mathstat.probdist.SeqMatrixNormalDist;
import cern.colt.matrix.linalg.SingularValueDecomposition;

import cern.jet.random.Uniform;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.stat.Gamma;

public class FLGFVModel extends Model {

	private static Uniform unifGen = new Uniform(
			new MersenneTwister(new Date()));
	
	private static Double2D rmvt(Double2D m, Double2D sg, double df) {
		int h = m.size();
		int d = m.numCols();
		double[] rg = new double[h];
		cern.jet.random.Gamma gammaDist = new cern.jet.random.Gamma(df/2, 2, new MersenneTwister(new Date()));
		for (int i=0; i<h; i++) {
			rg[i] = gammaDist.nextDouble();
		}
		MVNormalDist mvnDist = null;
		double[][] rtkp = new double[h][d];
		for (int i=0; i<h; i++) {
			try {
				mvnDist = new MVNormalDist(m.get(i), sg);
				rtkp[i] = mvnDist.variateFast().value();
			} catch (ProbDistParmException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return new Double2D(rtkp);
	}
	
	private static double aTFnNoexp(Double1D[] aK, Double1D m, Double3D sgK, double beta, int d, Integer1D active, double lambda) {
		double s = 0;
		for (Integer0D k : active) {
			Double1D am = aK[k.value()].minus(m);
			s += am.mult(sgK.get(k.value()).mult(1/beta)).mult(am);
		}
		return(1+s/2);	
	}
	
	private static double aTFn(Double1D[] aK, Double1D m, Double3D sgK, double beta, int d, Integer1D active, double lambda) {
		int K = active.size();	
		return Math.pow(aTFnNoexp(aK,m,sgK,beta,d,active,lambda), -d*K-lambda);
	}
	
	private static double logADistMCInt(Double1D[] beK, double[] bsum, int d, Double3D sgK, Integer1D active, Double1D m, double lambda, double beta) 
			throws ProbDistParmException {
		int MCNSAMP = 10;
		// step 1: generate samples from mvn
		int K = sgK.size();
		MVNormalDist mvnDist;
		Double1D[][] samples = new Double1D[MCNSAMP][K];
		for (Integer0D k : active) {
			mvnDist = new MVNormalDist(beK[k.value()].mult(1/bsum[k.value()]), sgK.get(k.value()).mult(1/Math.pow(bsum[k.value()], 2)));
			for (int i=0; i<MCNSAMP; i++) {
				samples[i][k.value()] = mvnDist.variateFast();
			}
		}
		// step 2: compute average of f(n) over all samples
		double logBsum = 0;
		for (Integer0D k : active) {
			logBsum += Math.log(bsum[k.value()]);
		}
		logBsum *= d/2;
		double v = 0;
		for (int i=0; i<MCNSAMP; i++) {
			v += aTFn(samples[i], m, sgK, beta, d, active, lambda);
		}
		return Math.log(v/MCNSAMP) + logBsum;
	}
	
	class PostGammaDist extends ProbDistMC<Integer1D> {
		public PostGammaDist(Numeric<?>... fixed)
				throws ProbDistParmException {
			super(fixed);
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

		protected Double0D getP() {
			return (Double0D) this.fixedParms[0];
		}

		protected Double0D getLambda() {
			return (Double0D) this.fixedParms[1];
		}

		protected Double0D getBeta() {
			return (Double0D) this.fixedParms[2];
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
			this.fixedParmNames = new ArrayList<String>(3);
			this.fixedParmNames.add("p");
			this.fixedParmNames.add("lambda");
			this.fixedParmNames.add("beta");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(6);
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
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					3);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			return; // FIXME??
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
			int d = FLGFVModel.this.dims;
			Double1D[] be = new Double1D[K];
			double[] bsum = new double[K];
			for (int i = 0; i < h; i++) {
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
					double[][] A = this.getA().get(k.value()).value();
					double[][] redA = new double[rowsRed][d];
					int[] B = bK[k.value()].sum().value();
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
				// bsum[k] := \sum_j (b_{i,j})
				for (Integer0D k : activeZ) {
					double bs = 0;
					int[][] bKval = bK[k.value()].value();
					for (int j=0; j<bKval.length; j++) {
						bs += bKval[j][i];
					}
					bsum[k.value()] = bs;
				}
				// draw gamma
				// p(g=1|X)
				double logP = Math.log(this.getP().value());
				Double2D sgK;
				for (Integer0D k : activeZ) {
					sgK = this.getSg().get(k.value());
					logP += Math.log(Gamma.gamma(d*Ka/2 + this.getLambda().value()));
					logP -= Math.log(Gamma.gamma(this.getLambda().value()));
					logP -= (this.getLambda().value()+1)*Math.log(this.getBeta().value());
					int nK = xK[k.value()].numRows();
					logP -= nK*Math.log(this.getSg().get(k.value()).det());
					logP -= (d/2)*nK*Math.log(2*Math.PI);
				}
				logP += FLGFVModel.logADistMCInt(be, bsum, d, this.getSg(), activeZ, this.getM().get(i), this.getLambda().value(), this.getBeta().value());
				// p(g=0|X)
				// 1. 1-p
				double logQ = Math.log(1-this.getP().value());
				// 2. exp{\sum_k be[k]^T \Sigma_k^{-1} be[k]}
				for (Integer0D k : activeZ) {
					sgK = this.getSg().get(k.value());	
					int nK = xK[k.value()].size();
					if (nK > 0) {
						logQ -= 0.5*be[k.value()].mult(sgK.inverse()).mult(be[k.value()]);
					}
				}
				// Normalizing constant
				double logZ = Math.log(Math.exp(logP)+Math.exp(logQ));
				double u = Math.log(FLGFVModel.unifGen.nextDouble());
				if (u < logP-logZ) {
					gamma[i] = 1;
				} else {
					gamma[i] = 0;
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
	
	class PriorSgDist extends ProbDistMC<Double3D> {
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

	class PostSgDist extends ProbDistMC<Double3D> {
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

		private Double3D iwVariates() throws ProbDistParmException {
			if (this.iwDists == null) {
				this.iwDists = new SeqInverseWishartDist(this.getPostPsi(),
						this.getPostKappa());
				return this.iwDists.variate();
			} else {
				return this.iwDists.variate(this.getPostPsi(),
						this.getPostKappa());
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
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("Psi");
			this.fixedParmNames.add("kappa");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
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
			for (Integer0D k : active) {
				this.getPostPsi().add(this.getPsi());
				Integer1D zK = this.getZ().which(k);
				Double2D xK = this.getSamplerData().getAll(zK.value());
				this.getPostKappa().add(this.getKappa().plus(xK.size()));
			}
		}

		@Override
		protected Double3D genVariate() throws ProbDistParmException {
			int d = FLGFVModel.this.dims;
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

	class PriorADist extends ProbDistMC<Double3D> {
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
					FLGFVModel.this.dims);
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

	class PostADist extends ProbDistMC<Double3D> {
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
			this.chainParmNames.add("gamma");
			this.chainParmNames.add("A0");
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
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		// It would be far more efficient to sample from MVT and
		// accept/reject based on MVN!
		@Override
		protected Double3D genVariate() throws ProbDistParmException {
			int[] gamma =  this.getGamma().value();
			int h = gamma.length;
			Integer1D active = this.getZ().items();
			int K = active.size();
			int d = FLGFVModel.this.dims;
			Integer1D[] zK = new Integer1D[this.getSg().size()];
			Double2D[] xK = new Double2D[this.getSg().size()];
			Integer2D[] bK = new Integer2D[this.getSg().size()];
			for (Integer0D k : active) {
				zK[k.value()] = this.getZ().which(k);
				xK[k.value()] = this.getSamplerData().getAll(zK[k.value()].value());
				bK[k.value()] = this.getB().getAll(zK[k.value()].value());
			}
			Double3D aNew = null;
			try {
				aNew = (Double3D) this.getA().clone();
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			Double2D[] aNewA = (Double2D[]) aNew.toArray();
			double[] z = new double[d];
			for (int i=0; i<d; i++) {
				z[0] = 0;
			}
			Double1D zero = new Double1D(z);
			for (int i=0; i<h; i++) {
				if (gamma[i] == 0) {
					for (Integer0D k : active) {
						aNewA[k.value()].set(i, zero);
					}
				} 
			}					
			double exp = -(d*K + this.getLambda().value());
			// be[k] := \sum_j (x_j - a_{0,k} - \sum_{h-i} b_{h,j} * a_{h,k})
			Double1D[] be = new Double1D[this.getSg().size()];
			// bsum[k] := \sum_j (b_{i,j})
			double[] bsum = new double[this.getSg().size()];
			Integer1D B;
			Double1D ba;
			for (int i=0; i<h; i++) {
				if (gamma[i] == 1) {
					for (Integer0D k : active) {
						aNewA[k.value()].set(i, zero);
						// \sum_j x_j
						be[k.value()] = xK[k.value()].sum();
						// \sum_j a_{0,k}
						int nK = xK[k.value()].numRows();
						be[k.value()] = be[k.value()].minus(this.getA0().get(k.value()).mult(nK));	
						// \sum_j \sum_{h-i} b_{h,j} * a_{h,k}
						B = bK[k.value()].sum();
						ba = B.mult(aNewA[k.value()]);
						be[k.value()] = be[k.value()].minus(ba);
						// (\sum_j (b_{i,j}))^2
						double bs = 0;
						int[][] bKval = bK[k.value()].value();
						for (int j=0; j<bKval.length; j++) {
							bs += bKval[j][i];
						}
						bsum[k.value()] = bs;
					}
					double logP;
					double logNum;
					double logDen;
					Double1D mi = this.getM().get(i);
					Double1D[] aKCurr = new Double1D[this.getSg().size()];
					Double1D[] aKProp = new Double1D[this.getSg().size()];
					MVNormalDist mvnDist;
					// sample from relevant gaussians
					for (Integer0D k : active) {
						aKCurr[k.value()] = this.getA().get(k.value()).get(i);
						mvnDist = new MVNormalDist(be[k.value()].mult(1/bsum[k.value()]), this.getSg().get(k.value()).mult(1/Math.pow(bsum[k.value()], 2)));
						aKProp[k.value()] = mvnDist.variateFast();
					}
					logP = Math.log(FLGFVModel.unifGen.nextDouble());
					logNum = Math.log(aTFnNoexp(aKProp, mi, this.getSg(), this.getBeta().value(), d, active, this.getLambda().value()));
					logDen = Math.log(aTFnNoexp(aKCurr, mi, this.getSg(), this.getBeta().value(), d, active, this.getLambda().value()));
					if (logP < exp*(logNum-logDen)) {
						for (Integer0D k : active) {
							aNewA[k.value()].set(i, aKProp[k.value()]);
						}
					} else {
						for (Integer0D k : active) {
							aNewA[k.value()].set(i, aKCurr[k.value()]);
						}
					}
				}
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
	
	
	class PriorA0Dist extends ProbDistMC<Double2D> {
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

	class PostA0Dist extends ProbDistMC<Double2D> {
		public PostA0Dist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
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

		protected Integer1D getGamma() {
			return (Integer1D) this.chainParms[4];
		}
		
		protected Double2D getA0() {
			return (Double2D) this.chainParms[5];
		}

		protected Double1D getM0() {
			return (Double1D) this.chainParms[6];
		}
		
		protected Double0D getBeta() {
			return (Double0D) this.fixedParms[0];
		}
		
		protected Double0D getLambda() {
			return (Double0D) this.fixedParms[1];
		}
		
		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(5);
			this.chainParmNames.add("A");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.chainParmNames.add("Gamma");
			this.chainParmNames.add("A0");
			this.chainParmNames.add("M0");
			this.fixedParmNames = new ArrayList<String>(1);
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
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					4);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double1D.class);
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
		protected Double2D genVariate() throws ProbDistParmException {
			int h = this.getA0().size();
			int K = this.getA().size();
			int d = this.getA().numCols();
			Integer1D active = this.getZ().items();
			int Ka = active.size();
			Integer1D[] zK = new Integer1D[K];
			Double2D[] xK = new Double2D[K];
			int[] nK = new int[K];
			Integer2D[] bK = new Integer2D[K];
			Double1D[] be = new Double1D[K];
			int gamma[] = this.getGamma().value();
			Integer1D B;
			for (Integer0D k : active) {
				zK[k.value()] = this.getZ().which(k);
				xK[k.value()] = this.getSamplerData().getAll(zK[k.value()].value());
				nK[k.value()] = xK[k.value()].size();
				bK[k.value()] = this.getB().getAll(zK[k.value()].value());
				// \sum_j x_j
				be[k.value()] = xK[k.value()].sum();
				for (int i=0; i<h; i++) {
					if (gamma[i] == 1) {
						// \sum_j \sum_{i} b_{i,j} * a_{i,k}
						B = bK[k.value()].sum();
						be[k.value()] = be[k.value()].minus(B.mult(this.getA().get(k.value())));
					}	
				}
			}
			// 1. sample from matrix normal centred at current point
			MVNormalDist mvnDist;
			Double1D[] aProp = new Double1D[K];
			Double1D[] aCurr = (Double1D[]) this.getA0().toArray();
			for (Integer0D k : active) {
				mvnDist = new MVNormalDist(this.getA0().get(k.value()), this.getSg().get(k.value()).mult(1/Math.pow(nK[k.value()], 2)));
				aProp[k.value()] = mvnDist.variateFast();
			}	
			// 2a. compute (unnormalized) matrix T density of proposal times matrix gaussian density of old
			// 2b. compute (unnormalized) matrix T density of old times matrix gaussian density of proposal
			double logNum = 0;
			double logDen = 0;
			double exp = -(d*K + this.getLambda().value());
			Double1D ambK;
			for (Integer0D k : active) {
				logNum += exp*Math.log(aTFnNoexp(aProp, this.getM0(), this.getSg(), this.getBeta().value(), d, active, this.getLambda().value()));
				ambK = aCurr[k.value()].minus(be[k.value()].mult(1/nK[k.value()]));
				logNum += ambK.mult(this.getSg().get(k.value()).mult(1/Math.pow(nK[k.value()], 2))).mult(ambK);
				logDen += exp*Math.log(aTFnNoexp(aCurr, this.getM0(), this.getSg(), this.getBeta().value(), d, active, this.getLambda().value()));
				ambK = aProp[k.value()].minus(be[k.value()].mult(1/nK[k.value()]));
				logDen += ambK.mult(this.getSg().get(k.value()).mult(1/Math.pow(nK[k.value()], 2))).mult(ambK);
			}
			// 3. accept or reject
			double logP = Math.log(FLGFVModel.unifGen.nextDouble());
			if (logP < logNum - logDen) {
				return new Double2D(aProp);
			} else {
				return this.getA0();
			}
		}

		@Override
		protected double getDensity(Double2D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}
		
	// Neal's (1998) Algorithm 8
	class CPDist extends ProbDistMC<Integer1D> {
		public CPDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double1D postProb;
		private CategoricalDist catDist;
		private HashMap<Integer0D, MVNormalDist> mvnDists; // Shouldn't you be
		// using (and
		// hacking)
		// SeqProbDist for
		// this? FIXME

		private ProbDist<Double2D> baseSg;
		private ProbDist<Double2D> baseA;
		private ProbDist<Double1D> baseA0;
		
		public ProbDist<Double2D> getBaseSg() {
			return this.baseSg;
		}

		public void setBaseSg(ProbDist<Double2D> baseSg) {
			this.baseSg = baseSg;
		}

		public ProbDist<Double2D> getBaseA() {
			return this.baseA;
		}

		public void setBaseA(ProbDist<Double2D> baseA) {
			this.baseA = baseA;
		}

		public ProbDist<Double1D> getBaseA0() {
			return this.baseA0;
		}

		public void setBaseA0(ProbDist<Double1D> baseA0) {
			this.baseA0 = baseA0;
		}

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
				double[] zeroD = new double[FLGFVModel.this.dims];
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
		
		protected Double3D getSg() {
			return (Double3D) this.chainParms[0];
		}

		protected Double3D getA() {
			return (Double3D) this.chainParms[1];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[2];
		}

		protected Integer2D getB() {
			return (Integer2D) this.chainParms[3];
		}

		protected Double2D getM() {
			return (Double2D) this.chainParms[7];
		}

		protected Double1D getM0() {
			return (Double1D) this.chainParms[8];
		}

		protected Double0D getAl() {
			return (Double0D) this.chainParms[4];
		}

		protected Integer1D getGamma() {
			return (Integer1D) this.chainParms[5];
		}
		
		protected Double2D getA0() {
			return (Double2D) this.chainParms[6];
		}
		
		protected Double2D getPsi() {
			return (Double2D) this.fixedParms[0];
		}

		protected Double0D getKappa() {
			return (Double0D) this.fixedParms[1];
		}

		protected Double0D getLambda() {
			return (Double0D) this.fixedParms[2];
		}

		protected Double0D getBeta() {
			return (Double0D) this.fixedParms[3];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(7);
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("A");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.chainParmNames.add("al");
			this.chainParmNames.add("Gamma");
			this.chainParmNames.add("A0");
			this.chainParmNames.add("M");
			this.chainParmNames.add("M0");
			this.fixedParmNames = new ArrayList<String>(5);
			this.fixedParmNames.add("Psi");
			this.fixedParmNames.add("kappa");
			this.fixedParmNames.add("lambda");
			this.fixedParmNames.add("beta");
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
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					7);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					5);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			this.mvnDists = null;
		}

		@Override
		protected Integer1D genVariate() throws ProbDistParmException {
			// FIXME - This sampling scheme destructively modifies
			// A, Sg, Omega - HACK
			int N = this.getSamplerData().size();
			Integer1D Z = null;
			try {
				Z = (Integer1D) this.getZ().clone();
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace(); // FIXME
			}
			Integer0D unused = new Integer0D(-1);
			int h = this.getA().get(0).size();
			double[][] gm = new double[h][h];
			for (int q=0; q<h; q++) {
				gm[q][q] = this.getGamma().get(q).value();
			}
			Double2D gammat = new Double2D(gm);
			Double1D zero = new Double1D(new double[FLGFVModel.this.dims]);
			for (int i = 0; i < N; i++) {
				// Sample a new category if necessary, remove the current
				// category from consideration
				Integer1D active = Z.items();
				int maxActive = active.get(active.size() - 1).value();
				Integer0D zedI = Z.get(i);
				Z.set(i, unused);
				Integer1D zZedI = Z.which(zedI);
				int nZZedI = zZedI.size();
				if (nZZedI > 0) { // Sample new category
					Integer0D zedNew = Z.minNotIn(0, maxActive);
					this.resetPostMvn(zedNew);
					Double2D sgNew = this.getBaseSg().variateFast();
					// sigma[k] = sg[k]*lambda/(beta)
					Double2D aNew = FLGFVModel.rmvt(this.getM(), sgNew.mult(this.getLambda().value()/this.getBeta().value()), 2*this.getLambda().value());
					// Eek! What a silly waste of time
					for (int j=0; j<h; j++) {
						if (this.getGamma().get(j).value()==0) {
							aNew.set(j, zero);
						}
					}
					Double1D a0New = this.getBaseA0().variateFast(this.getM0(), sgNew);
					// Add it
					this.getSg().set(zedNew.value(), sgNew);
					this.getA().set(zedNew.value(), aNew);
					this.getA0().set(zedNew.value(), a0New);
					active.add(zedNew);
				}
				// Select a category for this point
				double[] logP = new double[active.size()]; // FIXME - Probably
				// you could speed
				// this up by a lot
				double maxLogP = -Double.MAX_VALUE; // FIXME??
				int ai = 0;
				for (Integer0D k : active) {
					Integer1D zK = Z.which(k);
					Double1D b = this.getB().get(i).mult(gammat);
					Double1D aTB = b.mult(this.getA().get(k.value()));
					Double1D xIA = this.getSamplerData().get(i).minus(aTB).minus(this.getA0().get(k.value()));
					int nZK = zK.size();
					this.setUpPostMvn(k);
					double logPrior;
					if (nZK == 0) {
						logPrior = Math.log(this.getAl().value());
					} else {
						logPrior = Math.log(nZK);
					}
					logP[ai] = this.mvnLogDensity(k, xIA) + logPrior;
					if (logP[ai] > maxLogP) {
						maxLogP = logP[ai];
					}
					ai++;
				}
				Double1D p = (new Double1D(logP)).minus(maxLogP).exp();
				this.setPostProb(p);
				int az = this.catVariate().value();
				Z.set(i, active.get(az));
			}
			return Z;
		}

		@Override
		protected double getDensity(Integer1D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PostXDist extends ProbDistMC<Double0D> {
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

	class PostAlDist extends ProbDistMC<Double0D> {
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

	class PostBDist extends ProbDistMC<Integer2D> {
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
	
	class PostMDist extends ProbDistMC<Double2D> {
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
			this.chainParmNames.add("gamma");
			this.chainParmNames.add("Z");
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
			int d = this.getM().numCols();
			Integer1D active = this.getZ().items();
			int K = active.size();
			int[] gamma = this.getGamma().value();
			int h = gamma.length;
			// 1. sample from matrix normal centred at current value
			Double1D[] mNew = null;
			try {
				mNew = (Double1D[]) ((Double2D) this.getM().clone()).toArray();
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			MVNormalDist mvnDist = null;
			for (int i=0; i<h; i++) {
				if (gamma[i] == 1) {
					mvnDist = new MVNormalDist(this.getM().get(i), this.getS());
					mNew[i] = mvnDist.variateFast();
				}
			}
			// 2a. compute (unnormalized) matrix T density of proposal times matrix gaussian density of old
			// 2b. compute (unnormalized) matrix T density of old times matrix gaussian density of proposal
			Double1D[][] aK = new Double1D[h][this.getSg().size()];
			for (Integer0D k : active) {
				for (int i=0; i<h; i++) {
					if (gamma[i] == 1) {
						aK[i][k.value()] = this.getA().get(k.value()).get(i);
					}
				}
			}
			double exp = -(d*K + this.getLambda().value());
			double logNum = 0;
			double logDen = 0;
			for (int i=0; i<h; i++) {
				if (gamma[i] == 1) {
					logNum += exp*Math.log(aTFnNoexp(aK[i], mNew[i], this.getSg(), this.getBeta().value(), d, active, this.getLambda().value()));
					logNum += this.getM().get(i).mult(this.getS()).mult(this.getM().get(i));
					logDen += exp*Math.log(aTFnNoexp(aK[i], this.getM().get(i), this.getSg(), this.getBeta().value(), d, active, this.getLambda().value()));
					logDen += mNew[i].mult(this.getS()).mult(mNew[i]);
				}
			}
			// 3. accept or reject
			double logP = Math.log(FLGFVModel.unifGen.nextDouble());
			if (logP < logNum - logDen) {
				return new Double2D(mNew);
			} else {
				return this.getM();
			}
		}

		@Override
		protected double getDensity(Double2D pt) throws ProbDistParmException {
			// TODO Auto-generated method stub
			return 0;
		}
	}
	
	class PostM0Dist extends ProbDistMC<Double1D> {
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
			this.chainParmNames.add("A0");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("M");
			this.chainParmNames.add("gamma");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("M0");
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
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double1D.class);
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
			int d = this.getM().numCols();
			Integer1D active = this.getZ().items();
			int K = active.size();
			// 1. sample from matrix normal centred at current value
			Double1D mNew = null;
			MVNormalDist mvnDist = new MVNormalDist(this.getM0(), this.getS());
			mNew = mvnDist.variateFast();
			// 2a. compute (unnormalized) matrix T density of proposal times matrix gaussian density of old
			// 2b. compute (unnormalized) matrix T density of old times matrix gaussian density of proposal
			double exp = -(d*K + this.getLambda().value());
			double logNum = exp*Math.log(aTFnNoexp((Double1D[]) this.getA0().toArray(), mNew, this.getSg(), this.getBeta().value(), d, active, this.getLambda().value()));
			logNum += this.getM0().mult(this.getS()).mult(this.getM0());
			double logDen = exp*Math.log(aTFnNoexp((Double1D[]) this.getA0().toArray(), this.getM0(), this.getSg(), this.getBeta().value(), d, active, this.getLambda().value()));
			logDen += mNew.mult(this.getS()).mult(mNew);
			// 3. accept or reject
			double logP = Math.log(FLGFVModel.unifGen.nextDouble());
			if (logP < logNum - logDen) {
				return mNew;
			} else {
				return this.getM0();
			}
		}

		@Override
		protected double getDensity(Double1D pt) throws ProbDistParmException {
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
					CPDist postZ = (CPDist) rvZ.getPosterior();
					apost += postZ.getBaseSg().logDensity(Sg);
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
	
	public FLGFVModel() {
		super();
	}

	public FLGFVModel(Map<String, Numeric<? extends Numeric<?>>> hypers,
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
		PostMDist postM = new PostMDist((Double2D) this.getHyper("S"), (Double0D) this.getHyper("beta"), (Double0D) this.getHyper("lambda"));
		Double2D defaultM = (Double2D) init.get("M"); // FIXME
		RandomVar<Double2D> rvM = new RandomVar<Double2D>("M", priorM, postM,
				defaultM);
		this.params.put("M", rvM);

		// A0
		PriorA0Dist priorA0 = new PriorA0Dist();
		PostA0Dist postA0 = new PostA0Dist((Double0D) this.getHyper("beta"), (Double0D) this.getHyper("lambda"));
		Double2D defaultA0 = (Double2D) init.get("A0"); // FIXME
		RandomVar<Double2D> rvA0 = new RandomVar<Double2D>("A0", priorA0, postA0,
				defaultA0);
		this.params.put("A0", rvA0);
		
		// A
		PriorADist priorA = new PriorADist();
		PostADist postA = new PostADist((Double0D) this.getHyper("lambda"), (Double0D) this.getHyper("beta"));
		Double3D defaultA = (Double3D) init.get("A"); // FIXME
		RandomVar<Double3D> rvA = new RandomVar<Double3D>("A", priorA, postA,
				defaultA);
		this.params.put("A", rvA);

		// Z
		ProbDist<Integer1D> priorZ = new Deterministic<Integer1D>(init.get("Z"));
		CPDist postZ = new CPDist((Double2D) this.getHyper("Psi"),
				(Double0D) this.getHyper("kappa"),
				(Double0D) this.getHyper("lambda"), (Double0D) this.getHyper("beta")); // FIXME
		// --
		// Hack!
		postZ.setBaseSg(new InverseWishartDist(this.getHyper("Psi"), this
				.getHyper("kappa"))); // FIXME -- Eek! Hack!
		postZ.setBaseA(new MatrixNormalDist()); // FIXME -- Eek! Hack!
		postZ.setBaseA0(new MVNormalDist()); // FIXME -- Eek! Hack!
		Integer1D defaultZ = (Integer1D) init.get("Z");
		RandomVar<Integer1D> rvZ = new RandomVar<Integer1D>("Z", priorZ, postZ,
				defaultZ);
		this.params.put("Z", rvZ); // FIXME
		
		// Gamma
		ProbDist<Integer1D> priorGamma = new Deterministic<Integer1D>(init.get("Gamma"));
		PostGammaDist postGamma = new PostGammaDist(hypers.get("p"), hypers.get("lambda"), hypers.get("beta"));
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
		ProbDistMC<Double0D> postAl = new PostAlDist(
				(Double0D) this.getHyper("ala"),
				(Double0D) this.getHyper("alb"));
		Double0D defaultAl = (Double0D) init.get("al"); // FIXME
		RandomVar<Double0D> rvAl = new RandomVar<Double0D>("al", priorAl,
				postAl, defaultAl);
		this.params.put("al", rvAl); // FIXME
	}

	@Override
	public ChainLink getInitialLink() {
		return new ChainLink(this.getParam("Z"), this.getParam("B"), this.getParam("Sg"),
				this.getParam("Gamma"), this.getParam("A"), this.getParam("x"), this.getParam("A0"),
				this.getParam("al"), this.getParam("M"), this.getParam("M0"));
	}
}