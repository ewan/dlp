package org.jqgibbs.models;

import java.util.ArrayList;
import java.util.Arrays;
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
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.mathstat.RandomVar;
import org.jqgibbs.mathstat.probdist.BetaDist;
import org.jqgibbs.mathstat.probdist.CategoricalDistInitializeByP;
import org.jqgibbs.mathstat.probdist.GammaDist;
import org.jqgibbs.mathstat.probdist.InverseWishartDist;
import org.jqgibbs.mathstat.probdist.MVNormalDist;
import org.jqgibbs.mathstat.probdist.MatrixNormalDist;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistInitializeByChain;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;

import cern.jet.stat.Gamma;

public class FLGFDModel extends Model {

	private static double gammaD(double x, int d) { // FIXME
		if (d == 1) {
			return Gamma.gamma(x);
		} else { // FIXME check d0
			return Math.pow(Math.PI, (d - 1) / 2) * FLGFDModel.gammaD(x, d - 1)
					* Gamma.gamma(x + (1 - d) / 2);
		}
	}

	class PostOmegaDist extends ProbDistInitializeByChain<Double2D> {
		private Double2D postPhi;
		private Double0D postLambda;
		private InverseWishartDist iwDist;
		private ProbDist<Double2D> baseSg;
		
		// Fixed parameters
		private Double2D phi;
		private Double0D lambda;
		
		// Chain parameters
		private Integer1D Z;
		private Double2D M;
		private Double3D Sg;
		private Double3D A;

		public PostOmegaDist(Double2D phi, Double0D lambda) throws ProbDistParmException {
			this.phi = phi;
			this.lambda = lambda;
		}
		
		protected void initializeChainParms(ChainLink l) {
			this.Z = (Integer1D) l.get("Z").getNumericValue();
			this.M = (Double2D) l.get("M").getNumericValue();
			this.Sg = (Double3D) l.get("Sg").getNumericValue();
			this.A = (Double3D) l.get("A").getNumericValue();
		}
		
		public ProbDist<Double2D> getBaseSg() {
			return this.baseSg;
		}

		public void setBaseSg(ProbDist<Double2D> baseSg) {
			this.baseSg = baseSg;
		}

		protected Double2D getPostPhi() {
			return this.postPhi;
		}

		protected Double0D getPostLambda() {
			return this.postLambda;
		}

		protected void setPostPhi(Double2D postPhi) {
			this.postPhi = postPhi;
		}

		protected void setPostLambda(Double0D postLambda) {
			this.postLambda = postLambda;
		}

		protected Double2D getPhi() {
			return phi;
		}

		protected Double0D getLambda() {
			return lambda;
		}

		protected Integer1D getZ() {
			return Z;
		}

		protected Double2D getM() {
			return M;
		}

		protected Double3D getSg() {
			return Sg;
		}

		protected Double3D getA() {
			return A;
		}

		private Double2D iwVariate() throws ProbDistParmException {
			if (this.iwDist == null) {
				this.iwDist = new InverseWishartDist(this.getPostPhi(), this
						.getPostLambda());
				return this.iwDist.variate();
			} else {
				return this.iwDist.variate(this.getPostPhi(), this
						.getPostLambda());
			}
		}

		protected void setUpFromChainParms() {
			Double2D postPhi = this.getPhi();
			Double0D postLambda = this.getLambda();
			Integer1D active = this.getZ().items();
			int d = FLGFDModel.this.dims;
			for (Integer0D k : active) {
				Double2D AM = this.getA().get(k.value()).minus(this.getM());
				Double2D sgInv = this.getSg().get(k.value()).inverse();
				postLambda = postLambda.plus(d);
				postPhi = postPhi.plus(AM.mult(sgInv).mult(AM.transpose()));
			}
			this.setPostPhi(postPhi);
			this.setPostLambda(postLambda);
		}

		@Override
		protected Double2D genVariate() throws ProbDistParmException {
			boolean goodSample = false;
			Double2D o = null;
			while (!goodSample) {
				o = this.iwVariate();
				goodSample = true;
				Integer1D active = this.getZ().items();
				for (Integer0D k : active) {
					Double2D sg = this.getSg().get(k.value());
					Double2D sgKOmega = sg.kron(o);
					if (!sgKOmega.isWellConditioned()) {
						goodSample = false;
						System.err.println("Warning: threw out IW sample for Omega that had ill-conditioned Kronecker product with one or more current Sgs");
						break;
					}
				}
				if (goodSample) {
					Double2D sg = this.getBaseSg().variateFast();
					Double2D sgKOmega = sg.kron(o);
					if (!sgKOmega.isWellConditioned()) {
						goodSample = false;
						System.err.println("Warning: threw out IW sample for Omega that had ill-conditioned Kronecker product with a sample from the Sg prior");
					}
				}
			}
			return o;
		}

		@Override
		protected double getDensity(Double2D pt) throws ProbDistParmException {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PriorSgDist extends ProbDistInitializeByChain<Double3D> {
		// Fixed parameters
		private Double2D psi;
		private Double0D kappa;
		
		// Chain parameters
		private Integer1D Z;
		
		private List<Double2D> repPsi;
		private List<Double0D> repKappa;

		public PriorSgDist(Double2D psi, Double0D kappa) throws ProbDistParmException {
			this.psi = psi;
			this.kappa = kappa;
		}
		
		protected void initializeChainParms(ChainLink l) {
			this.Z = (Integer1D) l.get("Z").getNumericValue();
		}
		
		protected void setUpFromChainParms() {
			// empty
		}

		private List<Double2D> getRepPsi() {
			if (this.repPsi == null) {
				this.repPsi = new ArrayList<Double2D>();
				Integer1D active = this.getZ().items();
				for (int i = 0; i < active.size(); i++) {
					this.repPsi.add(this.getPsi());
				}
			}
			return this.repPsi;
		}

		private List<Double0D> getRepKappa() {
			if (this.repKappa == null) {
				this.repKappa = new ArrayList<Double0D>();
				Integer1D active = this.getZ().items();
				for (int i = 0; i < active.size(); i++) {
					this.repKappa.add(this.getKappa());
				}
			}
			return this.repKappa;
		}

		private Double2D getPsi() {
			return psi;
		}

		private Double0D getKappa() {
			return kappa;
		}

		private Integer1D getZ() {
			return Z;
		}

		private Double3D iwVariates() throws ProbDistParmException {
			return InverseWishartDist.variates(this.getRepPsi(), this.getRepKappa());
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
		// Fixed parameters
		private Double2D Psi;
		private Double0D kappa;
		
		// Chain parameters
		private Integer1D Z;
		private Double2D B;
		private Double2D Omega;
		private Double3D A;
		private Double2D M;
		
		private List<Double2D> postPsi;
		private List<Double0D> postKappa;

		public PostSgDist(Double2D Psi, Double0D kappa) throws ProbDistParmException {
			this.Psi = Psi;
			this.kappa = kappa;
		}
		
		protected void initializeChainParms(ChainLink l) {
			this.Z = (Integer1D) l.get("Z").getNumericValue();
			this.B = (Double2D) l.get("B").getNumericValue();
			this.Omega = (Double2D) l.get("Omega").getNumericValue();
			this.A = (Double3D) l.get("A").getNumericValue();
			this.M = (Double2D) l.get("M").getNumericValue();
		}

		protected List<Double2D> getPostPsi() {
			return this.postPsi;
		}

		protected List<Double0D> getPostKappa() {
			return this.postKappa;
		}

		protected void setPostPsi(List<Double2D> postPsi) {
			this.postPsi = postPsi;
		}

		protected void setPostKappa(List<Double0D> postKappa) {
			this.postKappa = postKappa;
		}

		protected Double2D getPsi() {
			return this.Psi;
		}

		protected Double0D getKappa() {
			return kappa;
		}

		protected Integer1D getZ() {
			return Z;
		}

		protected Double2D getB() {
			return B;
		}

		protected Double2D getOmega() {
			return Omega;
		}

		protected Double3D getA() {
			return A;
		}

		protected Double2D getM() {
			return M;
		}

		private Double3D iwVariates() throws ProbDistParmException {
			return InverseWishartDist.variates(this.getPostPsi(), this.getPostKappa());
		}
		
		protected void setUpFromChainParms() {
			this.setPostPsi(new ArrayList<Double2D>());
			this.setPostKappa(new ArrayList<Double0D>());
			Integer1D active = this.getZ().items();
			int h = this.getA().size();
			Double2D omegaInv = this.getOmega().inverse();
			for (Integer0D k : active) {
				Integer1D zK = this.getZ().which(k);
				Double2D xK = this.getSamplerData().getAll(zK.value());
				Integer0D nK = new Integer0D(xK.size());

				this.getPostKappa().add(this.getKappa().plus(nK));

				Double2D B = this.getB().getAll(zK.value());
				Double2D BTB = B.transpose().mult(B);
				Double2D MTOmegaInv = this.getM().transpose().mult(omegaInv);
				Double2D MTOmegaInvM = MTOmegaInv.mult(this.getM());
				Double2D XTB = B.transpose().mult(xK).transpose();
				Double2D Mp = XTB.plus(MTOmegaInv);
				Double2D postOmega = omegaInv.plus(BTB).inverse();
				Double2D MpPOmegaMpT = Mp.mult(postOmega).mult(Mp.transpose());
				Double2D XTX = xK.transpose().mult(xK);
				this.getPostPsi().add(
						this.getPsi().plus(XTX).plus(MTOmegaInvM).minus(
								MpPOmegaMpT));
			}
		}

		@Override
		protected Double3D genVariate() throws ProbDistParmException {
			int d = FLGFDModel.this.dims;
			double[][][] v = new double[this.getA().size()][d][d];
			Double3D activeUpdates = this.iwVariates();
			int i = 0;
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], d);
				i++;
			}
			return new Double3D(v);
		}

		@Override
		protected double getDensity(Double3D pt) throws ProbDistParmException {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PriorADist extends ProbDistInitializeByChain<Double3D> {
		// Chain parameters
		private Double2D Omega;
		private Double3D Sg;
		private Double2D M;
		private Integer1D Z;
		private Double2D B;
		
		public PriorADist() throws ProbDistParmException {
			// Empty (no fixed parameters)
		}
		
		protected void initializeChainParms(ChainLink l) {
			this.Omega = (Double2D) l.get("Omega").getNumericValue();
			this.Sg = (Double3D) l.get("Sg").getNumericValue();
			this.M = (Double2D) l.get("M").getNumericValue();
			this.Z = (Integer1D) l.get("Z").getNumericValue();
			this.B = (Double2D) l.get("B").getNumericValue();
		}
		
		protected void setUpFromChainParms() {
			// empty
		}

		private Double2D getOmega() {
			return Omega;
		}

		private Double3D getActiveSg() {
			Integer1D active = this.getZ().items();
			return Sg.getAll(active.value());
		}

		private Double3D getRepM() {
			Integer1D active = this.getZ().items();
			Double3D repM = new Double3D(active.size(), this.getB().numCols(),
					FLGFDModel.this.dims);
			for (int i = 0; i < active.size(); i++) {
				repM.add(M);
			}
			return repM;
		}

		private Integer1D getZ() {
			return Z;
		}

		private Double2D getB() {
			return B;
		}

		private Double3D matnVariates() throws ProbDistParmException {
			// This needs to be fixed!!!
			List<Double2D> repM = new ArrayList<Double2D>();
			for(Double2D m : getRepM()) {
				repM.add(m);
			}
			List<Double2D> Omega = new ArrayList<Double2D>();
			Omega.add(getOmega());
			return MatrixNormalDist.variates(repM, getActiveSg(), Omega);
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
		private List<Double2D> postMs;
		private List<Double2D> postOmegas;
		
		// Chain parameters
		private Double2D Omega;
		private Double3D Sg;
		private Double2D M;
		private Integer1D Z;
		private Double2D B;

		public PostADist() throws ProbDistParmException {
			// Empty (no fixed parameters)
		}
		
		protected void initializeChainParms(ChainLink l) {
			this.Omega = (Double2D) l.get("Omega").getNumericValue();
			this.Sg = (Double3D) l.get("Sg").getNumericValue();
			this.M = (Double2D) l.get("M").getNumericValue();
			this.Z = (Integer1D) l.get("Z").getNumericValue();
			this.B = (Double2D) l.get("B").getNumericValue();
		}

		private Double3D matnVariates() throws ProbDistParmException {		
			Integer1D active = this.getZ().items();
			Double3D sgActive = this.getSg().getAll(active.value());
			return MatrixNormalDist.variates(getPostM(), sgActive, getPostOmegas());
//			if (this.matnDists == null) {
//				this.matnDists = new SeqMatrixNormalDist(new Double2D(this.getPostM().toArray(d2d)),
//						sgActive, new Double2D(this.getPostOmegas().toArray(d2d)));
//				// FIXME - Will this work if all arguments are not
//				// Lists?
//				return this.matnDists.variate();
//			} else {
//				return this.matnDists.variate(new Double2D(this.getPostM().toArray(d2d)), sgActive, new Double2D(this
//						.getPostOmegas().toArray(d2d)));
//			}
		}

		protected List<Double2D> getPostM() {
			return this.postMs;
		}

		protected List<Double2D> getPostOmegas() {
			return this.postOmegas;
		}

		protected void setPostM(List<Double2D> postMs) {
			this.postMs = postMs;
		}

		protected void setPostOmegas(List<Double2D> postOmegas) {
			this.postOmegas = postOmegas;
		}

		protected Double2D getOmega() {
			return Omega;
		}

		protected Double3D getSg() {
			return Sg;
		}

		protected Double2D getM() {
			return M;
		}

		protected Integer1D getZ() {
			return Z;
		}

		protected Double2D getB() {
			return B;
		}
		
		protected void setUpFromChainParms() {
			this.setPostOmegas(new ArrayList<Double2D>());
			this.setPostM(new ArrayList<Double2D>());
			Integer1D active = this.getZ().items();
			Double2D omegaInv = this.getOmega().inverse();
			for (Integer0D k : active) {
				Integer1D zK = this.getZ().which(k);
				Double2D xK = this.getSamplerData().getAll(zK.value());
				Double2D B = this.getB().getAll(zK.value());
				Double2D BTB = B.transpose().mult(B);
				Double2D postOmegaInv = omegaInv.plus(BTB);
				Double2D postOmega = postOmegaInv.inverse();
				this.getPostOmegas().add(postOmega);
				Double2D BTX = B.transpose().mult(xK);
				Double2D omegaInvM = omegaInv.mult(this.getM());
				Double2D postM = postOmega.mult(BTX.plus(omegaInvM));
				this.getPostM().add(postM);
			}
		}

		@Override
		protected Double3D genVariate() throws ProbDistParmException {
			int K = this.getSg().size();
			int h = this.getB().numCols();
			int d = FLGFDModel.this.dims;
			double[][][] v = new double[K][h][d];
			boolean goodSample = false;
			while (!goodSample) {
				Double3D activeUpdates = this.matnVariates();
				goodSample = true;
				int i = 0;
				Integer1D active = this.getZ().items();
				for (Integer0D k : active) {
					v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], h); // FIXME??
					// Check
					// dims
					for (int m = 0; m < h; m++) {
						for (int n = 0; n < d; n++) {
							if (Double.isInfinite(v[k.value()][m][n])
									|| Double.isNaN(v[k.value()][m][n])) {
								System.err
										.println("Warning: bad sample for A: "
												+ activeUpdates.get(i));
								System.err.println("Sg: "
										+ this.getSg().get(k.value()));
								// System.err.println("Omega: " +
								// this.getOmega().get(k.value()));
								System.err.println("M: " + this.getM());
								goodSample = false;
								break;
							}
						}
						if (!goodSample) {
							break;
						}
					}
					if (!goodSample) {
						break;
					}
					i++;
				}
			}
			return new Double3D(v);
		}

		@Override
		protected double getDensity(Double3D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class CRPDist extends ProbDistInitializeByChain<Integer1D> {
		// Fixed parameters
		private Double2D Psi;
		private Double0D kappa;
		private Double2D Phi;
		private Double0D lambda;
		
		// Chain parameters
		private Double2D Omega;
		private Double3D Sg;
		private Double3D A;
		private Integer1D Z;
		private Double2D B;
		private Double2D M;
		private Double0D al;
		
		public CRPDist(Double2D Psi, Double0D kappa, Double2D Phi, Double0D lambda) throws ProbDistParmException {
			this.Psi = Psi;
			this.kappa = kappa;
			this.Phi = Phi;
			this.lambda = lambda;
		}
		
		protected void initializeChainParms(ChainLink l) {
			this.Omega = (Double2D) l.get("Omega").getNumericValue();
			this.Sg = (Double3D) l.get("Sg").getNumericValue();
			this.A = (Double3D) l.get("A").getNumericValue();
			this.Z = (Integer1D) l.get("Z").getNumericValue();
			this.B = (Double2D) l.get("B").getNumericValue();
			this.M = (Double2D) l.get("M").getNumericValue();
			this.al = (Double0D) l.get("al").getNumericValue();
		}

		protected Double1D postProb;
		private CategoricalDistInitializeByP catDist;
		private HashMap<Integer0D, MVNormalDist> mvnDists; // Shouldn't you be
		// using (and
		// hacking)
		// SeqProbDist for
		// this? FIXME
		private InverseWishartDist sgDist;
		private MatrixNormalDist aDist;

		private Double2D omegaInv;

		private ProbDist<Double2D> baseSg;
		private ProbDist<Double2D> baseA;

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

		protected Integer0D catVariate() throws ProbDistParmException {
			if (this.catDist == null) {
				this.catDist = new CategoricalDistInitializeByP(this
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
				double[] zeroD = new double[FLGFDModel.this.dims];
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

		protected double marginalNew(Double1D pt, Double1D b) {
			Double2D omegaN = this.getOmegaInv().plus(b.outer(b)).inverse();
			Double2D mH = this.getM().transpose().mult(this.getOmegaInv());
			Double2D mHX = mH.plus(pt.outer(b));
			Double2D mO = this.getM().transpose().mult(this.getOmegaInv())
					.mult(this.getM());
			Double2D mHXO = mHX.mult(omegaN).mult(mHX.transpose());
			Double2D psiN = this.getPsi().plus(
					pt.outer(pt).plus(mO).minus(mHXO));
			int d = this.getSamplerData().numCols();
			int h = this.getB().numCols();
			double kN = this.getKappa().value() + 1;
			double lmd = Math.log(FLGFDModel.gammaD(kN / 2, d));
			lmd -= Math.log(2) * d / 2;
			lmd -= Math.log(2 * Math.PI) * d / 2;
			lmd += Math.log(omegaN.det()) * d / 2;
			lmd -= Math.log(this.getOmega().det()) * d / 2;
			lmd -= Math.log(psiN.det()) * kN / 2;
			lmd += Math.log(this.getPsi().det()) * this.getKappa().value()/2;
			lmd -= Math.log(FLGFDModel.gammaD(this.getKappa().value()/2, d));
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

		protected List<Double2D> getVariateDPPost(Double1D pt,
				Double1D b) throws ProbDistParmException {
			Double2D omegaN = this.getOmega().inverse().plus(b.outer(b)).inverse();
			Double2D mH = this.getM().transpose().mult( this.getOmega().inverse());
			Double2D mHX = mH.plus(pt.outer(b));
			Double2D mHXO = mHX.mult(omegaN).mult(mHX.transpose());
			Double2D mO = this.getM().transpose().mult(
					this.getOmega().inverse()).mult(this.getM());
			Double2D psiN = this.getPsi().plus(
					pt.outer(pt).plus(mO).minus(mHXO));
			double kN = this.getKappa().value() + 1;
			Double2D sgNew = null;
			Double2D aNew = null;
			boolean goodSample = false;
			while (!goodSample) {
				goodSample = true;
				sgNew = this.sgVariate(psiN, kN);
				aNew = this.aVariate(omegaN.mult(mHX.transpose()), sgNew, omegaN);
				double[][] a = aNew.value();
				for (int m = 0; m < aNew.numRows(); m++) {
					for (int n = 0; n < aNew.numCols(); n++) {
						if (Double.isInfinite(a[m][n]) || Double.isNaN(a[m][n])) {
							System.err
									.println("Warning: bad sample for new A: "
											+ aNew);
							System.err.println("New Sg: " + sgNew);
							System.err.println("Omega: " + omegaN);
							System.err.println("mean: "
									+ omegaN.mult(mHX.transpose()));
							goodSample = false;
							break;
						}
					}
					if (!goodSample) {
						break;
					}
				}
			}
			return new ArrayList<Double2D>(Arrays.asList(sgNew, aNew));
		}

		protected Double2D getOmega() {
			return Omega;
		}

		protected Double3D getSg() {
			return Sg;
		}

		protected Double3D getA() {
			return A;
		}

		protected Integer1D getZ() {
			return Z;
		}

		protected Double2D getB() {
			return B;
		}

		protected Double2D getM() {
			return M;
		}

		protected Double0D getAl() {
			return al;
		}

		protected Double2D getPsi() {
			return Psi;
		}

		protected Double0D getKappa() {
			return kappa;
		}

		// protected Double2D getPhi() {
		public Double2D getPhi() { // FIXME
			return Phi;
		}

		protected Double0D getLambda() {
			return lambda;
		}

		protected Double2D getOmegaInv() {
			return this.omegaInv;
		}
		
		protected void setUpFromChainParms() {
			this.mvnDists = null;
			this.omegaInv = this.getOmega().inverse();
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
			for (int i = 0; i < N; i++) {
				Double1D Xi = this.getSamplerData().get(i);
				Integer1D active = Z.items();
				Z.set(i, unused);
				double[] logP = new double[active.size() + 1];
				double maxLogP = -Double.MAX_VALUE;
				int ai = 0;
				// Existing categories
				Double1D b = this.getB().get(i);
				for (Integer0D k : active) {
					Integer1D zK = Z.which(k);
					int nZK = zK.size();
					if (nZK > 0) {
						Double1D aTB = b.mult(this.getA().get(k.value()));
						Double1D xIA = this.getSamplerData().get(i).minus(aTB);
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
				logP[newc] = this.marginalNew(Xi, b) + Math.log(this.getAl().value());
				if (logP[ai] == Double.POSITIVE_INFINITY) {
					logP[ai] = Double.MAX_VALUE;
				}
				if (logP[newc] > maxLogP) {
					maxLogP = logP[newc];
				}
				// Select a category for this point
				Double1D p = (new Double1D(logP)).minus(maxLogP).exp();
				this.setPostProb(p);
				if(i < 8) {
					System.out.println(p);
				}
				int az = this.catVariate().value();
				// Sample new category if necessary
				Integer0D zedNew = null;
				if (az == newc) {
					int maxActive = active.get(active.size() - 1).value();
					zedNew = Z.minNotIn(0, maxActive);
					this.resetPostMvn(zedNew);
					List<Double2D> newcParms = this.getVariateDPPost(Xi, b);
					Double2D newSg = newcParms.get(0);
					Double2D newA = newcParms.get(1);
					this.getA().set(zedNew.value(), newA);
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

	class PostMDist extends ProbDistInitializeByChain<Double2D> {
		protected MVNormalDist mvnDist;
		protected Double1D vecW;
		protected Double2D SInv;
		protected Double1D postVecMu;
		protected Double2D postKronSg;
		
		// Fixed parameters
		private Double2D W;
		private Double2D S;
		private Double2D Phi;
		
		// Chain parameters
		private Integer1D Z;
		private Double2D B;
		private Double2D Omega;
		private Double3D Sg;
		private Double3D A;

		public PostMDist(Double2D W, Double2D S) throws ProbDistParmException {
			this.W = W;
			this.S = S;
		}
		
		protected void initializeChainParms(ChainLink l) {
			this.Z = (Integer1D) l.get("Z").getNumericValue();
			this.B = (Double2D) l.get("B").getNumericValue();
			this.Omega = (Double2D) l.get("Omega").getNumericValue();
			this.Sg = (Double3D) l.get("Sg").getNumericValue();
			this.A = (Double3D) l.get("A").getNumericValue();
		}

		protected Integer1D getZ() {
			return Z;
		}

		protected Double2D getB() {
			return B;
		}

		protected Double2D getOmega() {
			return Omega;
		}

		protected Double3D getSg() {
			return Sg;
		}

		protected Double3D getA() {
			return A;
		}

		protected Double2D getPhi() {
			return (Double2D) Phi; // when does this get set?
		}

		protected Double2D getSInv() {
			if (this.SInv == null) {
				this.SInv = S.inverse();
			}
			return this.SInv;
		}

		protected Double1D getVecW() {
			if (this.vecW == null) {
				this.vecW = W.colVec();
			}
			return this.vecW;
		}

		protected Double2D matnVariate() throws ProbDistParmException {
			int h = this.getA().numRows();
			Double1D vecVariate;
			if (this.mvnDist == null) {
				this.mvnDist = new MVNormalDist(this.getPostVecMu(), this
						.getPostKronSg());
				vecVariate = this.mvnDist.variate();

			} else {
				vecVariate = this.mvnDist.variate(this.getPostVecMu(), this
						.getPostKronSg());
			}
			return vecVariate.toDouble2D(h);
		}
		
		protected void setUpFromChainParms() {
			int h = this.getB().numCols();
			Double2D sInvKIdh = this.getSInv().kron(Double2D.ident(h));
			Double2D postKronSgInv = sInvKIdh;
			Double1D postVecMup = sInvKIdh.mult(this.getVecW());
			Integer1D active = this.getZ().items();
			Double2D omegaInv = this.getOmega().inverse();
			// Double2D omegaInv = this.getPhi().inverse();
			for (Integer0D k : active) {
				Double2D sgInv = this.getSg().get(k.value()).inverse();
				Double2D sgInvKOmegaInv = sgInv.kron(omegaInv);
				postKronSgInv = postKronSgInv.plus(sgInvKOmegaInv);
				Double1D vecA = this.getA().get(k.value()).colVec();
				postVecMup = postVecMup.plus(sgInvKOmegaInv.mult(vecA));
			}
			Double2D postKronSg = postKronSgInv.inverse();
			this.setPostKronSg(postKronSg);
			this.setPostVecMu(postKronSg.mult(postVecMup));
		}

		protected Double1D getPostVecMu() {
			return this.postVecMu;
		}

		protected void setPostVecMu(Double1D mu) {
			this.postVecMu = mu;
		}

		protected Double2D getPostKronSg() {
			return this.postKronSg;
		}

		protected void setPostKronSg(Double2D sg) {
			this.postKronSg = sg;
		}

		@Override
		protected Double2D genVariate() throws ProbDistParmException {
			return this.matnVariate();
		}

		@Override
		protected double getDensity(Double2D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PostXDist extends ProbDistInitializeByChain<Double0D> {
		private BetaDist betaDist;
		private Double0D postXA;
		private Double0D postXB;
		
		// Chain parameters
		private Double0D al;

		public PostXDist() throws ProbDistParmException {
			// Empty (no fixed parameters)
		}
		
		protected void initializeChainParms(ChainLink l) {
			this.al = (Double0D) l.get("al").getNumericValue();
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
			return al;
		}

		
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
		private CategoricalDistInitializeByP catDist;
		private Double0D postAlA;
		private Double0D postAlB;
		private Double1D postMix;
		
		// Fixed parameters
		private Double0D ala;
		private Double0D alb;
		
		// Chain parameters
		private Double0D x;
		private Integer1D Z;

		public PostAlDist(Double0D ala, Double0D alb) throws ProbDistParmException {
			this.ala = ala;
			this.alb = alb;
		}
		
		protected void initializeChainParms(ChainLink l) {
			this.x = (Double0D) l.get("x").getNumericValue();
			this.Z = (Integer1D) l.get("Z").getNumericValue();
		}

		private Double0D gammaVariate() throws ProbDistParmException {
			if (this.gammaDist == null) {
				this.gammaDist = new GammaDist(this.getPostAlA(), this
						.getPostAlB());
				return this.gammaDist.variate();
			} else {
				return this.gammaDist.variate(this.getPostAlA(), this
						.getPostAlB());
			}
		}

		private Integer0D catVariate() throws ProbDistParmException {
			if (this.catDist == null) {
				this.catDist = new CategoricalDistInitializeByP(this
						.getPostMix());
				return this.catDist.variate();
			} else {
				return this.catDist.variate(this.getPostMix());
			}
		}

		protected Double0D getX() {
			return x;
		}

		protected Integer1D getZ() {
			return Z;
		}

		protected Double0D getAlA() {
			return ala;
		}

		protected Double0D getAlB() {
			return alb;
		}
		
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

	class PostBDist extends ProbDistInitializeByChain<Double2D> {
		// Chain parameters
		private Double2D B;
		
		public PostBDist() throws ProbDistParmException {
			// Empty (no fixed parameters)
		}
		
		protected void initializeChainParms(ChainLink l) {
			this.B = (Double2D) l.get("B").getNumericValue();
		}
		
		protected void setUpFromChainParms() {
			// empty
		}
		
		private Double2D getB() {
			return B;
		}

		@Override
		protected Double2D genVariate() throws ProbDistParmException {
			return this.getB();
		}

		@Override
		protected double getDensity(Double2D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	@SuppressWarnings("unchecked")
	public static ChainLink pointEstimate(List<ChainLink> c, Double2D d)
			throws ProbDistParmException {
		// MAP
		double max = Double.NEGATIVE_INFINITY;
		int argmax = c.size() - 1;
		for (int i = 0; i < c.size(); i++) {
			double apost = 0;
			ChainLink cli = c.get(i);
			// Omega
			RandomVar<Double2D> rvOmega = (RandomVar<Double2D>) cli
					.get("Omega");
			Double2D Omega = rvOmega.getNumericValue();
			apost += rvOmega.getPrior().logDensity(Omega);
			// M
			RandomVar<Double2D> rvM = (RandomVar<Double2D>) cli.get("M");
			Double2D M = rvM.getNumericValue();
			apost += rvM.getPrior().logDensity(M);
			// Hyperparameters
			RandomVar<Double0D> rvAl = (RandomVar<Double0D>) cli.get("al");
			Double0D al = rvAl.getNumericValue();
			apost += rvAl.getPrior().logDensity(al);
			// Model parameters and likelihood
			RandomVar<Double3D> rvSg = (RandomVar<Double3D>) cli.get("Sg");
			RandomVar<Double3D> rvA = (RandomVar<Double3D>) cli.get("A");
			// B
			RandomVar<Double2D> rvB = (RandomVar<Double2D>) cli.get("B");
			Double2D B = rvB.getNumericValue();
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
					CRPDist postZ = (CRPDist) rvZ.getPosterior();
					apost += postZ.getBaseSg().logDensity(Sg);
					apost += new MatrixNormalDist(M, Sg, postZ.getPhi()).logDensity(A);
					// Likelihood
					Integer1D whichz = Z.which(z);
					Double2D db = B.getAll(whichz.value());
					Double2D dz = d.getAll(whichz.value());
					for (int j = 0; j < dz.size(); j++) {
						Double1D pt = dz.get(j);
						Double1D b = db.get(j);
						Double1D mean = b.mult(A);
						Double2D cov = Sg;
						apost += new MVNormalDist(mean, cov).logDensity(pt);
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

	public FLGFDModel() {
		super();
	}

	public FLGFDModel(Map<String, Numeric> hypers,
			Map<String, Numeric> init, int dims)
			throws ProbDistParmException {
		super(hypers, dims);
		// Set up parameters
		this.params = new HashMap<String, RandomVar<? extends Numeric>>();

		// Omega
		ProbDist<Double2D> priorOmega = new InverseWishartDist((Double2D) this
				.getHyper("Phi"), (Double0D) this.getHyper("lambda"));
		PostOmegaDist postOmega = new PostOmegaDist((Double2D) this
				.getHyper("Phi"), (Double0D) this.getHyper("lambda"));
		postOmega.setBaseSg(new InverseWishartDist((Double2D) this.getHyper("Psi"), (Double0D) this
				.getHyper("kappa"))); // FIXME -- Eek! Hack!
		Double2D defaultOmega = (Double2D) init.get("Omega"); // FIXME
		RandomVar<Double2D> rvOmega = new RandomVar<Double2D>("Omega",
				priorOmega, postOmega, defaultOmega);
		this.params.put("Omega", rvOmega);

		// Sg
		PriorSgDist priorSg = new PriorSgDist((Double2D) this.getHyper("Psi"),
				(Double0D) this.getHyper("kappa"));
		PostSgDist postSg = new PostSgDist((Double2D) this.getHyper("Psi"),
				(Double0D) this.getHyper("kappa"));
		Double3D defaultSg = (Double3D) init.get("Sg"); // FIXME
		RandomVar<Double3D> rvSg = new RandomVar<Double3D>("Sg", priorSg,
				postSg, defaultSg);
		this.params.put("Sg", rvSg);

		// A
		PriorADist priorA = new PriorADist();
		PostADist postA = new PostADist();
		Double3D defaultA = (Double3D) init.get("A"); // FIXME
		RandomVar<Double3D> rvA = new RandomVar<Double3D>("A", priorA, postA,
				defaultA);
		this.params.put("A", rvA);

		// Z
		ProbDist<Integer1D> priorZ = new Deterministic<Integer1D>(init.get("Z"));
		CRPDist postZ = new CRPDist((Double2D) this.getHyper("Psi"),
				(Double0D) this.getHyper("kappa"), (Double2D) this
						.getHyper("Phi"), (Double0D) this.getHyper("lambda")); // FIXME
		// --
		// Hack!
		postZ.setBaseSg(new InverseWishartDist((Double2D) this.getHyper("Psi"), (Double0D) this
				.getHyper("kappa"))); // FIXME -- Eek! Hack!
		postZ.setBaseA(new MatrixNormalDist()); // FIXME -- Eek! Hack!
		Integer1D defaultZ = (Integer1D) init.get("Z");
		RandomVar<Integer1D> rvZ = new RandomVar<Integer1D>("Z", priorZ, postZ,
				defaultZ);
		this.params.put("Z", rvZ); // FIXME

		// B
		ProbDist<Double2D> priorB = new Deterministic<Double2D>(init.get("B"));
		PostBDist postB = new PostBDist();
		Double2D defaultB = (Double2D) init.get("B"); // FIXME
		RandomVar<Double2D> rvB = new RandomVar<Double2D>("B", priorB, postB,
				defaultB);
		this.params.put("B", rvB); // FIXME

		// M
		Double2D defaultM = (Double2D) init.get("M"); // FIXME
		ProbDist<Double2D> priorM = new MatrixNormalDist((Double2D) this
				.getHyper("W"), (Double2D) this.getHyper("S"), Double2D
				.ident(defaultM.numRows()));
		PostMDist postM = new PostMDist((Double2D) this.getHyper("W"),
				(Double2D) this.getHyper("S"));
		RandomVar<Double2D> rvM = new RandomVar<Double2D>("M", priorM, postM,
				defaultM);
		this.params.put("M", rvM);

		// x
		ProbDist<Double0D> priorX = new BetaDist(
				(Double0D) this.getHyper("xa"), (Double0D) this.getHyper("xb"));
		PostXDist postX = new PostXDist();
		Double0D defaultX = (Double0D) init.get("x"); // FIXME
		RandomVar<Double0D> rvX = new RandomVar<Double0D>("x", priorX, postX,
				defaultX);
		this.params.put("x", rvX); // FIXME

		// al
		ProbDist<Double0D> priorAl = new GammaDist((Double0D) this
				.getHyper("ala"), (Double0D) this.getHyper("alb"));
		ProbDist<Double0D> postAl = new PostAlDist(
				(Double0D) this.getHyper("ala"), (Double0D) this
						.getHyper("alb"));
		Double0D defaultAl = (Double0D) init.get("al"); // FIXME
		RandomVar<Double0D> rvAl = new RandomVar<Double0D>("al", priorAl,
				postAl, defaultAl);
		this.params.put("al", rvAl); // FIXME
	}

	@Override
	public ChainLink getInitialLink() {
		return new ChainLink(this.getParam("Z"), this.getParam("B"), this
				.getParam("A"), this.getParam("Sg"), this.getParam("Omega"),
				this.getParam("M"), this.getParam("x"), this.getParam("al"));
	}
}
