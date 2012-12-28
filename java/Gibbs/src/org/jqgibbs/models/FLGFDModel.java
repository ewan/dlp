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
import org.jqgibbs.mathstat.probdist.CategoricalDist;
import org.jqgibbs.mathstat.probdist.GammaDist;
import org.jqgibbs.mathstat.probdist.InverseWishartDist;
import org.jqgibbs.mathstat.probdist.MVNormalDist;
import org.jqgibbs.mathstat.probdist.MatrixNormalDist;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistMC;

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

	class PostOmegaDist extends ProbDistMC<Double2D> {
		protected ProbDist<Double2D> baseSgBROKEN;

		private Double2D postPhi;
		private Double0D postLambda;
		private InverseWishartDist iwDist = new InverseWishartDist();

		// Fixed parameters
		private Double2D phi;
		private Double0D lambda;

		// Chain parameters
		private Integer1D Z;
		private Double2D M;
		private Double3D Sg;
		private Double3D A;

		public PostOmegaDist(Double2D phi, Double0D lambda) {
			this.phi = phi;
			this.lambda = lambda;
		}

		@Override
		public void setMCState(ChainLink l) {
			this.Z = (Integer1D) l.get("Z").getNumericValue();
			this.M = (Double2D) l.get("M").getNumericValue();
			this.Sg = (Double3D) l.get("Sg").getNumericValue();
			this.A = (Double3D) l.get("A").getNumericValue();
			this.setUpFromChainParms();
			this.initialized = true;
		}

		protected void setUpFromChainParms() {
			this.postPhi = this.phi;
			this.postLambda = this.lambda;
			Integer1D active = this.Z.items();
			int d = FLGFDModel.this.dims;
			for (Integer0D k : active) {
				Double2D AM = this.A.get(k.value()).minus(this.M);
				Double2D sgInv = this.Sg.get(k.value()).inverse();
				this.postLambda = this.postLambda.plus(d);
				this.postPhi = this.postPhi.plus(AM.mult(sgInv).mult(
						AM.transpose()));
			}
		}

		@Override
		protected Double2D genVariate() {
			boolean goodSample = false;
			Double2D o = null;
			while (!goodSample) {
				this.iwDist.setParms(this.postPhi, this.postLambda);
				o = this.iwDist.variate();
				goodSample = true;
				Integer1D active = this.Z.items();
				for (Integer0D k : active) {
					Double2D sg = this.Sg.get(k.value());
					Double2D sgKOmega = sg.kron(o);
					if (!sgKOmega.isWellConditioned()) {
						goodSample = false;
						System.err
								.println("Warning: threw out IW sample for Omega that had ill-conditioned Kronecker product with one or more current Sgs");
						break;
					}
				}
				if (goodSample) {
					Double2D sg = this.baseSgBROKEN.variate();
					Double2D sgKOmega = sg.kron(o);
					if (!sgKOmega.isWellConditioned()) {
						goodSample = false;
						System.err
								.println("Warning: threw out IW sample for Omega that had ill-conditioned Kronecker product with a sample from the Sg prior");
					}
				}
			}
			return o;
		}
	}

	class PostSgDist extends ProbDistMC<Double3D> {
		// Fixed parameters
		private Double2D Psi;
		private Double0D kappa;
		private Double2D data;

		// Chain parameters
		private Integer1D Z;
		private Double2D B;
		private Double2D Omega;
		private Double3D A;
		private Double2D M;

		private List<Double2D> postPsi;
		private List<Double0D> postKappa;

		public PostSgDist(Double2D Psi, Double0D kappa, Double2D data) {
			this.Psi = Psi;
			this.kappa = kappa;
			this.data = data;
		}

		@Override
		public void setMCState(ChainLink l) {
			this.Z = (Integer1D) l.get("Z").getNumericValue();
			this.B = (Double2D) l.get("B").getNumericValue();
			this.Omega = (Double2D) l.get("Omega").getNumericValue();
			this.A = (Double3D) l.get("A").getNumericValue();
			this.M = (Double2D) l.get("M").getNumericValue();
			this.setUpFromParms();
			this.initialized = true;
		}

		protected void setUpFromParms() {
			this.postPsi = new ArrayList<Double2D>();
			this.postKappa = new ArrayList<Double0D>();
			Integer1D active = this.Z.items();
			Double2D omegaInv = this.Omega.inverse();
			for (Integer0D k : active) {
				Integer1D zK = this.Z.which(k);
				Double2D xK = this.data.getAll(zK.value());
				Integer0D nK = new Integer0D(xK.size());

				this.postKappa.add(this.kappa.plus(nK));

				Double2D B = this.B.getAll(zK.value());
				Double2D BTB = B.transpose().mult(B);
				Double2D MTOmegaInv = this.M.transpose().mult(omegaInv);
				Double2D MTOmegaInvM = MTOmegaInv.mult(this.M);
				Double2D XTB = B.transpose().mult(xK).transpose();
				Double2D Mp = XTB.plus(MTOmegaInv);
				Double2D postOmega = omegaInv.plus(BTB).inverse();
				Double2D MpPOmegaMpT = Mp.mult(postOmega).mult(Mp.transpose());
				Double2D XTX = xK.transpose().mult(xK);
				this.postPsi.add(this.Psi.plus(XTX).plus(MTOmegaInvM)
						.minus(MpPOmegaMpT));
			}
		}

		@Override
		protected Double3D genVariate() {
			int d = FLGFDModel.this.dims;
			double[][][] v = new double[this.A.size()][d][d];
			Double3D activeUpdates = InverseWishartDist.variates(this.postPsi,
					this.postKappa);
			int i = 0;
			Integer1D active = this.Z.items();
			for (Integer0D k : active) {
				v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], d);
				i++;
			}
			return new Double3D(v);
		}
	}

	class PostADist extends ProbDistMC<Double3D> {
		private List<Double2D> postMs;
		private List<Double2D> postOmegas;

		// Fixed parameters
		private Double2D data;

		// Chain parameters
		private Double2D Omega;
		private Double3D Sg;
		private Double2D M;
		private Integer1D Z;
		private Double2D B;

		public PostADist(Double2D data) {
			this.data = data;
		}

		@Override
		public void setMCState(ChainLink l) {
			this.Omega = (Double2D) l.get("Omega").getNumericValue();
			this.Sg = (Double3D) l.get("Sg").getNumericValue();
			this.M = (Double2D) l.get("M").getNumericValue();
			this.Z = (Integer1D) l.get("Z").getNumericValue();
			this.B = (Double2D) l.get("B").getNumericValue();
			this.setUpFromParms();
			this.initialized = true;
		}

		private Double3D matnVariates() {
			Integer1D active = this.Z.items();
			Double3D sgActive = this.Sg.getAll(active.value());
			return MatrixNormalDist.variates(this.postMs, sgActive,
					this.postOmegas);
		}

		protected void setUpFromParms() {
			this.postOmegas = new ArrayList<Double2D>();
			this.postMs = new ArrayList<Double2D>();
			Integer1D active = this.Z.items();
			Double2D omegaInv = this.Omega.inverse();
			for (Integer0D k : active) {
				Integer1D zK = this.Z.which(k);
				Double2D xK = this.data.getAll(zK.value());
				Double2D B = this.B.getAll(zK.value());
				Double2D BTB = B.transpose().mult(B);
				Double2D postOmegaInv = omegaInv.plus(BTB);
				Double2D postOmega = postOmegaInv.inverse();
				this.postOmegas.add(postOmega);
				Double2D BTX = B.transpose().mult(xK);
				Double2D omegaInvM = omegaInv.mult(this.M);
				Double2D postM = postOmega.mult(BTX.plus(omegaInvM));
				this.postMs.add(postM);
			}
		}

		@Override
		protected Double3D genVariate() {
			int K = this.Sg.size();
			int h = this.B.numCols();
			int d = FLGFDModel.this.dims;
			double[][][] v = new double[K][h][d];
			boolean goodSample = false;
			while (!goodSample) {
				Double3D activeUpdates = this.matnVariates();
				goodSample = true;
				int i = 0;
				Integer1D active = this.Z.items();
				for (Integer0D k : active) {
					v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], h); // FIXME??
					for (int m = 0; m < h; m++) {
						for (int n = 0; n < d; n++) {
							if (Double.isInfinite(v[k.value()][m][n])
									|| Double.isNaN(v[k.value()][m][n])) {
								System.err
										.println("Warning: bad sample for A: "
												+ activeUpdates.get(i));
								System.err.println("Sg: "
										+ this.Sg.get(k.value()));
								System.err.println("M: " + this.M);
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
	}

	class CRPDist extends ProbDistMC<Integer1D> {
		private final Double1D ZERO_VEC = new Double1D(
				new double[FLGFDModel.this.dims]);

		// Fixed parameters
		private Double2D Psi;
		private Double0D kappa;
		private Double2D data;

		// Chain parameters
		private Double2D Omega;
		private Double3D Sg;
		private Double3D A;
		private Integer1D Z;
		private Double2D B;
		private Double2D M;
		private Double0D al;

		protected Double1D postProb;
		private CategoricalDist catDist = new CategoricalDist();
		private MVNormalDist mvnDist = new MVNormalDist();
		private InverseWishartDist sgDist = new InverseWishartDist();
		private MatrixNormalDist aDist = new MatrixNormalDist();

		private Double2D omegaInv;
		private Integer0D currMVNDist;

		public CRPDist(Double2D Psi, Double0D kappa, Double2D data) {
			this.Psi = Psi;
			this.kappa = kappa;
			this.data = data;
		}

		@Override
		public void setMCState(ChainLink l) {
			this.Omega = (Double2D) l.get("Omega").getNumericValue();
			this.Sg = (Double3D) l.get("Sg").getNumericValue();
			this.A = (Double3D) l.get("A").getNumericValue();
			this.Z = (Integer1D) l.get("Z").getNumericValue();
			this.B = (Double2D) l.get("B").getNumericValue();
			this.M = (Double2D) l.get("M").getNumericValue();
			this.al = (Double0D) l.get("al").getNumericValue();
			this.setUpFromParms();
			this.initialized = true;
		}

		protected Integer0D catVariate() {
			this.catDist.setParms(this.postProb);
			return this.catDist.variate();
		}

		protected double mvnLogDensity(Integer0D k, Double1D pt) {
			if (this.currMVNDist == null || !this.currMVNDist.equals(k)) {
				this.mvnDist.setParms(this.ZERO_VEC, this.Sg.get(k.value()));
			}
			return this.mvnDist.logDensity(pt);
		}
		
		protected void resetCurrMVNDist() {
			this.currMVNDist = null;
		}

		protected double marginalNew(Double1D pt, Double1D b) {
			Double2D omegaN = this.omegaInv.plus(b.outer(b)).inverse();
			Double2D mH = this.M.transpose().mult(this.omegaInv);
			Double2D mHX = mH.plus(pt.outer(b));
			Double2D mO = this.M.transpose().mult(this.omegaInv).mult(this.M);
			Double2D mHXO = mHX.mult(omegaN).mult(mHX.transpose());
			Double2D psiN = this.Psi.plus(pt.outer(pt).plus(mO).minus(mHXO));
			int d = this.data.numCols();
			double kN = this.kappa.value() + 1;
			double lmd = Math.log(FLGFDModel.gammaD(kN / 2, d));
			lmd -= Math.log(2) * d / 2;
			lmd -= Math.log(2 * Math.PI) * d / 2;
			lmd += Math.log(omegaN.det()) * d / 2;
			lmd -= Math.log(this.Omega.det()) * d / 2;
			lmd -= Math.log(psiN.det()) * kN / 2;
			lmd += Math.log(this.Psi.det()) * this.kappa.value() / 2;
			lmd -= Math.log(FLGFDModel.gammaD(this.kappa.value() / 2, d));
			return lmd;
		}

		protected Double2D sgVariate(Double2D psi, double k) {
			this.sgDist.setParms(psi, new Double0D(k));
			return this.sgDist.variate();
		}

		protected Double2D aVariate(Double2D M, Double2D sg, Double2D omega) {
			this.aDist.setParms(M, sg, omega);
			return this.aDist.variate();
		}

		protected List<Double2D> getVariateDPPost(Double1D pt, Double1D b) {
			Double2D omegaN = this.omegaInv.plus(b.outer(b)).inverse();
			Double2D mH = this.M.transpose().mult(this.omegaInv);
			Double2D mHX = mH.plus(pt.outer(b));
			Double2D mHXO = mHX.mult(omegaN).mult(mHX.transpose());
			Double2D mO = this.M.transpose().mult(this.omegaInv).mult(this.M);
			Double2D psiN = this.Psi.plus(pt.outer(pt).plus(mO).minus(mHXO));
			double kN = this.kappa.value() + 1;
			Double2D sgNew = null;
			Double2D aNew = null;
			boolean goodSample = false;
			while (!goodSample) {
				goodSample = true;
				sgNew = this.sgVariate(psiN, kN);
				aNew = this.aVariate(omegaN.mult(mHX.transpose()), sgNew,
						omegaN);
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

		protected void setUpFromParms() {
			this.omegaInv = this.Omega.inverse();
		}

		@Override
		protected Integer1D genVariate() {
			// FIXME - This sampling scheme destructively modifies
			// A, Sg - HACK
			int N = this.Z.size();
			Integer1D Z = null;
			try {
				Z = (Integer1D) this.Z.clone();
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace(); // FIXME
			}
			Integer0D unused = new Integer0D(-1);
			for (int i = 0; i < N; i++) {
				Double1D Xi = this.data.get(i);
				Integer1D active = Z.items();
				Z.set(i, unused);
				double[] logP = new double[active.size() + 1];
				double maxLogP = -Double.MAX_VALUE;
				int ai = 0;
				// Existing categories
				Double1D b = this.B.get(i);
				for (Integer0D k : active) {
					Integer1D zK = Z.which(k);
					int nZK = zK.size();
					if (nZK > 0) {
						Double1D aTB = b.mult(this.A.get(k.value()));
						Double1D xIA = this.data.get(i).minus(aTB);
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
				logP[newc] = this.marginalNew(Xi, b)
						+ Math.log(this.al.value());
				if (logP[ai] == Double.POSITIVE_INFINITY) {
					logP[ai] = Double.MAX_VALUE;
				}
				if (logP[newc] > maxLogP) {
					maxLogP = logP[newc];
				}
				// Select a category for this point
				Double1D p = (new Double1D(logP)).minus(maxLogP).exp();
				this.postProb = p;
				int az = this.catVariate().value();
				// Sample new category if necessary
				Integer0D zedNew = null;
				if (az == newc) {
					int maxActive = active.get(active.size() - 1).value();
					zedNew = Z.minNotIn(0, maxActive);
					this.resetCurrMVNDist();
					List<Double2D> newcParms = this.getVariateDPPost(Xi, b);
					Double2D newSg = newcParms.get(0);
					Double2D newA = newcParms.get(1);
					this.A.set(zedNew.value(), newA);
					this.Sg.set(zedNew.value(), newSg);
				} else {
					zedNew = active.get(az);
				}
				Z.set(i, zedNew);
			}
			return Z;
		}
	}

	class PostMDist extends ProbDistMC<Double2D> {
		// Fixed parameters
		protected Double1D vecW;
		protected Double2D SInv;

		protected MVNormalDist mvnDist = new MVNormalDist();
		protected Double1D postVecMu;
		protected Double2D postKronSg;

		// Chain parameters
		private Integer1D Z;
		private Double2D B;
		private Double2D Omega;
		private Double3D Sg;
		private Double3D A;

		public PostMDist(Double2D W, Double2D S) {
			this.vecW = W.colVec();
			this.SInv = S.inverse();
		}

		@Override
		public void setMCState(ChainLink l) {
			this.Z = (Integer1D) l.get("Z").getNumericValue();
			this.B = (Double2D) l.get("B").getNumericValue();
			this.Omega = (Double2D) l.get("Omega").getNumericValue();
			this.Sg = (Double3D) l.get("Sg").getNumericValue();
			this.A = (Double3D) l.get("A").getNumericValue();
			this.setUpFromParms();
			this.initialized = true;
		}

		protected void setUpFromParms() {
			Double2D sInvKIdh = this.SInv.kron(Double2D.ident(this.B.numCols()));
			Double2D postKronSgInv = sInvKIdh;
			Double1D postVecMup = sInvKIdh.mult(this.vecW);
			Integer1D active = this.Z.items();
			Double2D omegaInv = this.Omega.inverse();
			for (Integer0D k : active) {
				Double2D sgInv = this.Sg.get(k.value()).inverse();
				Double2D sgInvKOmegaInv = sgInv.kron(omegaInv);
				postKronSgInv = postKronSgInv.plus(sgInvKOmegaInv);
				Double1D vecA = this.A.get(k.value()).colVec();
				postVecMup = postVecMup.plus(sgInvKOmegaInv.mult(vecA));
			}
			this.postKronSg = postKronSgInv.inverse();
			this.postVecMu = postKronSg.mult(postVecMup);
		}

		@Override
		protected Double2D genVariate() {
			this.mvnDist.setParms(this.postVecMu, this.postKronSg);
			return this.mvnDist.variate().toDouble2D(this.B.numCols());
		}
	}

	class PostXDist extends ProbDistMC<Double0D> {
		private BetaDist betaDist = new BetaDist();
		private Double0D postXA;
		private Double0D postXB;

		// Chain parameters
		private Double0D al;
		private Integer1D Z;

		@Override
		public void setMCState(ChainLink l) {
			this.al = (Double0D) l.get("al").getNumericValue();
			this.Z = (Integer1D) l.get("Z").getNumericValue();
			this.setUpFromParms();
			this.initialized = true;
		}

		protected void setUpFromParms() {
			this.postXA = this.al.plus(1);
			this.postXB = new Double0D((double) this.Z.size());
			this.betaDist.setParms(this.postXA, this.postXB);
		}

		@Override
		protected Double0D genVariate() {
			return this.betaDist.variate();
		}

		@Override
		protected double getLogDensity(Double0D pt) {
			return this.betaDist.logDensity(pt);
		}
	}

	class PostAlDist extends ProbDistMC<Double0D> {
		private GammaDist gammaDist = new GammaDist();
		private CategoricalDist catDist = new CategoricalDist();
		private Double0D postAlA;
		private Double0D postAlB;

		// Fixed parameters
		private Double0D ala;
		private Double0D alb;

		// Chain parameters
		private Double0D x;
		private Integer1D Z;

		public PostAlDist(Double0D ala, Double0D alb) {
			this.ala = ala;
			this.alb = alb;
		}

		@Override
		public void setMCState(ChainLink l) {
			this.x = (Double0D) l.get("x").getNumericValue();
			this.Z = (Integer1D) l.get("Z").getNumericValue();
			this.setUpFromParms();
			this.initialized = true;
		}

		protected void setUpFromParms() {
			int K = this.Z.items().size();
			int N = this.Z.size();
			double logX = Math.log(this.x.value());
			double p1 = this.ala.value() + K - 1;
			double p2 = N * (this.alb.value() - logX);

			this.catDist.setParms(new Double1D(p1, p2));
			int i = this.catDist.variate().value();

			if (i == 0) {
				this.postAlA = this.ala.plus(K);
			} else {
				this.postAlA = this.ala.plus(K - 1);
			}
			this.postAlB = this.alb.plus(-logX);
			this.gammaDist.setParms(this.postAlA, this.postAlB);
		}

		@Override
		protected Double0D genVariate() {
			return this.gammaDist.variate();
		}

		@Override
		protected double getLogDensity(Double0D pt) {
			// FIXME NOW!!!
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	// FIXME
	class PostBDist extends ProbDistMC<Double2D> {
		// Chain parameters
		private Double2D B;

		public PostBDist() { }

		public void setMCState(ChainLink l) {
			this.B = (Double2D) l.get("B").getNumericValue();
			this.initialized = true;
		}
		
		@Override
		protected Double2D genVariate() {
			return this.B;
		}
	}

	@SuppressWarnings("unchecked")
	public static ChainLink pointEstimate(List<ChainLink> c, Double2D d,
			Map<String, Numeric> hypers) {
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
			InverseWishartDist iwDistOmega = new InverseWishartDist(
					(Double2D) hypers.get("Phi"),
					(Double0D) hypers.get("lambda"));
			apost += iwDistOmega.logDensity(Omega); // TODO - check
			// M
			RandomVar<Double2D> rvM = (RandomVar<Double2D>) cli.get("M");
			Double2D M = rvM.getNumericValue();
			int h = M.numCols();
			Double2D omega = Double2D.ident(h);
			MatrixNormalDist mnDist = new MatrixNormalDist(
					(Double2D) hypers.get("W"), (Double2D) hypers.get("S"),
					omega);
			apost += mnDist.logDensity(M); // TODO - check (a) dims; (b)
											// breakage
			// Hyperparameters
			RandomVar<Double0D> rvAl = (RandomVar<Double0D>) cli.get("al");
			Double0D al = rvAl.getNumericValue();
			GammaDist gammaDist = new GammaDist((Double0D) hypers.get("ala"),
					(Double0D) hypers.get("alb"));
			apost += gammaDist.logDensity(al); // TODO - check
			// Model parameters and likelihood
			RandomVar<Double3D> rvSg = (RandomVar<Double3D>) cli.get("Sg");
			RandomVar<Double3D> rvA = (RandomVar<Double3D>) cli.get("A");
			// B
			RandomVar<Double2D> rvB = (RandomVar<Double2D>) cli.get("B");
			Double2D B = rvB.getNumericValue();
			// skip
			// Z FIXME
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
			InverseWishartDist iwDistSg = new InverseWishartDist(
					(Double2D) hypers.get("Psi"),
					(Double0D) hypers.get("kappa"));
			for (Integer0D z : active) {
				try {
					// Model parameters
					Double2D Sg = ((Double3D) rvSg.getNumericValue()).get(z
							.value());
					Double2D A = ((Double3D) rvA.getNumericValue()).get(z
							.value());
					apost += iwDistSg.logDensity(Sg);
					apost += new MatrixNormalDist(M, Sg, Omega).logDensity(A); // (b)
																				// there
																				// was
																				// a
																				// bug
																				// here
																				// -
																				// should
																				// be
																				// Omega
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
					throw new RuntimeException(e);
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

	public FLGFDModel(Map<String, Numeric> hypers, Map<String, Numeric> init,
			Double2D data, int dims) {
		super(hypers, dims);
		// Set up parameters
		this.params = new HashMap<String, RandomVar<? extends Numeric>>();

		// Omega
		PostOmegaDist postOmega = new PostOmegaDist(
				(Double2D) this.getHyper("Phi"),
				(Double0D) this.getHyper("lambda"));
		postOmega.baseSgBROKEN = new InverseWishartDist(
				(Double2D) this.getHyper("Psi"),
				(Double0D) this.getHyper("kappa"));
		Double2D defaultOmega = (Double2D) init.get("Omega");
		RandomVar<Double2D> rvOmega = new RandomVar<Double2D>("Omega",
				postOmega, defaultOmega);
		this.params.put("Omega", rvOmega);

		// Sg
		PostSgDist postSg = new PostSgDist((Double2D) this.getHyper("Psi"),
				(Double0D) this.getHyper("kappa"), data);
		Double3D defaultSg = (Double3D) init.get("Sg");
		RandomVar<Double3D> rvSg = new RandomVar<Double3D>("Sg", postSg,
				defaultSg);
		this.params.put("Sg", rvSg);

		// A
		PostADist postA = new PostADist(data);
		Double3D defaultA = (Double3D) init.get("A");
		RandomVar<Double3D> rvA = new RandomVar<Double3D>("A", postA, defaultA);
		this.params.put("A", rvA);

		// Z
		CRPDist postZ = new CRPDist((Double2D) this.getHyper("Psi"),
				(Double0D) this.getHyper("kappa"), data);
		Integer1D defaultZ = (Integer1D) init.get("Z");
		RandomVar<Integer1D> rvZ = new RandomVar<Integer1D>("Z", postZ,
				defaultZ);
		this.params.put("Z", rvZ);

		// B
		PostBDist postB = new PostBDist();
		Double2D defaultB = (Double2D) init.get("B");
		RandomVar<Double2D> rvB = new RandomVar<Double2D>("B", postB, defaultB);
		this.params.put("B", rvB);

		// M
		Double2D defaultM = (Double2D) init.get("M");
		PostMDist postM = new PostMDist((Double2D) this.getHyper("W"),
				(Double2D) this.getHyper("S"));
		RandomVar<Double2D> rvM = new RandomVar<Double2D>("M", postM, defaultM);
		this.params.put("M", rvM);

		// x
		PostXDist postX = new PostXDist();
		Double0D defaultX = (Double0D) init.get("x");
		RandomVar<Double0D> rvX = new RandomVar<Double0D>("x", postX, defaultX);
		this.params.put("x", rvX);

		// al
		ProbDistMC<Double0D> postAl = new PostAlDist(
				(Double0D) this.getHyper("ala"),
				(Double0D) this.getHyper("alb"));
		Double0D defaultAl = (Double0D) init.get("al");
		RandomVar<Double0D> rvAl = new RandomVar<Double0D>("al", postAl,
				defaultAl);
		this.params.put("al", rvAl);
	}

	@Override
	public ChainLink getInitialLink() {
		return new ChainLink(this.getParam("Z"), this.getParam("B"),
				this.getParam("A"), this.getParam("Sg"),
				this.getParam("Omega"), this.getParam("M"), this.getParam("x"),
				this.getParam("al"));
	}
}
