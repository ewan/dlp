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
import org.jqgibbs.mathstat.Integer2D;
import org.jqgibbs.mathstat.ListSequence;
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
import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;
import org.jqgibbs.mathstat.probdist.SeqInverseWishartDist;
import org.jqgibbs.mathstat.probdist.SeqMatrixNormalDist;
import cern.colt.matrix.linalg.SingularValueDecomposition;

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

		public PostOmegaDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
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
			return (Double2D) this.fixedParms[0];
		}

		protected Double0D getLambda() {
			return (Double0D) this.fixedParms[1];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[0];
		}

		protected Double2D getM() {
			return (Double2D) this.chainParms[1];
		}

		protected Double3D getSg() {
			return (Double3D) this.chainParms[2];
		}

		protected Double3D getA() {
			return (Double3D) this.chainParms[3];
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

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(4);
			this.chainParmNames.add("Z");
			this.chainParmNames.add("M");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("A");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("Phi");
			this.fixedParmNames.add("lambda");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					4);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
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
				this.iwDists = new SeqInverseWishartDist(this.getRepPsi(), this
						.getRepKappa());
				return this.iwDists.variate();
			} else {
				return this.iwDists.variate(this.getRepKappa(), this
						.getRepKappa());
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

		protected Double2D getB() {
			return (Double2D) this.chainParms[1];
		}

		protected Double2D getOmega() {
			return (Double2D) this.chainParms[2];
		}

		protected Double3D getA() {
			return (Double3D) this.chainParms[3];
		}

		protected Double2D getM() {
			return (Double2D) this.chainParms[4];
		}

		private Double3D iwVariates() throws ProbDistParmException {
			if (this.iwDists == null) {
				this.iwDists = new SeqInverseWishartDist(this.getPostPsi(),
						this.getPostKappa());
				return this.iwDists.variate();
			} else {
				return this.iwDists.variate(this.getPostPsi(), this
						.getPostKappa());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(4);
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.chainParmNames.add("Omega");
			this.chainParmNames.add("A");
			this.chainParmNames.add("M");
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
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					5);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPostPsi(new ListSequence<Double2D>());
			this.setPostKappa(new ListSequence<Double0D>());
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
					FLGFDModel.this.dims);
			for (int i = 0; i < active.size(); i++) {
				repM.add(M);
			}
			return repM;
		}

		private Integer1D getZ() {
			return (Integer1D) this.chainParms[3];
		}

		private Double2D getB() {
			return (Double2D) this.chainParms[4];
		}

		private Double3D matnVariates() throws ProbDistParmException {
			if (this.matnDists == null) {
				this.matnDists = new SeqMatrixNormalDist(this.getRepM(), this
						.getActiveSg(), this.getOmega());
				return this.matnDists.variate();
			} else {
				return this.matnDists.variate(this.getRepM(), this
						.getActiveSg(), this.getOmega());
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
			this.chainParmClasses.add(Double2D.class);
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
		private SeqMatrixNormalDist matnDists;

		public PostADist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private Double3D matnVariates() throws ProbDistParmException {
			Integer1D active = this.getZ().items();
			Double3D sgActive = this.getSg().getAll(active.value());
			if (this.matnDists == null) {
				this.matnDists = new SeqMatrixNormalDist(this.getPostM(),
						sgActive, this.getPostOmegas());
				// FIXME - Will this work if all arguments are not
				// ListSequences?
				return this.matnDists.variate();
			} else {
				return this.matnDists.variate(this.getPostM(), sgActive, this
						.getPostOmegas());
			}
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

		protected Double2D getOmega() {
			return (Double2D) this.chainParms[0];
		}

		protected Double3D getSg() {
			return (Double3D) this.chainParms[1];
		}

		protected Double2D getM() {
			return (Double2D) this.chainParms[2];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[3];
		}

		protected Double2D getB() {
			return (Double2D) this.chainParms[4];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(5);
			this.chainParmNames.add("Omega");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("M");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.fixedParmNames = new ArrayList<String>(0);
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					5);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					0);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPostOmegas(new ListSequence<Double2D>());
			this.setPostM(new ListSequence<Double2D>());
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
		public CRPDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
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

		protected ListSequence<Double2D> getVariateDPPost(Double1D pt,
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
			return new ListSequence<Double2D>(sgNew, aNew);
		}

		protected Double2D getOmega() {
			return (Double2D) this.chainParms[0];
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

		protected Double2D getB() {
			return (Double2D) this.chainParms[4];
		}

		protected Double2D getM() {
			return (Double2D) this.chainParms[5];
		}

		protected Double0D getAl() {
			return (Double0D) this.chainParms[6];
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

		protected Double2D getOmegaInv() {
			return this.omegaInv;
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(7);
			this.chainParmNames.add("Omega");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("A");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.chainParmNames.add("M");
			this.chainParmNames.add("al");
			this.fixedParmNames = new ArrayList<String>(4);
			this.fixedParmNames.add("Psi");
			this.fixedParmNames.add("kappa");
			this.fixedParmNames.add("Phi");
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
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					7);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					4);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
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
				int az = this.catVariate().value();
				// Sample new category if necessary
				Integer0D zedNew = null;
				if (az == newc) {
					int maxActive = active.get(active.size() - 1).value();
					zedNew = Z.minNotIn(0, maxActive);
					this.resetPostMvn(zedNew);
					ListSequence<Double2D> newcParms = this.getVariateDPPost(Xi, b);
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

		public PostMDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[0];
		}

		protected Double2D getB() {
			return (Double2D) this.chainParms[1];
		}

		protected Double2D getOmega() {
			return (Double2D) this.chainParms[2];
		}

		protected Double3D getSg() {
			return (Double3D) this.chainParms[3];
		}

		protected Double3D getA() {
			return (Double3D) this.chainParms[4];
		}

		protected Double2D getPhi() {
			return (Double2D) this.fixedParms[2];
		}

		protected Double2D getSInv() {
			if (this.SInv == null) {
				Double2D S = (Double2D) this.fixedParms[1];
				this.SInv = S.inverse();
			}
			return this.SInv;
		}

		protected Double1D getVecW() {
			if (this.vecW == null) {
				Double2D W = (Double2D) this.fixedParms[0];
				this.vecW = W.colVec();
			}
			return this.vecW;
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(5);
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.chainParmNames.add("Omega");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("A");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("W");
			this.fixedParmNames.add("S");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
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
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double2D.class);
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

		@Override
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
		private CategoricalDistInitializeByP catDist;
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

	class PostBDist extends ProbDistInitializeByChain<Double2D> {
		public PostBDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private Double2D getB() {
			return (Double2D) this.chainParms[0];
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
		MVNormalDist mvnDist = new MVNormalDist();
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
					apost += postZ.getBaseA().logDensity(A, M, Sg,
							postZ.getPhi());
					// Likelihood
					Integer1D whichz = Z.which(z);
					Double2D db = B.getAll(whichz.value());
					Double2D dz = d.getAll(whichz.value());
					for (int j = 0; j < dz.size(); j++) {
						Double1D pt = dz.get(j);
						Double1D b = db.get(j);
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

	// // Set up matrix for pairs
	// int n = d.size();
	// int nPairs = n * n - n * (n + 1) / 2;
	// Integer2D pairsSeq = new Integer2D(c.size(), nPairs);
	// // Set up matrix for B
	// ListSequence<Integer2D> BSeq = new ListSequence<Integer2D>();
	// // For each link in the chain,
	// for (int i = 0; i < c.size(); i++) {
	// // Save pairs
	// pairsSeq.set(i, ((Integer1D) c.get(i).get("Z").getNumericValue())
	// .samePairs());
	// // Save B
	// BSeq.set(i, (Integer2D) c.get(i).get("B").getNumericValue());
	// }
	// // Take consistent pairs
	// Double1D ppair = pairsSeq.mean();
	// int[] pairsi = new int[nPairs];
	// // CategoricalDistInitializeByP dist = new
	// // CategoricalDistInitializeByP(new Double1D(0.5, 0.5));
	// for (int i = 0; i < nPairs; i++) {
	// double p = ppair.get(i).value();
	// // pairsi[i] = dist.variateFast(new Double1D((1-p), p)).value();
	// if (p > 0.75) {
	// pairsi[i] = 1;
	// }
	// // pairsi[i] = dist.variateFast(new Double1D((1-p), p)).value();
	// }
	// Integer1D pairs = new Integer1D(pairsi);
	// // Take mode B
	// Integer2D B = ((Integer2D) c.get(c.size() - 1).get("B")
	// .getNumericValue()).cloneFromVector(BSeq.modeVector());
	// // Get point estimate for Z by doing assignment to graph components
	// Integer1D Z = new UndirectedGraph(pairs).components();
	// // Initialize
	// ChainLink last = c.get(c.size() - 1);
	// int nComps = Z.max().value() + 1;
	// int nObs = ((Double3D) last.get("A").getNumericValue()).get(0)
	// .numRows();
	// int nDims = ((Double3D) last.get("A").getNumericValue()).get(0)
	// .numCols();
	// // Z
	// RandomVar<Integer1D> rvZ = (RandomVar<Integer1D>) last.get("Z")
	// .cloneFromVector(Z.rowVec());
	// // B
	// RandomVar<Integer1D> rvB = (RandomVar<Integer1D>) last.get("B")
	// .cloneFromVector(B.rowVec());
	// // M
	// Double3D Ms = new Double3D(c.size(), nObs, nDims);
	// for (int i = 0; i < c.size(); i++) {
	// Ms.set(i, (Double2D) c.get(i).get("M").getNumericValue());
	// }
	// RandomVar<Double2D> rvM = (RandomVar<Double2D>) last.get("M")
	// .cloneFromVector(Ms.meanVector());
	// Double2D M = rvM.getNumericValue();
	// // Model parameters
	// Double3D AInit = new Double3D(nComps, nObs, nDims);
	// Double3D SgInit = new Double3D(nComps, nDims, nDims); //
	// Double3D OmegaInit = new Double3D(nComps, nObs, nObs); //
	// for (Integer0D z : Z.items()) {
	// Double2D Omegai = ((CPDist) rvZ.getPosterior()).getBaseOmega()
	// .variateFast();
	// Double2D Sgi = ((CPDist) rvZ.getPosterior()).getBaseSg()
	// .variateFast();
	// Double2D Ai = ((CPDist) rvZ.getPosterior()).getBaseA().variateFast(
	// M, Sgi, Omegai);
	// SgInit.set(z.value(), (Double2D) Sgi);
	// OmegaInit.set(z.value(), (Double2D) Omegai);
	// AInit.set(z.value(), (Double2D) Ai);
	// }
	// RandomVar<Double3D> rvSg= new RandomVar<Double3D>("Sg",
	// ((RandomVar<Double3D>) last.get("Sg")).getPrior(),
	// ((RandomVar<Double3D>) last.get("Sg")).getPosterior(), SgInit);
	// RandomVar<Double3D> rvOmega = new RandomVar<Double3D>("Omega",
	// ((RandomVar<Double3D>) last.get("Omega")).getPrior(),
	// ((RandomVar<Double3D>) last.get("Omega")).getPosterior(),
	// OmegaInit);
	// RandomVar<Double3D> rvA = new RandomVar<Double3D>("A",
	// ((RandomVar<Double3D>) last.get("A")).getPrior(),
	// ((RandomVar<Double3D>) last.get("A")).getPosterior(), AInit);
	// // DP hyperparameters
	// RandomVar<Double0D> rvX = (RandomVar<Double0D>) last.get("x");
	// RandomVar<Double0D> rvAl = (RandomVar<Double0D>) last.get("al");
	// // Sample new model parameters using current prior, observing Z and B
	// Map<String, AbstractSequence<? extends AbstractSequence<?, ?>, ? extends
	// Numeric<?>>> vars = new HashMap<String, AbstractSequence<? extends
	// AbstractSequence<?, ?>, ? extends Numeric<?>>>();
	// vars.put("A", AInit.sequence());
	// vars.put("Sg", SgInit.sequence());
	// vars.put("Omega", OmegaInit.sequence());
	// ChainLink curr = new ChainLink(rvZ, rvB, rvA, rvSg, rvOmega, rvM, rvX,
	// rvAl);
	// ChainLink next = null;
	// // Now draw a big sample
	// for (int i = 0; i < 2*c.size(); i++) {
	// try {
	// next = (ChainLink) curr.clone();
	// next.get("A").updatePosteriorFast(next, d);
	// next.get("Sg").updatePosteriorFast(next, d);
	// next.get("Omega").updatePosteriorFast(next, d);
	// if (i >= c.size()) {
	// ((ListSequence<Double3D>) vars.get("A"))
	// .add((Double3D) next.get("A").getNumericValue()
	// .clone());
	// ((ListSequence<Double3D>) vars.get("Sg"))
	// .add((Double3D) next.get("Sg").getNumericValue()
	// .clone());
	// ((ListSequence<Double3D>) vars.get("Omega"))
	// .add((Double3D) next.get("Omega").getNumericValue()
	// .clone());
	// }
	// curr = next;
	// } catch (CloneNotSupportedException e) {
	// throw new RuntimeException(e);
	// } catch (Exception e) {
	// System.err.println(e.getMessage());
	// StackTraceElement[] st = e.getStackTrace();
	// for (int j = 0; j < st.length; j++) {
	// System.err.println(st[j].toString());
	// }
	// throw new ProbDistParmException(e);
	// }
	// }
	// // Construct sample mean of model parameters
	// rvSg = (RandomVar<Double3D>) next.get("Sg")
	// .cloneFromVector(vars.get("Sg").meanVector());
	// rvOmega = (RandomVar<Double3D>) next.get("Omega")
	// .cloneFromVector(vars.get("Omega").meanVector());
	// rvA = (RandomVar<Double3D>) next.get("A")
	// .cloneFromVector(vars.get("A").meanVector());
	// // rvM = (RandomVar<Double2D>)
	// next.get("M").cloneFromVector(vars.get("M").meanVector());
	// ChainLink pte = new ChainLink(rvZ, rvB, rvA, rvSg, rvOmega, rvM, rvX,
	// rvAl);
	// // Now resample Z a few times to make the spurious categories go away
	// for (int i = 0; i < 100; i++) {
	// pte.get("Z").updatePosteriorFast(pte, d);
	// }
	// return pte;
	// }

	public FLGFDModel() {
		super();
	}

	public FLGFDModel(Map<String, Numeric<? extends Numeric<?>>> hypers,
			Map<String, Numeric<? extends Numeric<?>>> init, int dims)
			throws ProbDistParmException {
		super(hypers, dims);
		// Set up parameters
		this.params = new HashMap<String, RandomVar<? extends Numeric<?>>>();

		// Omega
		ProbDist<Double2D> priorOmega = new InverseWishartDist(this
				.getHyper("Phi"), (Double0D) this.getHyper("lambda"));
		PostOmegaDist postOmega = new PostOmegaDist((Double2D) this
				.getHyper("Phi"), (Double0D) this.getHyper("lambda"));
		postOmega.setBaseSg(new InverseWishartDist(this.getHyper("Psi"), this
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
		postZ.setBaseSg(new InverseWishartDist(this.getHyper("Psi"), this
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
		ProbDistInitializeByChain<Double0D> postAl = new PostAlDist(
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

// package org.jqgibbs.models;
//
// import java.util.ArrayList;
// import java.util.Arrays;
// import java.util.HashMap;
// import java.util.List;
// import java.util.Map;
//
// import org.jqgibbs.ChainLink;
// import org.jqgibbs.Model;
// import org.jqgibbs.mathstat.AbstractSequence;
// import org.jqgibbs.mathstat.Double0D;
// import org.jqgibbs.mathstat.Double1D;
// import org.jqgibbs.mathstat.Double2D;
// import org.jqgibbs.mathstat.Double3D;
// import org.jqgibbs.mathstat.Integer0D;
// import org.jqgibbs.mathstat.Integer1D;
// import org.jqgibbs.mathstat.Integer2D;
// import org.jqgibbs.mathstat.ListSequence;
// import org.jqgibbs.mathstat.Numeric;
// import org.jqgibbs.mathstat.RandomVar;
// import org.jqgibbs.mathstat.probdist.BetaDist;
// import org.jqgibbs.mathstat.probdist.CategoricalDistInitializeByP;
// import org.jqgibbs.mathstat.probdist.GammaDist;
// import org.jqgibbs.mathstat.probdist.InverseWishartDist;
// import org.jqgibbs.mathstat.probdist.MVNormalDist;
// import org.jqgibbs.mathstat.probdist.MatrixNormalDist;
// import org.jqgibbs.mathstat.probdist.ProbDist;
// import org.jqgibbs.mathstat.probdist.ProbDistInitializeByChain;
// import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
// import org.jqgibbs.mathstat.probdist.ProbDistParmException;
// import org.jqgibbs.mathstat.probdist.SeqInverseWishartDist;
// import org.jqgibbs.mathstat.probdist.SeqMVNormalDist;
// import org.jqgibbs.mathstat.probdist.SeqMatrixNormalDist;
//
// public class FLGFDModel extends Model {
//
// class PriorOmegaDist extends ProbDistInitializeByChain<Double2D> {
// private InverseWishartDist iwDist; // FIXME - Shouldn't this use the IID?
//
//
// public PriorOmegaDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// private Double2D getPhi() {
// return (Double2D) this.fixedParms[0];
// }
//
// private Double0D getLambda() {
// return (Double0D) this.fixedParms[1];
// }
//
// private Integer1D getZ() {
// return (Integer1D) this.chainParms[0];
// }
//
// private Double2D iwVariate() throws ProbDistParmException {
// if (this.iwDist == null) {
// this.iwDist = new InverseWishartDist(this.getPhi(), this.getLambda());
// return this.iwDist.variate();
// } else {
// return this.iwDist.variate(this.getPhi(), this.getLambda());
// }
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(1);
// this.chainParmNames.add("Z");
// this.fixedParmNames = new ArrayList<String>(2);
// this.fixedParmNames.add("Phi");
// this.fixedParmNames.add("lambda");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
// this.chainParmCheck.add(null); // FIXME
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 1);
// this.chainParmClasses.add(Integer1D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.fixedParmClasses.add(Double2D.class);
// this.fixedParmClasses.add(Double0D.class);
// }
//
// @Override
// protected void setUpFromChainParms() {
// return;
// }
//
// @Override
// protected Double2D genVariate() throws ProbDistParmException {
// return this.iwVariate();
// }
//
// @Override
// protected double getDensity(Double2D pt) {
// throw new UnsupportedOperationException(
// "Too lazy, come back later");
// }
// }
//
// class PostOmegaDist extends ProbDistInitializeByChain<Double2D> {
// private Double2D postPhi;
// private Double0D postLambda;
// private InverseWishartDist iwDist;
//
// public PostOmegaDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//		
// protected Double2D getPostPhi() {
// return this.postPhi;
// }
//
// protected Double0D getPostLambda() {
// return this.postLambda;
// }
//
// protected void setPostPhi(Double2D postPhi) {
// this.postPhi = postPhi;
// }
//
// protected void setPostLambda(Double0D postLambda) {
// this.postLambda = postLambda;
// }
//
// protected Double2D getPhi() {
// return (Double2D) this.fixedParms[0];
// }
//
// protected Double0D getLambda() {
// return (Double0D) this.fixedParms[1];
// }
//
// protected Integer1D getZ() {
// return (Integer1D) this.chainParms[0];
// }
//
// protected Double2D getM() {
// return (Double2D) this.chainParms[1];
// }
//		
// protected Double3D getSg() {
// return (Double3D) this.chainParms[2];
// }
//		
// protected Double3D getA() {
// return (Double3D) this.chainParms[3];
// }
//
// private Double2D iwVariate() throws ProbDistParmException {
// if (this.iwDist == null) {
// this.iwDist = new InverseWishartDist(this.getPostPhi(),
// this.getPostLambda());
// return this.iwDist.variate();
// } else {
// return this.iwDist.variate(this.getPostPhi(), this
// .getPostLambda());
// }
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(4);
// this.chainParmNames.add("Z");
// this.chainParmNames.add("M");
// this.chainParmNames.add("Sg");
// this.chainParmNames.add("A");
// this.fixedParmNames = new ArrayList<String>(2);
// this.fixedParmNames.add("Phi");
// this.fixedParmNames.add("lambda");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 4);
// this.chainParmClasses.add(Integer1D.class);
// this.chainParmClasses.add(Double2D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double3D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.fixedParmClasses.add(Double2D.class);
// this.fixedParmClasses.add(Double0D.class);
// }
//
// @Override
// protected void setUpFromChainParms() {
// this.setPostPhi(this.getPhi());
// Integer1D active = this.getZ().items();
// int d = FLGFDModel.this.dims;
// this.setPostLambda(this.getLambda().plus(d*active.size()));
// for (Integer0D k : active) {
// Double2D AM = this.getA().get(k.value()).minus(this.getM());
// Double2D sgInv = this.getSg().get(k.value()).inverse();
// this.setPostPhi(this.getPostPhi().plus(AM.mult(sgInv).mult(AM.transpose())));
// }
// }
//
// @Override
// protected Double2D genVariate() throws ProbDistParmException {
// return this.iwVariate();
// }
//
// @Override
// protected double getDensity(Double2D pt)
// throws ProbDistParmException {
// throw new UnsupportedOperationException(
// "Too lazy, come back later");
// }
// }
//
// class PriorSgDist extends ProbDistInitializeByChain<Double3D> {
// private SeqInverseWishartDist iwDists; // FIXME - Shouldn't this use the IID?
// private ListSequence<Double2D> repPsi;
// private ListSequence<Double0D> repKappa;
//
// public PriorSgDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//		
// private ListSequence<Double2D> getRepPsi() {
// if (this.repPsi == null) {
// this.repPsi = new ListSequence<Double2D>();
// Integer1D active = this.getZ().items();
// for (int i = 0; i<active.size(); i++) {
// this.repPsi.add(this.getPsi());
// }
// }
// return this.repPsi;
// }
//
// private ListSequence<Double0D> getRepKappa() {
// if (this.repKappa == null) {
// this.repKappa = new ListSequence<Double0D>();
// Integer1D active = this.getZ().items();
// for (int i = 0; i < active.size(); i++) {
// this.repKappa.add(this.getKappa());
// }
// }
// return this.repKappa;
// }
//
// private Double2D getPsi() {
// return (Double2D) this.fixedParms[0];
// }
//
// private Double0D getKappa() {
// return (Double0D) this.fixedParms[1];
// }
//
// private Integer1D getZ() {
// return (Integer1D) this.chainParms[0];
// }
//
// private Double3D iwVariates() throws ProbDistParmException {
// if (this.iwDists == null) {
// this.iwDists = new SeqInverseWishartDist(this.getRepPsi(),
// this.getRepKappa());
// return this.iwDists.variate();
// } else {
// return this.iwDists.variate(this.getRepKappa(), this.getRepKappa());
// }
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(1);
// this.chainParmNames.add("Z");
// this.fixedParmNames = new ArrayList<String>(2);
// this.fixedParmNames.add("Psi");
// this.fixedParmNames.add("kappa");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
// this.chainParmCheck.add(null); // FIXME
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 1);
// this.chainParmClasses.add(Integer1D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.fixedParmClasses.add(Double2D.class);
// this.fixedParmClasses.add(Double0D.class);
// }
//
// @Override
// protected void setUpFromChainParms() {
// return;
// }
//
// @Override
// protected Double3D genVariate() throws ProbDistParmException {
// return this.iwVariates();
// }
//
// @Override
// protected double getDensity(Double3D pt) {
// throw new UnsupportedOperationException(
// "Too lazy, come back later");
// }
// }
//	
// class PostSgDist extends ProbDistInitializeByChain<Double3D> {
// private ListSequence<Double2D> postPsi;
// private ListSequence<Double0D> postKappa;
// private SeqInverseWishartDist iwDists;
//
// public PostSgDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//		
// protected ListSequence<Double2D> getPostPsi() {
// return this.postPsi;
// }
//
// protected ListSequence<Double0D> getPostKappa() {
// return this.postKappa;
// }
//
// protected void setPostPsi(ListSequence<Double2D> postPsi) {
// this.postPsi = postPsi;
// }
//
// protected void setPostKappa(ListSequence<Double0D> postKappa) {
// this.postKappa = postKappa;
// }
//
// protected Double2D getPsi() {
// return (Double2D) this.fixedParms[0];
// }
//
// protected Double0D getKappa() {
// return (Double0D) this.fixedParms[1];
// }
//
// protected Integer1D getZ() {
// return (Integer1D) this.chainParms[0];
// }
//
// protected Integer2D getB() {
// return (Integer2D) this.chainParms[1];
// }
//
// protected Double2D getOmega() {
// return (Double2D) this.chainParms[2];
// }
//
// protected Double3D getA() {
// return (Double3D) this.chainParms[3];
// }
//
// protected Double2D getM() {
// return (Double2D) this.chainParms[4];
// }
//		
// private Double3D iwVariates() throws ProbDistParmException {
// if (this.iwDists == null) {
// this.iwDists = new SeqInverseWishartDist(this.getPostPsi(),
// this.getPostKappa());
// return this.iwDists.variate();
// } else {
// return this.iwDists.variate(this.getPostPsi(), this
// .getPostKappa());
// }
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(5);
// this.chainParmNames.add("Z");
// this.chainParmNames.add("B");
// this.chainParmNames.add("Omega");
// this.chainParmNames.add("A");
// this.chainParmNames.add("M");
// this.fixedParmNames = new ArrayList<String>(2);
// this.fixedParmNames.add("Psi");
// this.fixedParmNames.add("kappa");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 5);
// this.chainParmClasses.add(Integer1D.class);
// this.chainParmClasses.add(Integer2D.class);
// this.chainParmClasses.add(Double2D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double2D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.fixedParmClasses.add(Double2D.class);
// this.fixedParmClasses.add(Double0D.class);
// }
//
// @Override
// protected void setUpFromChainParms() {
// this.setPostPsi(new ListSequence<Double2D>());
// this.setPostKappa(new ListSequence<Double0D>());
// Integer1D active = this.getZ().items();
// int h = this.getA().size();
// for (Integer0D k : active) {
// Integer1D zK = this.getZ().which(k);
// Double2D xK = this.getSamplerData().getAll(zK.value());
// Integer0D nK = new Integer0D(xK.size());
//				
// this.getPostKappa().add(this.getKappa().plus(nK));
//				
// Integer2D B = this.getB().getAll(zK.value());
// Double2D BTB = B.transpose().mult(B);
// Double2D omegaInv = this.getOmega().inverse();
// Double2D MTOmegaInv = this.getM().transpose().mult(omegaInv);
// Double2D MTOmegaInvM = MTOmegaInv.mult(this.getM());
// Double2D XTB = B.transpose().mult(xK).transpose();
// Double2D Mp = XTB.plus(MTOmegaInv);
// Double2D postOmega = omegaInv.plus(BTB).inverse();
// Double2D MpPOmegaMpT = Mp.mult(postOmega).mult(Mp.transpose());
// Double2D XTX = xK.transpose().mult(xK);
// this.getPostPsi().add(this.getPsi().plus(XTX).plus(MTOmegaInvM).minus(MpPOmegaMpT));
// }
// }
//
// @Override
// protected Double3D genVariate() throws ProbDistParmException {
// int d = FLGFDModel.this.dims;
// double[][][] v = new double[this.getA().size()][d][d];
// Double3D activeUpdates = this.iwVariates();
// int i=0;
// Integer1D active = this.getZ().items();
// for (Integer0D k : active) {
// v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], d);
// i++;
// }
// return new Double3D(v);
// }
//
// @Override
// protected double getDensity(Double3D pt)
// throws ProbDistParmException {
// throw new UnsupportedOperationException(
// "Too lazy, come back later");
// }
// }
//	
// class PriorADist extends ProbDistInitializeByChain<Double3D> {
// private SeqMatrixNormalDist matnDists;
//		
// public PriorADist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//		
// private Double3D getRepOmega() {
// Double2D Omega = (Double2D) this.chainParms[0];
// Integer1D active = this.getZ().items();
// Double3D repOmega = new
// Double3D(active.size(),this.getB().numCols(),this.getB().numCols());
// for (int i = 0; i < active.size(); i++) {
// repOmega.add(Omega);
// }
// return repOmega;
// }
//
// private Double3D getActiveSg() {
// Double3D sg = (Double3D) this.chainParms[1];
// Integer1D active = this.getZ().items();
// return sg.getAll(active.value());
// }
//
// private Double3D getRepM() {
// Double2D M = (Double2D) this.chainParms[2];
// Integer1D active = this.getZ().items();
// Double3D repM = new
// Double3D(active.size(),this.getB().numCols(),FLGFDModel.this.dims);
// for (int i = 0; i < active.size(); i++) {
// repM.add(M);
// }
// return repM;
// }
//		
// private Integer1D getZ() {
// return (Integer1D) this.chainParms[3];
// }
//		
// private Integer2D getB() {
// return (Integer2D) this.chainParms[4];
// }
//
// private Double3D matnVariates() throws ProbDistParmException {
// if (this.matnDists == null) {
// this.matnDists = new SeqMatrixNormalDist(this.getRepM(), this.getActiveSg(),
// this.getRepOmega());
// return this.matnDists.variate();
// } else {
// return this.matnDists.variate(this.getRepM(), this.getActiveSg(),
// this.getRepOmega());
// }
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(5);
// this.chainParmNames.add("Omega");
// this.chainParmNames.add("Sg");
// this.chainParmNames.add("M");
// this.chainParmNames.add("Z");
// this.chainParmNames.add("B");
// this.fixedParmNames = new ArrayList<String>(0);
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
// this.chainParmCheck.add(null); // FIXME
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(5);
// this.chainParmClasses.add(Double2D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Integer1D.class);
// this.chainParmClasses.add(Integer2D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 0);
// }
//
// @Override
// protected void setUpFromChainParms() {
// return;
// }
//
// @Override
// protected Double3D genVariate() throws ProbDistParmException {
// return this.matnVariates();
// }
//
// @Override
// protected double getDensity(Double3D pt) {
// throw new UnsupportedOperationException(
// "Too lazy, come back later");
// }
// }
//
// class PostADist extends ProbDistInitializeByChain<Double3D> {
// private ListSequence<Double2D> postMs;
// private ListSequence<Double2D> postOmegas;
// private SeqMatrixNormalDist matnDists;
//		
// public PostADist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// private Double3D matnVariates() throws ProbDistParmException {
// Integer1D active = this.getZ().items();
// Double3D sgActive = this.getSg().getAll(active.value());
// if (this.matnDists == null) {
// this.matnDists = new SeqMatrixNormalDist(this.getPostM(), sgActive,
// this.getPostOmega());
// // FIXME - Will this work if all arguments are not ListSequences?
// return this.matnDists.variate();
// } else {
// return this.matnDists.variate(this.getPostM(), sgActive,
// this.getPostOmega());
// }
// }
//
// protected ListSequence<Double2D> getPostM() {
// return this.postMs;
// }
//		
// protected ListSequence<Double2D> getPostOmega() {
// return this.postOmegas;
// }
//
// protected void setPostM(ListSequence<Double2D> postMs) {
// this.postMs = postMs;
// }
//
// protected void setPostOmega(ListSequence<Double2D> postOmegas) {
// this.postOmegas = postOmegas;
// }
//
// protected Double2D getOmega() {
// return (Double2D) this.chainParms[0];
// }
//		
// protected Double3D getSg() {
// return (Double3D) this.chainParms[1];
// }
//		
// protected Double2D getM() {
// return (Double2D) this.chainParms[2];
// }
//
// protected Integer1D getZ() {
// return (Integer1D) this.chainParms[3];
// }
//
// protected Integer2D getB() {
// return (Integer2D) this.chainParms[4];
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(5);
// this.chainParmNames.add("Omega");
// this.chainParmNames.add("Sg");
// this.chainParmNames.add("M");
// this.chainParmNames.add("Z");
// this.chainParmNames.add("B");
// this.fixedParmNames = new ArrayList<String>(0);
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
// this.chainParmCheck.add(null); // FIXME
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(5);
// this.chainParmClasses.add(Double2D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double2D.class);
// this.chainParmClasses.add(Integer1D.class);
// this.chainParmClasses.add(Integer2D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(0);
// }
//
// @Override
// protected void setUpFromChainParms() {
// this.setPostOmega(new ListSequence<Double2D>());
// this.setPostM(new ListSequence<Double2D>());
// Integer1D active = this.getZ().items();
// for (Integer0D k : active) {
// Integer1D zK = this.getZ().which(k);
// Double2D xK = this.getSamplerData().getAll(zK.value());
// Integer2D B = this.getB().getAll(zK.value());
// Double2D BTB = B.transpose().mult(B);
// Double2D omegaInv = this.getOmega().inverse();
// Double2D postOmegaInv = omegaInv.plus(BTB);
// Double2D postOmega = postOmegaInv.inverse();
// this.getPostOmega().add(postOmega);
// Double2D BTX = B.transpose().mult(xK);
// Double2D omegaInvM = omegaInv.mult(this.getM());
// Double2D postM = postOmega.mult(BTX.plus(omegaInvM));
// this.getPostM().add(postM);
// }
// }
//
// @Override
// protected Double3D genVariate() throws ProbDistParmException {
// int K = this.getSg().size();
// int h = this.getB().numCols();
// int d = FLGFDModel.this.dims;
// double[][][] v = new double[K][h][d];
// Double3D activeUpdates = this.matnVariates();
// int i=0;
// Integer1D active = this.getZ().items();
// for (Integer0D k : active) {
// v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], h); // FIXME?? Check
// dims
// i++;
// }
// return new Double3D(v);
// }
//
// @Override
// protected double getDensity(Double3D pt) {
// throw new UnsupportedOperationException(
// "Too lazy, come back later");
// }
// }
//	
// class CPDist extends ProbDistInitializeByChain<Integer1D> {
// public CPDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// protected Double1D postProb;
// private CategoricalDistInitializeByP catDist;
// private HashMap<Integer0D,MVNormalDist> mvnDists; // Shouldn't you be using
// (and hacking) SeqProbDist for this? FIXME
//
// private ProbDist<Double2D> baseSg;
// private ProbDist<Double2D> baseA;
//
// public ProbDist<Double2D> getBaseSg() {
// return this.baseSg;
// }
//		
// public void setBaseSg(ProbDist<Double2D> baseSg) {
// this.baseSg = baseSg;
// }
//
// public ProbDist<Double2D> getBaseA() {
// return this.baseA;
// }
//		
// public void setBaseA(ProbDist<Double2D> baseA) {
// this.baseA = baseA;
// }
//
// protected Integer0D catVariate() throws ProbDistParmException {
// if (this.catDist == null) {
// this.catDist = new CategoricalDistInitializeByP(this
// .getPostProb());
// return this.catDist.variate();
// } else {
// return this.catDist.variate(this.getPostProb());
// }
// }
//
// protected double mvnLogDensity(Integer0D k, Double1D pt)
// throws ProbDistParmException {
// return this.mvnDists.get(k).logDensity(pt);
// }
//
// protected void setUpPostMvn(Integer0D k) throws ProbDistParmException {
// if (this.mvnDists == null) {
// this.mvnDists = new HashMap<Integer0D,MVNormalDist>();
// }
// if (!this.mvnDists.containsKey(k)) {
// double[] zeroD = new double[FLGFDModel.this.dims];
// Double1D zero = new Double1D(zeroD);
// this.mvnDists.put(k, new MVNormalDist(zero, this.getSg().get(k.value())));
// }
// // else do nothing
// }
//		
// protected void resetPostMvn(Integer0D k) {
// if (this.mvnDists != null && this.mvnDists.containsKey(k)) {
// this.mvnDists.remove(k);
// }
// }
//
// private Double1D getPostProb() {
// return this.postProb;
// }
//
// protected void setPostProb(Double1D prob) {
// this.postProb = prob;
// }
//		
// protected Double2D getOmega() {
// return (Double2D) this.chainParms[0];
// }
//		
// protected Double3D getSg() {
// return (Double3D) this.chainParms[1];
// }
//
// protected Double3D getA() {
// return (Double3D) this.chainParms[2];
// }
//
// protected Integer1D getZ() {
// return (Integer1D) this.chainParms[3];
// }
//
// protected Integer2D getB() {
// return (Integer2D) this.chainParms[4];
// }
//
// protected Double2D getM() {
// return (Double2D) this.chainParms[5];
// }
//
// protected Double0D getAl() {
// return (Double0D) this.chainParms[6];
// }
//
// protected Double2D getPsi() {
// return (Double2D) this.fixedParms[0];
// }
//
// protected Double0D getKappa() {
// return (Double0D) this.fixedParms[1];
// }
//		
// protected Double2D getPhi() {
// return (Double2D) this.fixedParms[2];
// }
//
// protected Double0D getLambda() {
// return (Double0D) this.fixedParms[3];
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(7);
// this.chainParmNames.add("Omega");
// this.chainParmNames.add("Sg");
// this.chainParmNames.add("A");
// this.chainParmNames.add("Z");
// this.chainParmNames.add("B");
// this.chainParmNames.add("M");
// this.chainParmNames.add("al");
// this.fixedParmNames = new ArrayList<String>(4);
// this.fixedParmNames.add("Psi");
// this.fixedParmNames.add("kappa");
// this.fixedParmNames.add("Phi");
// this.fixedParmNames.add("lambda");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(7);
// this.chainParmCheck.add(null); // FIXME
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(7);
// this.chainParmClasses.add(Double2D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Integer1D.class);
// this.chainParmClasses.add(Integer2D.class);
// this.chainParmClasses.add(Double2D.class);
// this.chainParmClasses.add(Double0D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(4);
// this.fixedParmClasses.add(Double2D.class);
// this.fixedParmClasses.add(Double0D.class);
// this.fixedParmClasses.add(Double2D.class);
// this.fixedParmClasses.add(Double0D.class);
// }
//
// @Override
// protected void setUpFromChainParms() {
// this.mvnDists = null;
// }
//
// @Override
// protected Integer1D genVariate() throws ProbDistParmException {
// // FIXME - This sampling scheme destructively modifies
// // A, Sg, Omega - HACK
// int N = this.getSamplerData().size();
// Integer1D Z = null;
// try {
// Z = (Integer1D) this.getZ().clone();
// } catch (CloneNotSupportedException e) {
// // TODO Auto-generated catch block
// e.printStackTrace(); // FIXME
// }
// Integer0D unused = new Integer0D(-1);
// for (int i = 0; i < N; i++) {
// // Sample a new category if necessary, remove the current
// // category from consideration
// Integer1D active = Z.items();
// int maxActive = active.get(active.size() - 1).value();
// Integer0D zedI = Z.get(i);
// Z.set(i, unused);
// Integer1D zZedI = Z.which(zedI);
// int nZZedI = zZedI.size();
// if (nZZedI > 0) { // Sample new category
// Integer0D zedNew = Z.minNotIn(0, maxActive);
// this.resetPostMvn(zedNew);
// Double2D sgNew = this.getBaseSg().variate(this.getPsi(), this.getKappa());
// Double2D aNew = this.getBaseA().variate(this.getM(), sgNew, this.getOmega());
// // Add it
// this.getSg().set(zedNew.value(), sgNew);
// this.getA().set(zedNew.value(), aNew);
// active.add(zedNew);
// }
// // Select a category for this point
// double[] logP = new double[active.size()]; // FIXME - Probably you could
// speed this up by a lot
// double maxLogP = -Double.MAX_VALUE; // FIXME??
// int ai = 0;
// for (Integer0D k : active) {
// Integer1D zK = Z.which(k);
// Integer1D b = this.getB().get(i);
// Double1D aTB = b.mult(this.getA().get(k.value()));
// Double1D xIA = this.getSamplerData().get(i).minus(aTB);
// int nZK = zK.size();
// this.setUpPostMvn(k);
// double logPrior;
// if (nZK == 0) {
// logPrior = Math.log(this.getAl().value());
// } else {
// logPrior = Math.log(nZK);
// }
// logP[ai] = this.mvnLogDensity(k, xIA) + logPrior;
// if (logP[ai] > maxLogP) {
// maxLogP = logP[ai];
// }
// ai++;
// }
// Double1D p = (new Double1D(logP)).minus(maxLogP).exp();
// this.setPostProb(p);
// int az = this.catVariate().value();
// Z.set(i, active.get(az));
// }
// return Z;
// }
//		
// @Override
// protected double getDensity(Integer1D pt) {
// throw new UnsupportedOperationException(
// "Too lazy, come back later");
// }
// }
//
// class PostMDist extends ProbDistInitializeByChain<Double2D> {
// protected MVNormalDist mvnDist;
// protected Double1D vecW;
// protected Double2D SInv;
// protected Double1D postVecMu;
// protected Double2D postKronSg;
//		
// public PostMDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// protected Integer1D getZ() {
// return (Integer1D) this.chainParms[0];
// }
//
// protected Integer2D getB() {
// return (Integer2D) this.chainParms[1];
// }
//
// protected Double2D getOmega() {
// return (Double2D) this.chainParms[2];
// }
//		
// protected Double3D getSg() {
// return (Double3D) this.chainParms[3];
// }
//
// protected Double3D getA() {
// return (Double3D) this.chainParms[4];
// }
//		
// protected Double2D getSInv() {
// if (this.SInv == null) {
// Double2D S = (Double2D) this.fixedParms[1];
// this.SInv = S.inverse();
// }
// return this.SInv;
// }
//		
// protected Double1D getVecW() {
// if (this.vecW == null) {
// Double2D W = (Double2D) this.fixedParms[0];
// this.vecW = W.colVec();
// }
// return this.vecW;
// }
//		
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(5);
// this.chainParmNames.add("Z");
// this.chainParmNames.add("B");
// this.chainParmNames.add("Omega");
// this.chainParmNames.add("Sg");
// this.chainParmNames.add("A");
// this.fixedParmNames = new ArrayList<String>(2);
// this.fixedParmNames.add("W");
// this.fixedParmNames.add("S");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 5);
// this.chainParmClasses.add(Integer1D.class);
// this.chainParmClasses.add(Integer2D.class);
// this.chainParmClasses.add(Double2D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double3D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.fixedParmClasses.add(Double2D.class);
// this.fixedParmClasses.add(Double2D.class);
// }
//
// protected Double2D matnVariate() throws ProbDistParmException {
// int h = this.getA().numRows();
// Double1D vecVariate;
// if (this.mvnDist == null) {
// this.mvnDist = new MVNormalDist(this.getPostVecMu(), this
// .getPostKronSg());
// vecVariate = this.mvnDist.variate();
//				
// } else {
// vecVariate = this.mvnDist.variate(this.getPostVecMu(), this
// .getPostKronSg());
// }
// return vecVariate.toDouble2D(h);
// }
//
// @Override
// protected void setUpFromChainParms() {
// int h = this.getB().numCols();
// Double2D sInvKIdh = this.getSInv().kron(Double2D.ident(h));
// Double2D postKronSgInv = sInvKIdh;
// Double1D postVecMup = sInvKIdh.mult(this.getVecW());
// Integer1D active = this.getZ().items();
// for (Integer0D k : active) {
// Double2D sgInv = this.getSg().get(k.value()).inverse();
// Double2D omegaInv = this.getOmega().inverse();
// Double2D sgInvKOmegaInv = sgInv.kron(omegaInv);
// postKronSgInv = postKronSgInv.plus(sgInvKOmegaInv);
// Double1D vecA = this.getA().get(k.value()).colVec();
// postVecMup = postVecMup.plus(sgInvKOmegaInv.mult(vecA));
// }
// Double2D postKronSg = postKronSgInv.inverse();
// this.setPostKronSg(postKronSg);
// this.setPostVecMu(postKronSg.mult(postVecMup));
// }
//
// protected Double1D getPostVecMu() {
// return this.postVecMu;
// }
//
// protected void setPostVecMu(Double1D mu) {
// this.postVecMu = mu;
// }
//
// protected Double2D getPostKronSg() {
// return this.postKronSg;
// }
//
// protected void setPostKronSg(Double2D sg) {
// this.postKronSg = sg;
// }
//
// @Override
// protected Double2D genVariate() throws ProbDistParmException {
// return this.matnVariate();
// }
//
// @Override
// protected double getDensity(Double2D pt) {
// throw new UnsupportedOperationException(
// "Too lazy, come back later");
// }
// }
//	
// class PostXDist extends ProbDistInitializeByChain<Double0D> {
// private BetaDist betaDist;
// private Double0D postXA;
// private Double0D postXB;
//	
// public PostXDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//		
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(1);
// this.chainParmNames.add("al");
// this.fixedParmNames = new ArrayList<String>(0);
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 1);
// this.chainParmClasses.add(Double0D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 0);
// }
//	
// private Double0D betaVariate() throws ProbDistParmException {
// if (this.betaDist == null) {
// this.betaDist = new BetaDist(this.getPostXA(), this
// .getPostXB());
// return this.betaDist.variate();
// } else {
// return this.betaDist.variate(this.getPostXA(), this
// .getPostXB());
// }
// }
//	
// protected Double0D getAl() {
// return (Double0D) this.chainParms[0];
// }
//	
// @Override
// protected void setUpFromChainParms() {
// this.setPostXA(this.getAl().plus(1));
// this.setPostXB(new Double0D((double) this.getSamplerData()
// .size()));
// }
//	
// private Double0D getPostXA() {
// return this.postXA;
// }
//	
// private void setPostXA(Double0D xA) {
// this.postXA = xA;
// }
//	
// private Double0D getPostXB() {
// return this.postXB;
// }
//	
// private void setPostXB(Double0D xB) {
// this.postXB = xB;
// }
//	
// @Override
// protected Double0D genVariate() throws ProbDistParmException {
// return this.betaVariate();
// }
//	
// @Override
// protected double getDensity(Double0D pt) {
// throw new UnsupportedOperationException(
// "Too lazy, come back later");
// }
// }
//	
// class PostAlDist extends ProbDistInitializeByChain<Double0D> {
// private GammaDist gammaDist;
// private CategoricalDistInitializeByP catDist;
// private Double0D postAlA;
// private Double0D postAlB;
// private Double1D postMix;
//
// public PostAlDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//		
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(2);
// this.chainParmNames.add("x");
// this.chainParmNames.add("Z");
// this.fixedParmNames = new ArrayList<String>(2);
// this.fixedParmNames.add("ala");
// this.fixedParmNames.add("alb");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.chainParmClasses.add(Double0D.class);
// this.chainParmClasses.add(Integer1D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.fixedParmClasses.add(Double0D.class);
// this.fixedParmClasses.add(Double0D.class);
// }
//		
// private Double0D gammaVariate() throws ProbDistParmException {
// if (this.gammaDist == null) {
// this.gammaDist = new GammaDist(this.getPostAlA(), this
// .getPostAlB());
// return this.gammaDist.variate();
// } else {
// return this.gammaDist.variate(this.getPostAlA(), this
// .getPostAlB());
// }
// }
//		
// private Integer0D catVariate() throws ProbDistParmException {
// if (this.catDist == null) {
// this.catDist = new CategoricalDistInitializeByP(this
// .getPostMix());
// return this.catDist.variate();
// } else {
// return this.catDist.variate(this.getPostMix());
// }
// }
//		
// protected Double0D getX() {
// return (Double0D) this.chainParms[0];
// }
//		
// protected Integer1D getZ() {
// return (Integer1D) this.chainParms[1];
// }
//		
// protected Double0D getAlA() {
// return (Double0D) this.fixedParms[0];
// }
//		
// protected Double0D getAlB() {
// return (Double0D) this.fixedParms[1];
// }
//		
// @Override
// protected void setUpFromChainParms() {
// int K = this.getZ().items().size();
// int N = this.getSamplerData().size();
// double logX = Math.log(this.getX().value());
// double p1 = this.getAlA().value() + K - 1;
// double p2 = N * (this.getAlB().value() - logX);
// this.setPostMix(new Double1D(p1, p2));
// int i;
// try {
// i = this.catVariate().value();
// } catch (ProbDistParmException e) {
// // FIXME
// throw new RuntimeException(e);
// }
// if (i == 0) {
// this.setPostAlA(this.getAlA().plus(K));
// } else {
// this.setPostAlA(this.getAlA().plus(K - 1));
// }
// this.setPostAlB(this.getAlB().plus(-logX)); // FIXME
// }
//		
// private Double0D getPostAlA() {
// return this.postAlA;
// }
//		
// private void setPostAlA(Double0D alA) {
// this.postAlA = alA;
// }
//		
// private Double0D getPostAlB() {
// return this.postAlB;
// }
//		
// private void setPostAlB(Double0D alB) {
// this.postAlB = alB;
// }
//		
// private Double1D getPostMix() {
// return this.postMix;
// }
//		
// private void setPostMix(Double1D mix) {
// this.postMix = mix;
// }
//		
// @Override
// protected Double0D genVariate() throws ProbDistParmException {
// // FIXME - create
// return this.gammaVariate();
// }
//		
// @Override
// protected double getDensity(Double0D pt) {
// throw new UnsupportedOperationException(
// "Too lazy, come back later");
// }
// }
//	
// class PostBDist extends ProbDistInitializeByChain<Integer2D> {
// public PostBDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//		
// private Integer2D getB() {
// return (Integer2D) this.chainParms[0];
// }
//		
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(1);
// this.chainParmNames.add("B");
// this.fixedParmNames = new ArrayList<String>(0);
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
// this.chainParmCheck.add(null); // FIXME
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(1);
// this.chainParmClasses.add(Integer2D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(0);
// }
//
// @Override
// protected void setUpFromChainParms() {
// return;
// }
//
// @Override
// protected Integer2D genVariate() throws ProbDistParmException {
// return this.getB();
// }
//
// @Override
// protected double getDensity(Integer2D pt) {
// throw new UnsupportedOperationException(
// "Too lazy, come back later");
// }
// }
//	
// @SuppressWarnings("unchecked")
// public static ChainLink pointEstimate(List<ChainLink> c, Double2D d)
// throws ProbDistParmException {
// // MAP
// double max = Double.NEGATIVE_INFINITY;
// int argmax = c.size() - 1;
// MVNormalDist mvnDist = new MVNormalDist();
// for (int i = 0; i < c.size(); i++) {
// double apost = 0;
// ChainLink cli = c.get(i);
// // Omega
// RandomVar<Double2D> rvOmega = (RandomVar<Double2D>) cli.get("Omega");
// Double2D Omega = rvOmega.getNumericValue();
// // apost += rvOmega.getPrior().logDensity(Omega);
// // M
// RandomVar<Double2D> rvM = (RandomVar<Double2D>) cli.get("M");
// Double2D M = rvM.getNumericValue();
// apost += rvM.getPrior().logDensity(M);
// // Hyperparameters
// RandomVar<Double0D> rvAl = (RandomVar<Double0D>) cli.get("al");
// Double0D al = rvAl.getNumericValue();
// apost += rvAl.getPrior().logDensity(al);
// // Model parameters and likelihood
// RandomVar<Double3D> rvSg = (RandomVar<Double3D>) cli.get("Sg");
// RandomVar<Double3D> rvA = (RandomVar<Double3D>) cli.get("A");
// // B
// RandomVar<Integer2D> rvB = (RandomVar<Integer2D>) cli.get("B");
// Integer2D B = rvB.getNumericValue();
// // skip
// // Z
// RandomVar<Integer1D> rvZ = (RandomVar<Integer1D>) cli.get("Z");
// Integer1D Z = rvZ.getNumericValue();
// for (int j = 0; j < Z.size(); j++) {
// int[] zIndices = new int[j];
// for (int k=0; k<j; k++) {
// zIndices[k] = k;
// }
// Integer1D ZPreJ = Z.getAll(zIndices);
// int nZ = ZPreJ.which(new Integer0D(Z.getValue(j))).size();
// if (nZ != 0) {
// apost += Math.log(nZ) - Math.log(ZPreJ.size() + al.value());
// } else {
// apost += Math.log(al.value())
// - Math.log(ZPreJ.size() + al.value());
// }
// }
// Integer1D active = Z.items();
// for (Integer0D z : active) {
// try {
// // Model parameters
// Double2D Sg = ((Double3D) rvSg.getNumericValue()).get(z.value());
// Double2D A = ((Double3D) rvA.getNumericValue()).get(z.value());
// CPDist postZ = (CPDist) rvZ.getPosterior();
// apost += postZ.getBaseSg().logDensity(Sg);
// apost += postZ.getBaseA().logDensity(A, M, Sg, postZ.getPhi());
// // Likelihood
// Integer1D whichz = Z.which(z);
// Integer2D db = B.getAll(whichz.value());
// Double2D dz = d.getAll(whichz.value());
// for (int j = 0; j < dz.size(); j++) {
// Double1D pt = dz.get(j);
// Integer1D b = db.get(j);
// Double1D mean = b.mult(A);
// Double2D cov = Sg;
// apost += mvnDist.logDensity(pt, mean, cov);
// }
// } catch (Exception e) {
// System.err.println(e.getMessage());
// StackTraceElement[] st = e.getStackTrace();
// for (int j = 0; j < st.length; j++) {
// System.err.println(st[j].toString());
// }
// throw new ProbDistParmException(e);
// }
// }
// // Set max
// if (apost > max) {
// argmax = i;
// max = apost;
// }
// }
// return c.get(argmax);
// }
//	
// public FLGFDModel(Map<String, Numeric<? extends Numeric<?>>> hypers,
// Map<String, Numeric<? extends Numeric<?>>> init, int dims)
// throws ProbDistParmException {
// super(hypers, dims);
// // Set up parameters
// this.params = new HashMap<String, RandomVar<? extends Numeric<?>>>();
//
// // Omega
// PriorOmegaDist priorOmega = new PriorOmegaDist((Double2D)
// this.getHyper("Phi"), (Double0D) this.getHyper("lambda"));
// PostOmegaDist postOmega = new PostOmegaDist((Double2D) this.getHyper("Phi"),
// (Double0D) this.getHyper("lambda"));
// Double2D defaultOmega = (Double2D) init.get("Omega"); // FIXME
// RandomVar<Double2D> rvOmega = new RandomVar<Double2D>("Omega", priorOmega,
// postOmega, defaultOmega);
// this.params.put("Omega", rvOmega);
//
// // Sg
// PriorSgDist priorSg = new PriorSgDist((Double2D) this.getHyper("Psi"),
// (Double0D) this.getHyper("kappa"));
// PostSgDist postSg = new PostSgDist((Double2D) this.getHyper("Psi"),
// (Double0D) this.getHyper("kappa"));
// Double3D defaultSg = (Double3D) init.get("Sg"); // FIXME
// RandomVar<Double3D> rvSg = new RandomVar<Double3D>("Sg", priorSg, postSg,
// defaultSg);
// this.params.put("Sg", rvSg);
//		
// // A
// PriorADist priorA = new PriorADist();
// PostADist postA = new PostADist();
// Double3D defaultA = (Double3D) init.get("A"); // FIXME
// RandomVar<Double3D> rvA = new RandomVar<Double3D>("A", priorA, postA,
// defaultA);
// this.params.put("A", rvA);
//
// // Z
// ProbDist<Integer1D> priorZ = new Deterministic<Integer1D>(init.get("Z"));
// CPDist postZ = new CPDist((Double2D) this.getHyper("Psi"),
// (Double0D) this.getHyper("kappa"), (Double2D) this.getHyper("Phi"),
// (Double0D) this.getHyper("lambda")); // FIXME -- Hack!
// postZ.setBaseSg(new InverseWishartDist()); // FIXME -- Eek! Hack!
// postZ.setBaseA(new MatrixNormalDist()); // FIXME -- Eek! Hack!
// Integer1D defaultZ = (Integer1D) init.get("Z");
// RandomVar<Integer1D> rvZ = new RandomVar<Integer1D>("Z", priorZ, postZ,
// defaultZ);
// this.params.put("Z", rvZ); // FIXME
//
// // B
// ProbDist<Integer2D> priorB = new Deterministic<Integer2D>(init.get("B"));
// PostBDist postB = new PostBDist();
// Integer2D defaultB = (Integer2D) init.get("B"); // FIXME
// RandomVar<Integer2D> rvB = new RandomVar<Integer2D>("B", priorB, postB,
// defaultB);
// this.params.put("B", rvB); // FIXME
//		
// // M
// Double2D defaultM = (Double2D) init.get("M"); // FIXME
// ProbDist<Double2D> priorM = new MatrixNormalDist((Double2D) this
// .getHyper("W"), (Double2D) this.getHyper("S"),
// Double2D.ident(defaultM.numRows()));
// PostMDist postM = new PostMDist((Double2D) this
// .getHyper("W"), (Double2D) this.getHyper("S"));
// RandomVar<Double2D> rvM = new RandomVar<Double2D>("M", priorM, postM,
// defaultM);
// this.params.put("M", rvM);
//
// // x
// ProbDist<Double0D> priorX = new BetaDist(
// (Double0D) this.getHyper("xa"), (Double0D) this.getHyper("xb"));
// PostXDist postX = new PostXDist();
// Double0D defaultX = (Double0D) init.get("x"); // FIXME
// RandomVar<Double0D> rvX = new RandomVar<Double0D>("x", priorX, postX,
// defaultX);
// this.params.put("x", rvX); // FIXME
//
// // al
// ProbDist<Double0D> priorAl = new GammaDist((Double0D) this
// .getHyper("ala"), (Double0D) this.getHyper("alb"));
// ProbDistInitializeByChain<Double0D> postAl = new PostAlDist((Double0D) this
// .getHyper("ala"), (Double0D) this.getHyper("alb"));
// Double0D defaultAl = (Double0D) init.get("al"); // FIXME
// RandomVar<Double0D> rvAl = new RandomVar<Double0D>("al", priorAl,
// postAl, defaultAl);
// this.params.put("al", rvAl); // FIXME
// }
//
// @Override
// public ChainLink getInitialLink() {
// return new ChainLink(this.getParam("Z"), this.getParam("B"),
// this.getParam("A"), this.getParam("Sg"), this.getParam("Omega"),
// this.getParam("M"), this.getParam("x"), this.getParam("al"));
// }
// }
//
// package org.jqgibbs.models;
//
// import java.util.ArrayList;
// import java.util.Arrays;
// import java.util.Collection;
// import java.util.HashMap;
// import java.util.List;
// import java.util.Map;
//
// import cern.colt.matrix.linalg.LUDecomposition;
// import cern.jet.stat.Gamma;
//
// import org.jqgibbs.ChainLink;
// import org.jqgibbs.Model;
// import org.jqgibbs.mathstat.AbstractSequence;
// import org.jqgibbs.mathstat.Double0D;
// import org.jqgibbs.mathstat.Double1D;
// import org.jqgibbs.mathstat.Double2D;
// import org.jqgibbs.mathstat.Double3D;
// import org.jqgibbs.mathstat.Integer0D;
// import org.jqgibbs.mathstat.Integer1D;
// import org.jqgibbs.mathstat.Integer2D;
// import org.jqgibbs.mathstat.ListSequence;
// import org.jqgibbs.mathstat.Numeric;
// import org.jqgibbs.mathstat.RandomVar;
// import org.jqgibbs.mathstat.UndirectedGraph;
// import org.jqgibbs.mathstat.probdist.BetaDist;
// import org.jqgibbs.mathstat.probdist.CategoricalDistInitializeByP;
// import org.jqgibbs.mathstat.probdist.GammaDist;
// import org.jqgibbs.mathstat.probdist.InverseWishartDist;
// import org.jqgibbs.mathstat.probdist.MVNormalDist;
// import org.jqgibbs.mathstat.probdist.MatrixNormalDist;
// import org.jqgibbs.mathstat.probdist.ProbDist;
// import org.jqgibbs.mathstat.probdist.ProbDistInitializeByChain;
// import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
// import org.jqgibbs.mathstat.probdist.ProbDistParmException;
// import org.jqgibbs.mathstat.probdist.SeqInverseWishartDist;
// import org.jqgibbs.mathstat.probdist.SeqMVNormalDist;
// import org.jqgibbs.mathstat.probdist.SeqMatrixNormalDist;
// import org.jqgibbs.mathstat.probdist.SeqWishartDist;
//
// public class FLGFDModel extends Model {
//	
// private static double gammaD(double x, int d) { // FIXME
// if (d == 1) {
// return Gamma.gamma(x);
// } else { // FIXME check d0
// return Math.pow(Math.PI, (d-1)/2)*FLGFDModel.gammaD(x,
// d-1)*Gamma.gamma(x+(1-d)/2);
// }
// }
//	
// class PriorOmegaDist extends ProbDistInitializeByChain<Double3D> {
// private SeqInverseWishartDist iwDists; // FIXME - Shouldn't this use the
// // IID?
// private ListSequence<Double2D> repPhi;
// private ListSequence<Double0D> repLambda;
//
// public PriorOmegaDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// private ListSequence<Double2D> getRepPhi() {
// if (this.repPhi == null) {
// this.repPhi = new ListSequence<Double2D>();
// Integer1D active = this.getZ().items();
// for (int i = 0; i < active.size(); i++) {
// this.repPhi.add(this.getPhi());
// }
// }
// return this.repPhi;
// }
//
// private ListSequence<Double0D> getRepLambda() {
// if (this.repLambda == null) {
// this.repLambda = new ListSequence<Double0D>();
// Integer1D active = this.getZ().items();
// for (int i = 0; i < active.size(); i++) {
// this.repLambda.add(this.getLambda());
// }
// }
// return this.repLambda;
// }
//
// private Double2D getPhi() {
// return (Double2D) this.fixedParms[0];
// }
//
// private Double0D getLambda() {
// return (Double0D) this.fixedParms[1];
// }
//
// private Integer1D getZ() {
// return (Integer1D) this.chainParms[0];
// }
//
// private Double3D iwVariates() throws ProbDistParmException {
// if (this.iwDists == null) {
// this.iwDists = new SeqInverseWishartDist(this.getRepPhi(), this
// .getRepLambda());
// return this.iwDists.variate();
// } else {
// return this.iwDists.variate(this.getRepPhi(), this
// .getRepLambda());
// }
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(1);
// this.chainParmNames.add("Z");
// this.fixedParmNames = new ArrayList<String>(2);
// this.fixedParmNames.add("Phi");
// this.fixedParmNames.add("lambda");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
// this.chainParmCheck.add(null); // FIXME
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 1);
// this.chainParmClasses.add(Integer1D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.fixedParmClasses.add(Double2D.class);
// this.fixedParmClasses.add(Double0D.class);
// }
//
// @Override
// protected void setUpFromChainParms() {
// return;
// }
//
// @Override
// protected Double3D genVariate() throws ProbDistParmException {
// return this.iwVariates();
// }
//
// @Override
// protected double getDensity(Double3D pt) {
// throw new UnsupportedOperationException("Too lazy, come back later");
// }
// }
//
// class PostOmegaDist extends ProbDistInitializeByChain<Double3D> {
// private ListSequence<Double2D> postPhis;
// private ListSequence<Double0D> postLambdas;
// private SeqInverseWishartDist iwDists;
//
// public PostOmegaDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// protected ListSequence<Double2D> getPostPhis() {
// return this.postPhis;
// }
//
// protected ListSequence<Double0D> getPostLambdas() {
// return this.postLambdas;
// }
//
// protected void setPostPhis(ListSequence<Double2D> postPhis) {
// this.postPhis = postPhis;
// }
//
// protected void setPostLambdas(ListSequence<Double0D> postLambdas) {
// this.postLambdas = postLambdas;
// }
//
// protected Double2D getPhi() {
// return (Double2D) this.fixedParms[0];
// }
//
// protected Double0D getLambda() {
// return (Double0D) this.fixedParms[1];
// }
//
// protected Integer1D getZ() {
// return (Integer1D) this.chainParms[0];
// }
//
// protected Double2D getM() {
// return (Double2D) this.chainParms[1];
// }
//
// protected Double3D getSg() {
// return (Double3D) this.chainParms[2];
// }
//
// protected Double3D getA() {
// return (Double3D) this.chainParms[3];
// }
//
// private Double3D iwVariates() throws ProbDistParmException {
// if (this.iwDists == null) {
// this.iwDists = new SeqInverseWishartDist(this.getPostPhis(),
// this.getPostLambdas());
// return this.iwDists.variate();
// } else {
// return this.iwDists.variate(this.getPostPhis(), this
// .getPostLambdas());
// }
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(4);
// this.chainParmNames.add("Z");
// this.chainParmNames.add("M");
// this.chainParmNames.add("Sg");
// this.chainParmNames.add("A");
// this.fixedParmNames = new ArrayList<String>(2);
// this.fixedParmNames.add("Phi");
// this.fixedParmNames.add("lambda");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 4);
// this.chainParmClasses.add(Integer1D.class);
// this.chainParmClasses.add(Double2D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double3D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.fixedParmClasses.add(Double2D.class);
// this.fixedParmClasses.add(Double0D.class);
// }
//
// @Override
// protected void setUpFromChainParms() {
// this.setPostPhis(new ListSequence<Double2D>());
// this.setPostLambdas(new ListSequence<Double0D>());
// Integer1D active = this.getZ().items();
// int d = FLGFDModel.this.dims;
// for (Integer0D k : active) {
// Double2D AM = this.getA().get(k.value()).minus(this.getM());
// Double2D sgInv = this.getSg().get(k.value()).inverse();
// this.getPostLambdas().add(this.getLambda().plus(d));
// this.getPostPhis()
// .add(
// this.getPhi().plus(
// AM.mult(sgInv).mult(AM.transpose())));
// }
// }
//
// @Override
// protected Double3D genVariate() throws ProbDistParmException {
// int h = this.getM().size();
// double[][][] v = new double[this.getA().size()][h][h];
// Double3D activeUpdates = this.iwVariates();
// int i = 0;
// Integer1D active = this.getZ().items();
// for (Integer0D k : active) {
// v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], h);
// i++;
// }
// return new Double3D(v);
// }
//
// @Override
// protected double getDensity(Double3D pt) throws ProbDistParmException {
// throw new UnsupportedOperationException("Too lazy, come back later");
// }
// }
//
// class PriorSgDist extends ProbDistInitializeByChain<Double3D> {
// private SeqInverseWishartDist iwDists; // FIXME - Shouldn't this use the
// // IID?
// private ListSequence<Double2D> repPsi;
// private ListSequence<Double0D> repKappa;
//
// public PriorSgDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// private ListSequence<Double2D> getRepPsi() {
// if (this.repPsi == null) {
// this.repPsi = new ListSequence<Double2D>();
// Integer1D active = this.getZ().items();
// for (int i = 0; i < active.size(); i++) {
// this.repPsi.add(this.getPsi());
// }
// }
// return this.repPsi;
// }
//
// private ListSequence<Double0D> getRepKappa() {
// if (this.repKappa == null) {
// this.repKappa = new ListSequence<Double0D>();
// Integer1D active = this.getZ().items();
// for (int i = 0; i < active.size(); i++) {
// this.repKappa.add(this.getKappa());
// }
// }
// return this.repKappa;
// }
//
// private Double2D getPsi() {
// return (Double2D) this.fixedParms[0];
// }
//
// private Double0D getKappa() {
// return (Double0D) this.fixedParms[1];
// }
//
// private Integer1D getZ() {
// return (Integer1D) this.chainParms[0];
// }
//
// private Double3D iwVariates() throws ProbDistParmException {
// if (this.iwDists == null) {
// this.iwDists = new SeqInverseWishartDist(this.getRepPsi(), this
// .getRepKappa());
// return this.iwDists.variate();
// } else {
// return this.iwDists.variate(this.getRepKappa(), this
// .getRepKappa());
// }
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(1);
// this.chainParmNames.add("Z");
// this.fixedParmNames = new ArrayList<String>(2);
// this.fixedParmNames.add("Psi");
// this.fixedParmNames.add("kappa");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
// this.chainParmCheck.add(null); // FIXME
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 1);
// this.chainParmClasses.add(Integer1D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.fixedParmClasses.add(Double2D.class);
// this.fixedParmClasses.add(Double0D.class);
// }
//
// @Override
// protected void setUpFromChainParms() {
// return;
// }
//
// @Override
// protected Double3D genVariate() throws ProbDistParmException {
// return this.iwVariates();
// }
//
// @Override
// protected double getDensity(Double3D pt) {
// throw new UnsupportedOperationException("Too lazy, come back later");
// }
// }
//
// class PostSgDist extends ProbDistInitializeByChain<Double3D> {
// private ListSequence<Double2D> postPsi;
// private ListSequence<Double0D> postKappa;
// private SeqInverseWishartDist iwDists;
//
// public PostSgDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// protected ListSequence<Double2D> getPostPsi() {
// return this.postPsi;
// }
//
// protected ListSequence<Double0D> getPostKappa() {
// return this.postKappa;
// }
//
// protected void setPostPsi(ListSequence<Double2D> postPsi) {
// this.postPsi = postPsi;
// }
//
// protected void setPostKappa(ListSequence<Double0D> postKappa) {
// this.postKappa = postKappa;
// }
//
// protected Double2D getPsi() {
// return (Double2D) this.fixedParms[0];
// }
//
// protected Double0D getKappa() {
// return (Double0D) this.fixedParms[1];
// }
//
// protected Integer1D getZ() {
// return (Integer1D) this.chainParms[0];
// }
//
// protected Integer2D getB() {
// return (Integer2D) this.chainParms[1];
// }
//
// protected Double3D getOmega() {
// return (Double3D) this.chainParms[2];
// }
//
// protected Double3D getA() {
// return (Double3D) this.chainParms[3];
// }
//
// protected Double2D getM() {
// return (Double2D) this.chainParms[4];
// }
//
// private Double3D iwVariates() throws ProbDistParmException {
// if (this.iwDists == null) {
// this.iwDists = new SeqInverseWishartDist(this.getPostPsi(),
// this.getPostKappa());
// return this.iwDists.variate();
// } else {
// return this.iwDists.variate(this.getPostPsi(), this
// .getPostKappa());
// }
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(5);
// this.chainParmNames.add("Z");
// this.chainParmNames.add("B");
// this.chainParmNames.add("Omega");
// this.chainParmNames.add("A");
// this.chainParmNames.add("M");
// this.fixedParmNames = new ArrayList<String>(2);
// this.fixedParmNames.add("Psi");
// this.fixedParmNames.add("kappa");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 5);
// this.chainParmClasses.add(Integer1D.class);
// this.chainParmClasses.add(Integer2D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double2D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.fixedParmClasses.add(Double2D.class);
// this.fixedParmClasses.add(Double0D.class);
// }
//
// @Override
// protected void setUpFromChainParms() {
// this.setPostPsi(new ListSequence<Double2D>());
// this.setPostKappa(new ListSequence<Double0D>());
// Integer1D active = this.getZ().items();
// int h = this.getA().size();
// for (Integer0D k : active) {
// Integer1D zK = this.getZ().which(k);
// Double2D xK = this.getSamplerData().getAll(zK.value());
// Integer0D nK = new Integer0D(xK.size());
//
// // this.getPostKappa().add(this.getKappa().plus(nK).plus(h));
// this.getPostKappa().add(this.getKappa().plus(nK));
//
// Integer2D B = this.getB().getAll(zK.value());
// Double2D BTB = B.transpose().mult(B);
// Double2D omegaInv = this.getOmega().get(k.value()).inverse();
// Double2D MTOmegaInv = this.getM().transpose().mult(omegaInv);
// Double2D MTOmegaInvM = MTOmegaInv.mult(this.getM());
// Double2D XTB = B.transpose().mult(xK).transpose();
// Double2D Mp = XTB.plus(MTOmegaInv);
// Double2D postOmega = omegaInv.plus(BTB).inverse();
// Double2D MpPOmegaMpT = Mp.mult(postOmega).mult(Mp.transpose());
// Double2D XTX = xK.transpose().mult(xK);
// this.getPostPsi().add(
// this.getPsi().plus(XTX).plus(MTOmegaInvM).minus(
// MpPOmegaMpT));
// }
// }
//
// @Override
// protected Double3D genVariate() throws ProbDistParmException {
// int d = FLGFDModel.this.dims;
// double[][][] v = new double[this.getA().size()][d][d];
// Double3D activeUpdates = this.iwVariates();
// int i = 0;
// Integer1D active = this.getZ().items();
// for (Integer0D k : active) {
// v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], d);
// i++;
// }
// return new Double3D(v);
// }
//
// @Override
// protected double getDensity(Double3D pt) throws ProbDistParmException {
// throw new UnsupportedOperationException("Too lazy, come back later");
// }
// }
//
// class PriorADist extends ProbDistInitializeByChain<Double3D> {
// private SeqMatrixNormalDist matnDists;
//
// public PriorADist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// private Double3D getActiveOmega() {
// Double3D omega = (Double3D) this.chainParms[0];
// Integer1D active = this.getZ().items();
// return omega.getAll(active.value());
// }
//
// private Double3D getActiveSg() {
// Double3D sg = (Double3D) this.chainParms[1];
// Integer1D active = this.getZ().items();
// return sg.getAll(active.value());
// }
//
// private Double3D getRepM() {
// Double2D M = (Double2D) this.chainParms[2];
// Integer1D active = this.getZ().items();
// Double3D repM = new Double3D(active.size(), this.getB().numCols(),
// FLGFDModel.this.dims);
// for (int i = 0; i < active.size(); i++) {
// repM.add(M);
// }
// return repM;
// }
//
// private Integer1D getZ() {
// return (Integer1D) this.chainParms[3];
// }
//
// private Integer2D getB() {
// return (Integer2D) this.chainParms[4];
// }
//
// private Double3D matnVariates() throws ProbDistParmException {
// if (this.matnDists == null) {
// this.matnDists = new SeqMatrixNormalDist(this.getRepM(), this
// .getActiveSg(), this.getActiveOmega());
// return this.matnDists.variate();
// } else {
// return this.matnDists.variate(this.getRepM(), this
// .getActiveSg(), this.getActiveOmega());
// }
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(5);
// this.chainParmNames.add("Omega");
// this.chainParmNames.add("Sg");
// this.chainParmNames.add("M");
// this.chainParmNames.add("Z");
// this.chainParmNames.add("B");
// this.fixedParmNames = new ArrayList<String>(0);
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
// this.chainParmCheck.add(null); // FIXME
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 5);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Integer1D.class);
// this.chainParmClasses.add(Integer2D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 0);
// }
//
// @Override
// protected void setUpFromChainParms() {
// return;
// }
//
// @Override
// protected Double3D genVariate() throws ProbDistParmException {
// return this.matnVariates();
// }
//
// @Override
// protected double getDensity(Double3D pt) {
// throw new UnsupportedOperationException("Too lazy, come back later");
// }
// }
//
// class PostADist extends ProbDistInitializeByChain<Double3D> {
// private ListSequence<Double2D> postMs;
// private ListSequence<Double2D> postSgs;
// private ListSequence<Double2D> postOmegas;
// private SeqMatrixNormalDist matnDists;
//
// public PostADist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// private Double3D matnVariates() throws ProbDistParmException {
// Integer1D active = this.getZ().items();
// Double3D sgActive = this.getSg().getAll(active.value());
// if (this.matnDists == null) {
// this.matnDists = new SeqMatrixNormalDist(this.getPostM(),
// sgActive, this.getPostOmega());
// // FIXME - Will this work if all arguments are not
// // ListSequences?
// return this.matnDists.variate();
// } else {
// return this.matnDists.variate(this.getPostM(), sgActive, this
// .getPostOmega());
// }
// }
//
// protected ListSequence<Double2D> getPostM() {
// return this.postMs;
// }
//
// protected ListSequence<Double2D> getPostOmega() {
// return this.postOmegas;
// }
//
// protected void setPostM(ListSequence<Double2D> postMs) {
// this.postMs = postMs;
// }
//
// protected void setPostOmega(ListSequence<Double2D> postOmegas) {
// this.postOmegas = postOmegas;
// }
//
// protected Double3D getOmega() {
// return (Double3D) this.chainParms[0];
// }
//
// protected Double3D getSg() {
// return (Double3D) this.chainParms[1];
// }
//
// protected Double2D getM() {
// return (Double2D) this.chainParms[2];
// }
//
// protected Integer1D getZ() {
// return (Integer1D) this.chainParms[3];
// }
//
// protected Integer2D getB() {
// return (Integer2D) this.chainParms[4];
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(5);
// this.chainParmNames.add("Omega");
// this.chainParmNames.add("Sg");
// this.chainParmNames.add("M");
// this.chainParmNames.add("Z");
// this.chainParmNames.add("B");
// this.fixedParmNames = new ArrayList<String>(0);
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
// this.chainParmCheck.add(null); // FIXME
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 5);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double2D.class);
// this.chainParmClasses.add(Integer1D.class);
// this.chainParmClasses.add(Integer2D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 0);
// }
//
// @Override
// protected void setUpFromChainParms() {
// this.setPostOmega(new ListSequence<Double2D>());
// this.setPostM(new ListSequence<Double2D>());
// Integer1D active = this.getZ().items();
// for (Integer0D k : active) {
// Integer1D zK = this.getZ().which(k);
// Double2D xK = this.getSamplerData().getAll(zK.value());
// Integer2D B = this.getB().getAll(zK.value());
// Double2D BTB = B.transpose().mult(B);
// Double2D omegaInv = this.getOmega().get(k.value()).inverse();
// Double2D postOmegaInv = omegaInv.plus(BTB);
// Double2D postOmega = postOmegaInv.inverse();
// this.getPostOmega().add(postOmega);
// Double2D BTX = B.transpose().mult(xK);
// Double2D omegaInvM = omegaInv.mult(this.getM());
// Double2D postM = postOmega.mult(BTX.plus(omegaInvM));
// this.getPostM().add(postM);
// }
// }
//
// @Override
// protected Double3D genVariate() throws ProbDistParmException {
// int K = this.getSg().size();
// int h = this.getB().numCols();
// int d = FLGFDModel.this.dims;
// double[][][] v = new double[K][h][d];
// boolean goodSample = false;
// while (!goodSample) {
// Double3D activeUpdates = this.matnVariates();
// goodSample = true;
// int i = 0;
// Integer1D active = this.getZ().items();
// for (Integer0D k : active) {
// v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], h); // FIXME??
// // Check
// // dims
// for (int m=0; m<h; m++) {
// for (int n=0; n<d; n++) {
// if (Double.isInfinite(v[k.value()][m][n]) ||
// Double.isNaN(v[k.value()][m][n])) {
// System.err.println("Warning: bad sample for A: " + activeUpdates.get(i));
// System.err.println("Sg: " + this.getSg().get(k.value()));
// System.err.println("Omega: " + this.getOmega().get(k.value()));
// System.err.println("M: " + this.getM());
// goodSample = false;
// break;
// }
// }
// if (!goodSample) {
// break;
// }
// }
// if (!goodSample) {
// break;
// }
// i++;
// }
// }
// return new Double3D(v);
// }
//
// @Override
// protected double getDensity(Double3D pt) {
// throw new UnsupportedOperationException("Too lazy, come back later");
// }
// }
//
//
// class CRPDist extends ProbDistInitializeByChain<Integer1D> {
// public CRPDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// protected Double1D postProb;
// private CategoricalDistInitializeByP catDist;
// private HashMap<Integer0D,MVNormalDist> mvnDists; // Shouldn't you be using
// (and hacking) SeqProbDist for this? FIXME
// private InverseWishartDist sgDist;
// private MatrixNormalDist aDist;
//		
// private Double3D savedOmegas;
//		
// private ProbDist<Double3D> baseOmegas;
// private ProbDist<Double2D> baseOmega;
// private ProbDist<Double2D> baseSg;
// private ProbDist<Double2D> baseA;
//		
// public ProbDist<Double2D> getBaseOmega() {
// return this.baseOmega;
// }
//		
// public void setBaseOmega(ProbDist<Double2D> baseOmega) {
// this.baseOmega = baseOmega;
// }
//
// public ProbDist<Double2D> getBaseSg() {
// return this.baseSg;
// }
//		
// public void setBaseSg(ProbDist<Double2D> baseSg) {
// this.baseSg = baseSg;
// }
//
// public ProbDist<Double2D> getBaseA() {
// return this.baseA;
// }
//		
// public void setBaseA(ProbDist<Double2D> baseA) {
// this.baseA = baseA;
// }
//		
// public ProbDist<Double3D> getBaseOmegas() throws ProbDistParmException {
// if (this.baseOmegas == null) {
// int h = this.getPhi().numCols();
// int N = this.getSamplerData().size();
// double[][][] phis = new double[N][h][h];
// double[] lambdas = new double[N];
// for (int i=0; i<N; i++) {
// phis[i] = Arrays.copyOf(this.getPhi().inverse().value(), h);
// lambdas[i] = this.getLambda().value();
// }
// this.baseOmegas = new SeqWishartDist(new Double3D(phis), new
// Double1D(lambdas));
// }
//			
// return this.baseOmegas;
// }
//
// protected Integer0D catVariate() throws ProbDistParmException {
// if (this.catDist == null) {
// this.catDist = new CategoricalDistInitializeByP(this
// .getPostProb());
// return this.catDist.variate();
// } else {
// return this.catDist.variate(this.getPostProb());
// }
// }
//
// protected double mvnLogDensity(Integer0D k, Double1D pt)
// throws ProbDistParmException {
// return this.mvnDists.get(k).logDensity(pt);
// }
//
// protected void setUpPostMvn(Integer0D k) throws ProbDistParmException {
// if (this.mvnDists == null) {
// this.mvnDists = new HashMap<Integer0D, MVNormalDist>();
// }
// if (!this.mvnDists.containsKey(k)) {
// double[] zeroD = new double[FLGFDModel.this.dims];
// Double1D zero = new Double1D(zeroD);
// this.mvnDists.put(k, new MVNormalDist(zero, this.getSg().get(
// k.value())));
// }
// // else do nothing
// }
//
// protected void resetPostMvn(Integer0D k) {
// if (this.mvnDists != null && this.mvnDists.containsKey(k)) {
// this.mvnDists.remove(k);
// }
// }
//
// private Double1D getPostProb() {
// return this.postProb;
// }
//
// protected void setPostProb(Double1D prob) {
// this.postProb = prob;
// }
//		
// protected double marginalNew(Double1D pt, Integer1D b, Double2D omegaInv) {
// // Double2D omegaN = this.getOmegaInv().plus(b.outer(b)).inverse();
// Double2D omegaN = omegaInv.plus(b.outer(b)).inverse();
// // Double2D mH = this.getM().transpose().mult(this.getOmegaInv());
// Double2D mH = this.getM().transpose().mult(omegaInv);
// Double2D mHX = mH.plus(pt.outer(b));
// // Double2D mO =
// this.getM().transpose().mult(this.getOmegaInv()).mult(this.getM());
// Double2D mO = this.getM().transpose().mult(omegaInv).mult(this.getM());
// Double2D mHXO = mHX.mult(omegaN).mult(mHX.transpose());
// Double2D psiN = this.getPsi().plus(pt.outer(pt).plus(mO).minus(mHXO));
// int d = this.getSamplerData().numCols();
// int h = this.getB().numCols();
// double kN = this.getKappa().value() + 1;
// double lmd = Math.log(FLGFDModel.gammaD(kN/2, d));
// lmd += Math.log(2)*kN*d/2;
// lmd += Math.log(2*Math.PI)*h*d/2;
// lmd += Math.log(omegaN.det())*d/2;
// lmd -= Math.log(psiN.det())*kN/2;
// return lmd;
// }
//		
// protected Double2D omegaVariate() throws ProbDistParmException {
// return this.getBaseOmega().variate();
// }
//		
// protected Double2D sgVariate(Double2D psi, double k) throws
// ProbDistParmException {
// if (this.sgDist == null) {
// this.sgDist = new InverseWishartDist(psi, new Double0D(k));
// return this.sgDist.variate();
// } else {
// return this.sgDist.variate(psi, new Double0D(k));
// }
// }
//		
// protected Double2D aVariate(Double2D M, Double2D sg, Double2D omega) throws
// ProbDistParmException {
// if (this.aDist == null) {
// this.aDist = new MatrixNormalDist(M, sg, omega);
// return this.aDist.variate();
// } else {
// return this.aDist.variate(M, sg, omega);
// }
// }
//		
// protected ListSequence<Double2D> getVariateDPPost(Double1D pt, Integer1D b)
// throws ProbDistParmException {
// Double2D omegaNew = this.omegaVariate();
// Double2D omegaInv = omegaNew.inverse();
// // Double2D omegaN = this.getOmega().inverse().plus(b.outer(b)).inverse();
// Double2D omegaN = omegaInv.plus(b.outer(b)).inverse();
// // Double2D mH = this.getM().transpose().mult(this.getOmega().inverse());
// Double2D mH = this.getM().transpose().mult(omegaInv);
// Double2D mHX = mH.plus(pt.outer(b));
// Double2D mHXO = mHX.mult(omegaN).mult(mHX.transpose());
// // Double2D mO =
// this.getM().transpose().mult(this.getOmega().inverse()).mult(this.getM());
// Double2D mO = this.getM().transpose().mult(omegaInv).mult(this.getM());
// Double2D psiN = this.getPsi().plus(pt.outer(pt).plus(mO).minus(mHXO));
// double kN = this.getKappa().value() + 1;
// Double2D sgNew = this.sgVariate(psiN, kN);
// Double2D aNew = this.aVariate(omegaN.mult(mHX.transpose()), sgNew, omegaN);
// return new ListSequence<Double2D>(omegaNew, sgNew, aNew);
// }
//		
// protected Double3D getOmega() {
// return (Double3D) this.chainParms[0];
// }
// //
// protected Double3D getSg() {
// return (Double3D) this.chainParms[1];
// // return (Double3D) this.chainParms[0];
// }
//
// protected Double3D getA() {
// return (Double3D) this.chainParms[2];
// // return (Double3D) this.chainParms[1];
// }
//
// protected Integer1D getZ() {
// return (Integer1D) this.chainParms[3];
// // return (Integer1D) this.chainParms[2];
// }
//
// protected Integer2D getB() {
// return (Integer2D) this.chainParms[4];
// // return (Integer2D) this.chainParms[3];
// }
//
// protected Double2D getM() {
// return (Double2D) this.chainParms[5];
// // return (Double2D) this.chainParms[4];
// }
//
// protected Double0D getAl() {
// return (Double0D) this.chainParms[6];
// // return (Double0D) this.chainParms[5];
// }
//
// protected Double2D getPsi() {
// return (Double2D) this.fixedParms[0];
// }
//
// protected Double0D getKappa() {
// return (Double0D) this.fixedParms[1];
// }
//
// protected Double0D getLambda() {
// return (Double0D) this.fixedParms[2];
// }
//		
// protected Double2D getPhi() {
// return (Double2D) this.fixedParms[3];
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(6);
// this.chainParmNames.add("Omega");
// this.chainParmNames.add("Sg");
// this.chainParmNames.add("A");
// this.chainParmNames.add("Z");
// this.chainParmNames.add("B");
// this.chainParmNames.add("M");
// this.chainParmNames.add("al");
// this.fixedParmNames = new ArrayList<String>(4);
// this.fixedParmNames.add("Psi");
// this.fixedParmNames.add("kappa");
// this.fixedParmNames.add("lambda");
// this.fixedParmNames.add("Phi");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(6);
// this.chainParmCheck.add(null); // FIXME
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(6);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Integer1D.class);
// this.chainParmClasses.add(Integer2D.class);
// this.chainParmClasses.add(Double2D.class);
// this.chainParmClasses.add(Double0D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(4);
// this.fixedParmClasses.add(Double2D.class);
// this.fixedParmClasses.add(Double0D.class);
// this.fixedParmClasses.add(Double0D.class);
// this.fixedParmClasses.add(Double2D.class);
// }
//
// @Override
// protected void setUpFromChainParms() {
// this.mvnDists = null;
// // this.omegaInv = this.getOmega().inverse();
// }
//
// @Override
// protected Integer1D genVariate() throws ProbDistParmException {
// // FIXME - This sampling scheme destructively modifies
// // A, Sg - HACK
// int N = this.getSamplerData().size();
// Integer1D Z = null;
// try {
// Z = (Integer1D) this.getZ().clone();
// } catch (CloneNotSupportedException e) {
// // TODO Auto-generated catch block
// e.printStackTrace(); // FIXME
// }
// Integer0D unused = new Integer0D(-1);
// if (this.savedOmegas == null) {
// this.savedOmegas = this.getBaseOmegas().variate();
// } else {
// this.savedOmegas = this.savedOmegas.cycle();
// }
// for (int i = 0; i < N; i++) {
// Double1D Xi = this.getSamplerData().get(i);
// Integer1D active = Z.items();
// Z.set(i, unused);
// double[] logP = new double[active.size()+1];
// double maxLogP = -Double.MAX_VALUE;
// int ai = 0;
// // Existing categories
// Integer1D b = this.getB().get(i);
// for (Integer0D k : active) {
// Integer1D zK = Z.which(k);
// int nZK = zK.size();
// if (nZK > 0) {
// Double1D aTB = b.mult(this.getA().get(k.value()));
// Double1D xIA = this.getSamplerData().get(i).minus(aTB);
// this.setUpPostMvn(k);
// logP[ai] = this.mvnLogDensity(k, xIA) + Math.log(nZK);
// if (logP[ai] == Double.POSITIVE_INFINITY) {
// logP[ai] = Double.MAX_VALUE;
// }
// if (logP[ai] > maxLogP) {
// maxLogP = logP[ai];
// }
// } else {
// logP[ai] = -Double.MAX_VALUE;
// }
// ai++;
// }
// // New category
// int newc = ai;
// logP[newc] = this.marginalNew(Xi, b, this.savedOmegas.get(i)) +
// Math.log(this.getAl().value());
// if (logP[ai] == Double.POSITIVE_INFINITY) {
// logP[ai] = Double.MAX_VALUE;
// }
// if (logP[newc] > maxLogP) {
// maxLogP = logP[newc];
// }
// // Select a category for this point
// Double1D p = (new Double1D(logP)).minus(maxLogP).exp();
// this.setPostProb(p);
// int az = this.catVariate().value();
// // Sample new category if necessary
// Integer0D zedNew = null;
// if (az == newc) {
// int maxActive = active.get(active.size() - 1).value();
// zedNew = Z.minNotIn(0, maxActive);
// this.resetPostMvn(zedNew);
// ListSequence<Double2D> newcParms = this.getVariateDPPost(Xi, b);
// Double2D newOmega = newcParms.get(0);
// Double2D newSg = newcParms.get(1);
// Double2D newA = newcParms.get(2);
// this.getA().set(zedNew.value(), newA);
// this.getSg().set(zedNew.value(), newSg);
// this.getOmega().set(zedNew.value(), newOmega);
// } else {
// zedNew = active.get(az);
// }
// Z.set(i, zedNew);
// }
// return Z;
// }
//
// @Override
// protected double getDensity(Integer1D pt) {
// throw new UnsupportedOperationException("Too lazy, come back later");
// }
// }
//
// class PostMDist extends ProbDistInitializeByChain<Double2D> {
// protected MVNormalDist mvnDist;
// protected Double1D vecW;
// protected Double2D SInv;
// protected Double1D postVecMu;
// protected Double2D postKronSg;
//
// public PostMDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// protected Integer1D getZ() {
// return (Integer1D) this.chainParms[0];
// }
//
// protected Integer2D getB() {
// return (Integer2D) this.chainParms[1];
// }
//
// protected Double3D getOmega() {
// return (Double3D) this.chainParms[2];
// }
//
// protected Double3D getSg() {
// return (Double3D) this.chainParms[3];
// }
//
// protected Double3D getA() {
// return (Double3D) this.chainParms[4];
// }
//
// protected Double2D getSInv() {
// if (this.SInv == null) {
// Double2D S = (Double2D) this.fixedParms[1];
// this.SInv = S.inverse();
// }
// return this.SInv;
// }
//
// protected Double1D getVecW() {
// if (this.vecW == null) {
// Double2D W = (Double2D) this.fixedParms[0];
// this.vecW = W.colVec();
// }
// return this.vecW;
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(5);
// this.chainParmNames.add("Z");
// this.chainParmNames.add("B");
// this.chainParmNames.add("Omega");
// this.chainParmNames.add("Sg");
// this.chainParmNames.add("A");
// this.fixedParmNames = new ArrayList<String>(2);
// this.fixedParmNames.add("W");
// this.fixedParmNames.add("S");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 5);
// this.chainParmClasses.add(Integer1D.class);
// this.chainParmClasses.add(Integer2D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double3D.class);
// this.chainParmClasses.add(Double3D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.fixedParmClasses.add(Double2D.class);
// this.fixedParmClasses.add(Double2D.class);
// }
//
// protected Double2D matnVariate() throws ProbDistParmException {
// int h = this.getA().numRows();
// Double1D vecVariate;
// if (this.mvnDist == null) {
// this.mvnDist = new MVNormalDist(this.getPostVecMu(), this
// .getPostKronSg());
// vecVariate = this.mvnDist.variate();
//
// } else {
// vecVariate = this.mvnDist.variate(this.getPostVecMu(), this
// .getPostKronSg());
// }
// return vecVariate.toDouble2D(h);
// }
//
// @Override
// protected void setUpFromChainParms() {
// int h = this.getB().numCols();
// Double2D sInvKIdh = this.getSInv().kron(Double2D.ident(h));
// Double2D postKronSgInv = sInvKIdh;
// Double1D postVecMup = sInvKIdh.mult(this.getVecW());
// Integer1D active = this.getZ().items();
// for (Integer0D k : active) {
// Double2D sgInv = this.getSg().get(k.value()).inverse();
// Double2D omegaInv = this.getOmega().get(k.value()).inverse();
// Double2D sgInvKOmegaInv = sgInv.kron(omegaInv);
// postKronSgInv = postKronSgInv.plus(sgInvKOmegaInv);
// Double1D vecA = this.getA().get(k.value()).colVec();
// postVecMup = postVecMup.plus(sgInvKOmegaInv.mult(vecA));
// }
// Double2D postKronSg = postKronSgInv.inverse();
// this.setPostKronSg(postKronSg);
// this.setPostVecMu(postKronSg.mult(postVecMup));
// }
//
// protected Double1D getPostVecMu() {
// return this.postVecMu;
// }
//
// protected void setPostVecMu(Double1D mu) {
// this.postVecMu = mu;
// }
//
// protected Double2D getPostKronSg() {
// return this.postKronSg;
// }
//
// protected void setPostKronSg(Double2D sg) {
// this.postKronSg = sg;
// }
//
// @Override
// protected Double2D genVariate() throws ProbDistParmException {
// return this.matnVariate();
// }
//
// @Override
// protected double getDensity(Double2D pt) {
// throw new UnsupportedOperationException("Too lazy, come back later");
// }
// }
//
// class PostXDist extends ProbDistInitializeByChain<Double0D> {
// private BetaDist betaDist;
// private Double0D postXA;
// private Double0D postXB;
//
// public PostXDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(1);
// this.chainParmNames.add("al");
// this.fixedParmNames = new ArrayList<String>(0);
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 1);
// this.chainParmClasses.add(Double0D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 0);
// }
//
// private Double0D betaVariate() throws ProbDistParmException {
// if (this.betaDist == null) {
// this.betaDist = new BetaDist(this.getPostXA(), this.getPostXB());
// return this.betaDist.variate();
// } else {
// return this.betaDist
// .variate(this.getPostXA(), this.getPostXB());
// }
// }
//
// protected Double0D getAl() {
// return (Double0D) this.chainParms[0];
// }
//
// @Override
// protected void setUpFromChainParms() {
// this.setPostXA(this.getAl().plus(1));
// this.setPostXB(new Double0D((double) this.getSamplerData().size()));
// }
//
// private Double0D getPostXA() {
// return this.postXA;
// }
//
// private void setPostXA(Double0D xA) {
// this.postXA = xA;
// }
//
// private Double0D getPostXB() {
// return this.postXB;
// }
//
// private void setPostXB(Double0D xB) {
// this.postXB = xB;
// }
//
// @Override
// protected Double0D genVariate() throws ProbDistParmException {
// return this.betaVariate();
// }
//
// @Override
// protected double getDensity(Double0D pt) {
// throw new UnsupportedOperationException("Too lazy, come back later");
// }
// }
//
// class PostAlDist extends ProbDistInitializeByChain<Double0D> {
// private GammaDist gammaDist;
// private CategoricalDistInitializeByP catDist;
// private Double0D postAlA;
// private Double0D postAlB;
// private Double1D postMix;
//
// public PostAlDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(2);
// this.chainParmNames.add("x");
// this.chainParmNames.add("Z");
// this.fixedParmNames = new ArrayList<String>(2);
// this.fixedParmNames.add("ala");
// this.fixedParmNames.add("alb");
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.chainParmCheck.add(null);
// this.chainParmCheck.add(null);
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
// this.fixedParmCheck.add(null);
// this.fixedParmCheck.add(null);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.chainParmClasses.add(Double0D.class);
// this.chainParmClasses.add(Integer1D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 2);
// this.fixedParmClasses.add(Double0D.class);
// this.fixedParmClasses.add(Double0D.class);
// }
//
// private Double0D gammaVariate() throws ProbDistParmException {
// if (this.gammaDist == null) {
// this.gammaDist = new GammaDist(this.getPostAlA(), this
// .getPostAlB());
// return this.gammaDist.variate();
// } else {
// return this.gammaDist.variate(this.getPostAlA(), this
// .getPostAlB());
// }
// }
//
// private Integer0D catVariate() throws ProbDistParmException {
// if (this.catDist == null) {
// this.catDist = new CategoricalDistInitializeByP(this
// .getPostMix());
// return this.catDist.variate();
// } else {
// return this.catDist.variate(this.getPostMix());
// }
// }
//
// protected Double0D getX() {
// return (Double0D) this.chainParms[0];
// }
//
// protected Integer1D getZ() {
// return (Integer1D) this.chainParms[1];
// }
//
// protected Double0D getAlA() {
// return (Double0D) this.fixedParms[0];
// }
//
// protected Double0D getAlB() {
// return (Double0D) this.fixedParms[1];
// }
//
// @Override
// protected void setUpFromChainParms() {
// int K = this.getZ().items().size();
// int N = this.getSamplerData().size();
// double logX = Math.log(this.getX().value());
// double p1 = this.getAlA().value() + K - 1;
// double p2 = N * (this.getAlB().value() - logX);
// this.setPostMix(new Double1D(p1, p2));
// int i;
// try {
// i = this.catVariate().value();
// } catch (ProbDistParmException e) {
// // FIXME
// throw new RuntimeException(e);
// }
// if (i == 0) {
// this.setPostAlA(this.getAlA().plus(K));
// } else {
// this.setPostAlA(this.getAlA().plus(K - 1));
// }
// this.setPostAlB(this.getAlB().plus(-logX)); // FIXME
// }
//
// private Double0D getPostAlA() {
// return this.postAlA;
// }
//
// private void setPostAlA(Double0D alA) {
// this.postAlA = alA;
// }
//
// private Double0D getPostAlB() {
// return this.postAlB;
// }
//
// private void setPostAlB(Double0D alB) {
// this.postAlB = alB;
// }
//
// private Double1D getPostMix() {
// return this.postMix;
// }
//
// private void setPostMix(Double1D mix) {
// this.postMix = mix;
// }
//
// @Override
// protected Double0D genVariate() throws ProbDistParmException {
// // FIXME - create
// return this.gammaVariate();
// }
//
// @Override
// protected double getDensity(Double0D pt) {
// throw new UnsupportedOperationException("Too lazy, come back later");
// }
// }
//
// class PostBDist extends ProbDistInitializeByChain<Integer2D> {
// public PostBDist(Numeric<?>... fixed) throws ProbDistParmException {
// super(fixed);
// }
//
// private Integer2D getB() {
// return (Integer2D) this.chainParms[0];
// }
//
// @Override
// protected void installParmChecks() {
// // Names
// this.chainParmNames = new ArrayList<String>(1);
// this.chainParmNames.add("B");
// this.fixedParmNames = new ArrayList<String>(0);
// // Checks
// this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
// this.chainParmCheck.add(null); // FIXME
// this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
// // Classes
// this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 1);
// this.chainParmClasses.add(Integer2D.class);
// this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
// 0);
// }
//
// @Override
// protected void setUpFromChainParms() {
// return;
// }
//
// @Override
// protected Integer2D genVariate() throws ProbDistParmException {
// return this.getB();
// }
//
// @Override
// protected double getDensity(Integer2D pt) {
// throw new UnsupportedOperationException("Too lazy, come back later");
// }
// }
//
// @SuppressWarnings("unchecked")
// public static ChainLink pointEstimate(List<ChainLink> c, Double2D d)
// throws ProbDistParmException {
// // MAP
// double max = Double.NEGATIVE_INFINITY;
// int argmax = c.size() - 1;
// MVNormalDist mvnDist = new MVNormalDist();
// for (int i = 0; i < c.size(); i++) {
// double apost = 0;
// ChainLink cli = c.get(i);
// // M
// RandomVar<Double2D> rvM = (RandomVar<Double2D>) cli.get("M");
// Double2D M = rvM.getNumericValue();
// apost += rvM.getPrior().logDensity(M);
// // Hyperparameters
// RandomVar<Double0D> rvAl = (RandomVar<Double0D>) cli.get("al");
// Double0D al = rvAl.getNumericValue();
// apost += rvAl.getPrior().logDensity(al);
// // Model parameters and likelihood
// RandomVar<Double3D> rvOmega = (RandomVar<Double3D>) cli
// .get("Omega");
// RandomVar<Double3D> rvSg = (RandomVar<Double3D>) cli.get("Sg");
// RandomVar<Double3D> rvA = (RandomVar<Double3D>) cli.get("A");
// // B
// RandomVar<Integer2D> rvB = (RandomVar<Integer2D>) cli.get("B");
// Integer2D B = rvB.getNumericValue();
// // skip
// // Z
// RandomVar<Integer1D> rvZ = (RandomVar<Integer1D>) cli.get("Z");
// Integer1D Z = rvZ.getNumericValue();
// for (int j = 0; j < Z.size(); j++) {
// int[] zIndices = new int[j];
// for (int k=0; k<j; k++) {
// zIndices[k] = k;
// }
// Integer1D ZPreJ = Z.getAll(zIndices);
// int nZ = ZPreJ.which(new Integer0D(Z.getValue(j))).size();
// if (nZ != 0) {
// apost += Math.log(nZ) - Math.log(ZPreJ.size() + al.value());
// } else {
// apost += Math.log(al.value())
// - Math.log(ZPreJ.size() + al.value());
// }
// }
// Integer1D active = Z.items();
// for (Integer0D z : active) {
// try {
// // Model parameters
// Double2D Omega = ((Double3D) rvOmega.getNumericValue()).get(z
// .value());
// Double2D Sg = ((Double3D) rvSg.getNumericValue())
// .get(z.value());
// Double2D A = ((Double3D) rvA.getNumericValue()).get(z.value());
// CRPDist postZ = (CRPDist) rvZ.getPosterior();
// apost += postZ.getBaseOmega().logDensity(Omega);
// apost += postZ.getBaseSg().logDensity(Sg);
// apost += postZ.getBaseA().logDensity(A, M, Sg, Omega);
// // Likelihood
// Integer1D whichz = Z.which(z);
// Integer2D db = B.getAll(whichz.value());
// Double2D dz = d.getAll(whichz.value());
// for (int j = 0; j < dz.size(); j++) {
// Double1D pt = dz.get(j);
// Integer1D b = db.get(j);
// Double1D mean = b.mult(A);
// Double2D cov = Sg;
// apost += mvnDist.logDensity(pt, mean, cov);
// }
// } catch (Exception e) {
//					
// System.err.println(e.getMessage());
// StackTraceElement[] st = e.getStackTrace();
// for (int j = 0; j < st.length; j++) {
// System.err.println(st[j].toString());
// }
// throw new ProbDistParmException(e);
// }
// }
// // Set max
// if (apost > max) {
// argmax = i;
// max = apost;
// }
// }
// return c.get(argmax);
// }
//
// // // Set up matrix for pairs
// // int n = d.size();
// // int nPairs = n * n - n * (n + 1) / 2;
// // Integer2D pairsSeq = new Integer2D(c.size(), nPairs);
// // // Set up matrix for B
// // ListSequence<Integer2D> BSeq = new ListSequence<Integer2D>();
// // // For each link in the chain,
// // for (int i = 0; i < c.size(); i++) {
// // // Save pairs
// // pairsSeq.set(i, ((Integer1D) c.get(i).get("Z").getNumericValue())
// // .samePairs());
// // // Save B
// // BSeq.set(i, (Integer2D) c.get(i).get("B").getNumericValue());
// // }
// // // Take consistent pairs
// // Double1D ppair = pairsSeq.mean();
// // int[] pairsi = new int[nPairs];
// // // CategoricalDistInitializeByP dist = new
// // // CategoricalDistInitializeByP(new Double1D(0.5, 0.5));
// // for (int i = 0; i < nPairs; i++) {
// // double p = ppair.get(i).value();
// // // pairsi[i] = dist.variateFast(new Double1D((1-p), p)).value();
// // if (p > 0.75) {
// // pairsi[i] = 1;
// // }
// // // pairsi[i] = dist.variateFast(new Double1D((1-p), p)).value();
// // }
// // Integer1D pairs = new Integer1D(pairsi);
// // // Take mode B
// // Integer2D B = ((Integer2D) c.get(c.size() - 1).get("B")
// // .getNumericValue()).cloneFromVector(BSeq.modeVector());
// // // Get point estimate for Z by doing assignment to graph components
// // Integer1D Z = new UndirectedGraph(pairs).components();
// // // Initialize
// // ChainLink last = c.get(c.size() - 1);
// // int nComps = Z.max().value() + 1;
// // int nObs = ((Double3D) last.get("A").getNumericValue()).get(0)
// // .numRows();
// // int nDims = ((Double3D) last.get("A").getNumericValue()).get(0)
// // .numCols();
// // // Z
// // RandomVar<Integer1D> rvZ = (RandomVar<Integer1D>) last.get("Z")
// // .cloneFromVector(Z.rowVec());
// // // B
// // RandomVar<Integer1D> rvB = (RandomVar<Integer1D>) last.get("B")
// // .cloneFromVector(B.rowVec());
// // // M
// // Double3D Ms = new Double3D(c.size(), nObs, nDims);
// // for (int i = 0; i < c.size(); i++) {
// // Ms.set(i, (Double2D) c.get(i).get("M").getNumericValue());
// // }
// // RandomVar<Double2D> rvM = (RandomVar<Double2D>) last.get("M")
// // .cloneFromVector(Ms.meanVector());
// // Double2D M = rvM.getNumericValue();
// // // Model parameters
// // Double3D AInit = new Double3D(nComps, nObs, nDims);
// // Double3D SgInit = new Double3D(nComps, nDims, nDims); //
// // Double3D OmegaInit = new Double3D(nComps, nObs, nObs); //
// // for (Integer0D z : Z.items()) {
// // Double2D Omegai = ((CPDist) rvZ.getPosterior()).getBaseOmega()
// // .variateFast();
// // Double2D Sgi = ((CPDist) rvZ.getPosterior()).getBaseSg()
// // .variateFast();
// // Double2D Ai = ((CPDist) rvZ.getPosterior()).getBaseA().variateFast(
// // M, Sgi, Omegai);
// // SgInit.set(z.value(), (Double2D) Sgi);
// // OmegaInit.set(z.value(), (Double2D) Omegai);
// // AInit.set(z.value(), (Double2D) Ai);
// // }
// // RandomVar<Double3D> rvSg= new RandomVar<Double3D>("Sg",
// // ((RandomVar<Double3D>) last.get("Sg")).getPrior(),
// // ((RandomVar<Double3D>) last.get("Sg")).getPosterior(), SgInit);
// // RandomVar<Double3D> rvOmega = new RandomVar<Double3D>("Omega",
// // ((RandomVar<Double3D>) last.get("Omega")).getPrior(),
// // ((RandomVar<Double3D>) last.get("Omega")).getPosterior(),
// // OmegaInit);
// // RandomVar<Double3D> rvA = new RandomVar<Double3D>("A",
// // ((RandomVar<Double3D>) last.get("A")).getPrior(),
// // ((RandomVar<Double3D>) last.get("A")).getPosterior(), AInit);
// // // DP hyperparameters
// // RandomVar<Double0D> rvX = (RandomVar<Double0D>) last.get("x");
// // RandomVar<Double0D> rvAl = (RandomVar<Double0D>) last.get("al");
// // // Sample new model parameters using current prior, observing Z and B
// // Map<String, AbstractSequence<? extends AbstractSequence<?, ?>, ? extends
// // Numeric<?>>> vars = new HashMap<String, AbstractSequence<? extends
// // AbstractSequence<?, ?>, ? extends Numeric<?>>>();
// // vars.put("A", AInit.sequence());
// // vars.put("Sg", SgInit.sequence());
// // vars.put("Omega", OmegaInit.sequence());
// // ChainLink curr = new ChainLink(rvZ, rvB, rvA, rvSg, rvOmega, rvM, rvX,
// // rvAl);
// // ChainLink next = null;
// // // Now draw a big sample
// // for (int i = 0; i < 2*c.size(); i++) {
// // try {
// // next = (ChainLink) curr.clone();
// // next.get("A").updatePosteriorFast(next, d);
// // next.get("Sg").updatePosteriorFast(next, d);
// // next.get("Omega").updatePosteriorFast(next, d);
// // if (i >= c.size()) {
// // ((ListSequence<Double3D>) vars.get("A"))
// // .add((Double3D) next.get("A").getNumericValue()
// // .clone());
// // ((ListSequence<Double3D>) vars.get("Sg"))
// // .add((Double3D) next.get("Sg").getNumericValue()
// // .clone());
// // ((ListSequence<Double3D>) vars.get("Omega"))
// // .add((Double3D) next.get("Omega").getNumericValue()
// // .clone());
// // }
// // curr = next;
// // } catch (CloneNotSupportedException e) {
// // throw new RuntimeException(e);
// // } catch (Exception e) {
// // System.err.println(e.getMessage());
// // StackTraceElement[] st = e.getStackTrace();
// // for (int j = 0; j < st.length; j++) {
// // System.err.println(st[j].toString());
// // }
// // throw new ProbDistParmException(e);
// // }
// // }
// // // Construct sample mean of model parameters
// // rvSg = (RandomVar<Double3D>) next.get("Sg")
// // .cloneFromVector(vars.get("Sg").meanVector());
// // rvOmega = (RandomVar<Double3D>) next.get("Omega")
// // .cloneFromVector(vars.get("Omega").meanVector());
// // rvA = (RandomVar<Double3D>) next.get("A")
// // .cloneFromVector(vars.get("A").meanVector());
// // // rvM = (RandomVar<Double2D>)
// // next.get("M").cloneFromVector(vars.get("M").meanVector());
// // ChainLink pte = new ChainLink(rvZ, rvB, rvA, rvSg, rvOmega, rvM, rvX,
// // rvAl);
// // // Now resample Z a few times to make the spurious categories go away
// // for (int i = 0; i < 100; i++) {
// // pte.get("Z").updatePosteriorFast(pte, d);
// // }
// // return pte;
// // }
//
// public FLGFDModel() {
// super();
// }
//
// public FLGFDModel(Map<String, Numeric<? extends Numeric<?>>> hypers,
// Map<String, Numeric<? extends Numeric<?>>> init, int dims)
// throws ProbDistParmException {
// super(hypers, dims);
// // Set up parameters
// this.params = new HashMap<String, RandomVar<? extends Numeric<?>>>();
//
// // Omega
// PriorOmegaDist priorOmega = new PriorOmegaDist((Double2D) this
// .getHyper("Phi"), (Double0D) this.getHyper("lambda"));
// PostOmegaDist postOmega = new PostOmegaDist((Double2D) this
// .getHyper("Phi"), (Double0D) this.getHyper("lambda"));
// Double3D defaultOmega = (Double3D) init.get("Omega"); // FIXME
// RandomVar<Double3D> rvOmega = new RandomVar<Double3D>("Omega",
// priorOmega, postOmega, defaultOmega);
// this.params.put("Omega", rvOmega);
//
// // Sg
// PriorSgDist priorSg = new PriorSgDist((Double2D) this.getHyper("Psi"),
// (Double0D) this.getHyper("kappa"));
// PostSgDist postSg = new PostSgDist((Double2D) this.getHyper("Psi"),
// (Double0D) this.getHyper("kappa"));
// Double3D defaultSg = (Double3D) init.get("Sg"); // FIXME
// RandomVar<Double3D> rvSg = new RandomVar<Double3D>("Sg", priorSg,
// postSg, defaultSg);
// this.params.put("Sg", rvSg);
//
// // A
// PriorADist priorA = new PriorADist();
// PostADist postA = new PostADist();
// Double3D defaultA = (Double3D) init.get("A"); // FIXME
// RandomVar<Double3D> rvA = new RandomVar<Double3D>("A", priorA, postA,
// defaultA);
// this.params.put("A", rvA);
//
// // Z
// ProbDist<Integer1D> priorZ = new Deterministic<Integer1D>(init.get("Z"));
// CRPDist postZ = new CRPDist((Double2D) this.getHyper("Psi"),
// (Double0D) this.getHyper("kappa"), (Double0D) this.getHyper("lambda"),
// (Double2D) this.getHyper("Phi")); // FIXME
// // --
// // Hack!
// postZ.setBaseOmega(new InverseWishartDist(this.getHyper("Phi"),
// this.getHyper("lambda")));
// postZ.setBaseSg(new InverseWishartDist(this.getHyper("Psi"), this
// .getHyper("kappa"))); // FIXME -- Eek! Hack!
// postZ.setBaseA(new MatrixNormalDist()); // FIXME -- Eek! Hack!
// Integer1D defaultZ = (Integer1D) init.get("Z");
// RandomVar<Integer1D> rvZ = new RandomVar<Integer1D>("Z", priorZ, postZ,
// defaultZ);
// this.params.put("Z", rvZ); // FIXME
//
// // B
// ProbDist<Integer2D> priorB = new Deterministic<Integer2D>(init.get("B"));
// PostBDist postB = new PostBDist();
// Integer2D defaultB = (Integer2D) init.get("B"); // FIXME
// RandomVar<Integer2D> rvB = new RandomVar<Integer2D>("B", priorB, postB,
// defaultB);
// this.params.put("B", rvB); // FIXME
//
// // M
// Double2D defaultM = (Double2D) init.get("M"); // FIXME
// ProbDist<Double2D> priorM = new MatrixNormalDist((Double2D) this
// .getHyper("W"), (Double2D) this.getHyper("S"), Double2D
// .ident(defaultM.numRows()));
// PostMDist postM = new PostMDist((Double2D) this.getHyper("W"),
// (Double2D) this.getHyper("S"));
// RandomVar<Double2D> rvM = new RandomVar<Double2D>("M", priorM, postM,
// defaultM);
// this.params.put("M", rvM);
//
// // x
// ProbDist<Double0D> priorX = new BetaDist(
// (Double0D) this.getHyper("xa"), (Double0D) this.getHyper("xb"));
// PostXDist postX = new PostXDist();
// Double0D defaultX = (Double0D) init.get("x"); // FIXME
// RandomVar<Double0D> rvX = new RandomVar<Double0D>("x", priorX, postX,
// defaultX);
// this.params.put("x", rvX); // FIXME
//
// // al
// ProbDist<Double0D> priorAl = new GammaDist((Double0D) this
// .getHyper("ala"), (Double0D) this.getHyper("alb"));
// ProbDistInitializeByChain<Double0D> postAl = new PostAlDist(
// (Double0D) this.getHyper("ala"), (Double0D) this
// .getHyper("alb"));
// Double0D defaultAl = (Double0D) init.get("al"); // FIXME
// RandomVar<Double0D> rvAl = new RandomVar<Double0D>("al", priorAl,
// postAl, defaultAl);
// this.params.put("al", rvAl); // FIXME
// }
//
// @Override
// public ChainLink getInitialLink() {
// return new ChainLink(this.getParam("Z"), this.getParam("B"), this
// .getParam("A"), this.getParam("Sg"), this.getParam("Omega"),
// this.getParam("M"), this.getParam("x"), this.getParam("al"));
// }
// }
