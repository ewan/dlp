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
import org.jqgibbs.mathstat.probdist.SeqMatrixNormalDist;
import cern.colt.matrix.linalg.SingularValueDecomposition;

import cern.jet.stat.Gamma;

public class MVRModel extends Model {

	class PostOmegaDist extends ProbDistMC<Double2D> {
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
		
		protected Double2D getM() {
			return (Double2D) this.fixedParms[2];
		}

		protected Double2D getSg() {
			return (Double2D) this.chainParms[0];
		}

		protected Double2D getA() {
			return (Double2D) this.chainParms[1];
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
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("A");
			this.fixedParmNames = new ArrayList<String>(3);
			this.fixedParmNames.add("Phi");
			this.fixedParmNames.add("lambda");
			this.fixedParmNames.add("M");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(3);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(3);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					3);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double2D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			Double2D postPhi = this.getPhi();
			Double0D postLambda = this.getLambda();
			int d = MVRModel.this.dims;
			Double2D AM = this.getA().minus(this.getM());
			Double2D sgInv = this.getSg().inverse();
			postLambda = postLambda.plus(d);
			postPhi = postPhi.plus(AM.mult(sgInv).mult(AM.transpose()));
			this.setPostPhi(postPhi);
			this.setPostLambda(postLambda);
		}

		@Override
		protected Double2D genVariate() throws ProbDistParmException {
			return this.iwVariate();
		}

		@Override
		protected double getDensity(Double2D pt) throws ProbDistParmException {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PriorSgDist extends ProbDistMC<Double2D> {
		private InverseWishartDist iwDist;

		public PriorSgDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private Double2D getPsi() {
			return (Double2D) this.fixedParms[0];
		}

		private Double0D getKappa() {
			return (Double0D) this.fixedParms[1];
		}

		private Double2D iwVariate() throws ProbDistParmException {
			if (this.iwDist == null) {
				this.iwDist = new InverseWishartDist(this.getPsi(), this.getKappa());
				return this.iwDist.variate();
			} else {
				return this.iwDist.variate(this.getKappa(), this.getKappa());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(0);
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("Psi");
			this.fixedParmNames.add("kappa");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(0);
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
		protected Double2D genVariate() throws ProbDistParmException {
			return this.iwVariate();
		}

		@Override
		protected double getDensity(Double2D pt) {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PostSgDist extends ProbDistMC<Double2D> {
		private Double2D postPsi;
		private Double0D postKappa;
		private InverseWishartDist iwDist;

		public PostSgDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double2D getPostPsi() {
			return this.postPsi;
		}

		protected Double0D getPostKappa() {
			return this.postKappa;
		}

		protected void setPostPsi(Double2D postPsi) {
			this.postPsi = postPsi;
		}

		protected void setPostKappa(Double0D postKappa) {
			this.postKappa = postKappa;
		}

		protected Double2D getPsi() {
			return (Double2D) this.fixedParms[0];
		}

		protected Double0D getKappa() {
			return (Double0D) this.fixedParms[1];
		}

		protected Double2D getM() {
			return (Double2D) this.fixedParms[2];
		}

		protected Double2D getB() {
			return (Double2D) this.chainParms[0];
		}

		protected Double2D getOmega() {
			return (Double2D) this.chainParms[1];
		}

		protected Double3D getA() {
			return (Double3D) this.chainParms[2];
		}

		private Double2D iwVariate() throws ProbDistParmException {
			if (this.iwDist == null) {
				this.iwDist = new InverseWishartDist(this.getPostPsi(), this.getPostKappa());
				return this.iwDist.variate();
			} else {
				return this.iwDist.variate(this.getPostPsi(), this.getPostKappa());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(3);
			this.chainParmNames.add("B");
			this.chainParmNames.add("Omega");
			this.chainParmNames.add("A");
			this.fixedParmNames = new ArrayList<String>(3);
			this.fixedParmNames.add("Psi");
			this.fixedParmNames.add("kappa");
			this.fixedParmNames.add("M");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(3);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(3);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(3);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					3);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double2D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			Double2D omegaInv = this.getOmega().inverse();
			Double2D X = this.getSamplerData();
			Integer0D N = new Integer0D(X.size());
			this.setPostKappa(this.getKappa().plus(N));
			Double2D B = this.getB();
			Double2D BTB = B.transpose().mult(B);
			Double2D MTOmegaInv = this.getM().transpose().mult(omegaInv);
			Double2D MTOmegaInvM = MTOmegaInv.mult(this.getM());
			Double2D XTB = B.transpose().mult(X).transpose();
			Double2D Mp = XTB.plus(MTOmegaInv);
			Double2D postOmega = omegaInv.plus(BTB).inverse();
			Double2D MpPOmegaMpT = Mp.mult(postOmega).mult(Mp.transpose());
			Double2D XTX = X.transpose().mult(X);
			this.setPostPsi(this.getPsi().plus(XTX).plus(MTOmegaInvM).minus(MpPOmegaMpT));
		}

		@Override
		protected Double2D genVariate() throws ProbDistParmException {
			return this.iwVariate();
		}

		@Override
		protected double getDensity(Double2D pt) throws ProbDistParmException {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PriorADist extends ProbDistMC<Double2D> {
		private MatrixNormalDist matnDist;

		public PriorADist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private Double2D getOmega() {
			return (Double2D) this.chainParms[0];
		}

		private Double2D getSg() {
			return (Double2D) this.chainParms[1];
		}

		private Double2D getM() {
			return (Double2D) this.fixedParms[0];
		}

		private Double2D matnVariate() throws ProbDistParmException {
			if (this.matnDist == null) {
				this.matnDist = new MatrixNormalDist(this.getM(), this.getSg(), this.getOmega());
				return this.matnDist.variate();
			} else {
				return this.matnDist.variate(this.getM(), this.getSg(), this.getOmega());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(4);
			this.chainParmNames.add("Omega");
			this.chainParmNames.add("Sg");
			this.fixedParmNames = new ArrayList<String>(1);
			this.fixedParmNames.add("M");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					4);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					1);
			this.fixedParmClasses.add(Double2D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
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

	class PostADist extends ProbDistMC<Double2D> {
		private Double2D postM;
		private Double2D postOmega;
		private MatrixNormalDist matnDist;

		public PostADist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private Double2D matnVariate() throws ProbDistParmException {
			if (this.matnDist == null) {
				this.matnDist = new MatrixNormalDist(this.getPostM(),
						this.getSg(), this.getPostOmega());
				return this.matnDist.variate();
			} else {
				return this.matnDist.variate(this.getPostM(), this.getSg(), this
						.getPostOmega());
			}
		}

		protected Double2D getPostM() {
			return this.postM;
		}

		protected Double2D getPostOmega() {
			return this.postOmega;
		}

		protected void setPostM(Double2D postM) {
			this.postM = postM;
		}

		protected void setPostOmega(Double2D postOmega) {
			this.postOmega = postOmega;
		}

		protected Double2D getOmega() {
			return (Double2D) this.chainParms[0];
		}

		protected Double2D getSg() {
			return (Double2D) this.chainParms[1];
		}
		
		protected Double2D getB() {
			return (Double2D) this.chainParms[2];
		}

		protected Double2D getM() {
			return (Double2D) this.fixedParms[0];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(5);
			this.chainParmNames.add("Omega");
			this.chainParmNames.add("Sg");
			this.chainParmNames.add("B");
			this.fixedParmNames = new ArrayList<String>(1);
			this.fixedParmNames.add("M");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					5);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(1);
			this.fixedParmClasses.add(Double2D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			Double2D omegaInv = this.getOmega().inverse();
			Double2D X = this.getSamplerData();
			Double2D B = this.getB();
			Double2D BTB = B.transpose().mult(B);
			Double2D postOmegaInv = omegaInv.plus(BTB);
			Double2D postOmega = postOmegaInv.inverse();
			this.setPostOmega(postOmega);
			Double2D BTX = B.transpose().mult(X);
			Double2D omegaInvM = omegaInv.mult(this.getM());
			Double2D postM = postOmega.mult(BTX.plus(omegaInvM));
			this.setPostM(postM);
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

	class PostBDist extends ProbDistMC<Double2D> {
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
			this.chainParmCheck.add(null); 
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
	public static ChainLink pointEstimate(List<ChainLink> c, Double2D d) throws ProbDistParmException {
		// MAP
		double max = Double.NEGATIVE_INFINITY;
		int argmax = c.size() - 1;
		MVNormalDist mvnDist = new MVNormalDist();
		for (int i = 0; i < c.size(); i++) {
			double apost = 0;
			ChainLink cli = c.get(i);
			// Omega
			RandomVar<Double2D> rvOmega = (RandomVar<Double2D>) cli.get("Omega");
			Double2D Omega = rvOmega.getNumericValue();
			apost += rvOmega.getPrior().logDensity(Omega);
			// Model parameters and likelihood
			RandomVar<Double2D> rvSg = (RandomVar<Double2D>) cli.get("Sg");
			RandomVar<Double2D> rvA = (RandomVar<Double2D>) cli.get("A");
			// Pull out B (not really random)
			RandomVar<Double2D> rvB = (RandomVar<Double2D>) cli.get("B");
			Double2D B = rvB.getNumericValue();
			// Model parameters
			Double2D Sg = (Double2D) rvSg.getNumericValue();
			Double2D A = (Double2D) rvA.getNumericValue();
			// Likelihood
			for (int j = 0; j < d.size(); j++) {
				Double1D pt = d.get(j);
				Double1D b = B.get(j);
				Double1D mean = b.mult(A);
				apost += mvnDist.logDensity(pt, mean, Sg);
			}
			// Set max
			if (apost > max) {
				argmax = i;
				max = apost;
			}
		}
		return c.get(argmax);
	}

	public MVRModel() {
		super();
	}

	public MVRModel(Map<String, Numeric<? extends Numeric<?>>> hypers,
			Map<String, Numeric<? extends Numeric<?>>> init, int dims)
			throws ProbDistParmException {
		super(hypers, dims);
		// Set up parameters
		this.params = new HashMap<String, RandomVar<? extends Numeric<?>>>();

		// Omega
		ProbDist<Double2D> priorOmega = new InverseWishartDist(this
				.getHyper("Phi"), (Double0D) this.getHyper("lambda"));
		PostOmegaDist postOmega = new PostOmegaDist((Double2D) this
				.getHyper("Phi"), (Double0D) this.getHyper("lambda"),
				(Double2D) this.getHyper("M"));
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
				(Double0D) this.getHyper("kappa"), (Double2D) this.getHyper("M"));
		Double2D defaultSg = (Double2D) init.get("Sg"); // FIXME
		RandomVar<Double2D> rvSg = new RandomVar<Double2D>("Sg", priorSg,
				postSg, defaultSg);
		this.params.put("Sg", rvSg);

		// A
		PriorADist priorA = new PriorADist((Double2D) this.getHyper("M"));
		PostADist postA = new PostADist((Double2D) this.getHyper("M"));
		Double2D defaultA = (Double2D) init.get("A"); // FIXME
		RandomVar<Double2D> rvA = new RandomVar<Double2D>("A", priorA, postA,
				defaultA);
		this.params.put("A", rvA);

		// B
		ProbDist<Double2D> priorB = new Deterministic<Double2D>(init.get("B"));
		PostBDist postB = new PostBDist();
		Double2D defaultB = (Double2D) init.get("B"); // FIXME
		RandomVar<Double2D> rvB = new RandomVar<Double2D>("B", priorB, postB,
				defaultB);
		this.params.put("B", rvB); // FIXME
	}

	@Override
	public ChainLink getInitialLink() {
		return new ChainLink(this.getParam("B"), this.getParam("A"), this.getParam("Sg"), this.getParam("Omega"));
	}
}
