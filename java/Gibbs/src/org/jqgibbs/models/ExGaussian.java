package org.jqgibbs.models;

import java.util.ArrayList;
import java.util.Arrays;
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
//import org.jqgibbs.mathstat.ListSequence;
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.mathstat.RandomVar;
import org.jqgibbs.mathstat.probdist.GammaDist;
import org.jqgibbs.mathstat.probdist.MVNormalDist;
import org.jqgibbs.mathstat.probdist.MatrixNormalDist;
import org.jqgibbs.mathstat.probdist.NormalDist;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistInitializeByChain;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;
import org.jqgibbs.mathstat.probdist.SeqInverseWishartDist;
import org.jqgibbs.mathstat.probdist.TruncatedNormalDist;
import org.jqgibbs.models.FLGFAModel.PostOmegaDist;
import org.jqgibbs.models.FLGFAModel.PriorOmegaDist;

public class ExGaussian extends Model {

	class PriorYDist extends ProbDistInitializeByChain<Double1D> {
		private GammaDist gammaDist;

		public PriorYDist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double0D getLambda() {
			return (Double0D) this.chainParms[0];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(1);
			this.chainParmNames.add("lambda");
			this.fixedParmNames = new ArrayList<String>(0);
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					1);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
					0);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		private Double1D expVariates() throws ProbDistParmException {
			int N = this.getSamplerData().size();
			if (this.gammaDist == null) {
				this.gammaDist = new GammaDist(new Double0D(1), this
						.getLambda());
				return (Double1D) this.gammaDist.variatesIID(N); // FIXME -
				// probably
				// unsafe
			} else {
				return (Double1D) this.gammaDist.variatesIID(N,
						new Double0D(1), this.getLambda());
			}
		}

		@Override
		protected Double1D genVariate() throws ProbDistParmException {
			return this.expVariates();
		}

		@Override
		protected double getDensity(Double1D pt) throws ProbDistParmException {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PostYDist extends ProbDistInitializeByChain<Double1D> {
		private Double1D postM;
		private Double0D postS;
		private TruncatedNormalDist truncatedNormalDist;

		public PostYDist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double1D getPostM() {
			return this.postM;
		}

		protected Double0D getPostS() {
			return this.postS;
		}

		protected void setPostM(Double1D postM) {
			this.postM = postM;
		}

		protected void setPostS(Double0D postS) {
			this.postS = postS;
		}

		protected Double0D getLambda() {
			return (Double0D) this.chainParms[0];
		}

		protected Double0D getMu() {
			return (Double0D) this.chainParms[1];
		}

		protected Double0D getSigma() {
			return (Double0D) this.chainParms[2];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(3);
			this.chainParmNames.add("lambda");
			this.chainParmNames.add("mu");
			this.chainParmNames.add("sigma");
			this.fixedParmNames = new ArrayList<String>(0);
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(3);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					3);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
					0);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPostS(this.getSigma());
			Double1D x = this.getSamplerData().getCol(0);
			Double1D m = x.minus(this.getMu().value()).minus(
					this.getLambda().mult(this.getSigma().pow(2)).value());
			this.setPostM(m);
		}

		private Double1D normalVariates() throws ProbDistParmException {
			double[] v = new double[this.getPostM().size()];
			for (int i=0; i < this.getPostM().size(); i++) { 
				if (this.truncatedNormalDist == null) {
					this.truncatedNormalDist = new TruncatedNormalDist(this.getPostM().get(i), this
							.getPostS(), new Double0D(0));
				}
				v[i] = this.truncatedNormalDist.variateFast().value();
			}
			return new Double1D(v);
		}

		@Override
		protected Double1D genVariate() throws ProbDistParmException {
			return this.normalVariates();
		}

		@Override
		protected double getDensity(Double1D pt) throws ProbDistParmException {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PostLambdaDist extends ProbDistInitializeByChain<Double0D> {
		private Double0D postK;
		private Double0D postTheta;
		private GammaDist gammaDist;

		public PostLambdaDist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double0D getPostK() {
			return this.postK;
		}

		protected Double0D getPostTheta() {
			return this.postTheta;
		}

		protected void setPostK(Double0D postK) {
			this.postK = postK;
		}

		protected void setPostTheta(Double0D postTheta) {
			this.postTheta = postTheta;
		}

		protected Double0D getK() {
			return (Double0D) this.fixedParms[0];
		}

		protected Double0D getTheta() {
			return (Double0D) this.fixedParms[1];
		}

		protected Double1D getY() {
			return (Double1D) this.chainParms[0];
		}

		private Double0D gammaVariate() throws ProbDistParmException {
			if (this.gammaDist == null) {
				this.gammaDist = new GammaDist(this.getPostK(), this
						.getPostTheta());
				return this.gammaDist.variate();
			} else {
				return this.gammaDist.variate(this.getPostK(), this
						.getPostTheta());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(1);
			this.chainParmNames.add("y");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("k_l");
			this.fixedParmNames.add("th_l");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					1);
			this.chainParmClasses.add(Double1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
					2);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPostK(this.getK().plus(this.getSamplerData().size()));
			this.setPostTheta(this.getTheta().plus(this.getY().sum()));
		}

		@Override
		protected Double0D genVariate() throws ProbDistParmException {
			return this.gammaVariate();
		}

		@Override
		protected double getDensity(Double0D pt) throws ProbDistParmException {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PriorMuDist extends ProbDistInitializeByChain<Double0D> {
		private NormalDist normalDist;

		public PriorMuDist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double0D getM() {
			return (Double0D) this.fixedParms[0];
		}

		protected Double0D getNu() {
			return (Double0D) this.fixedParms[1];
		}

		protected Double0D getSigma() {
			return (Double0D) this.chainParms[0];
		}

		private Double0D normalVariate() throws ProbDistParmException {
			if (this.normalDist == null) {
				this.normalDist = new NormalDist(this.getM(), this.getSigma()
						.divide(this.getNu().sqrt()));
				return this.normalDist.variate();
			} else {
				return this.normalDist.variate(this.getM(), this.getSigma()
						.divide(this.getNu().sqrt()));
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(2);
			this.chainParmNames.add("y");
			this.chainParmNames.add("sigma");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("k_l");
			this.fixedParmNames.add("th_l");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					1);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
					2);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		@Override
		protected Double0D genVariate() throws ProbDistParmException {
			return this.normalVariate();
		}

		@Override
		protected double getDensity(Double0D pt) throws ProbDistParmException {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PostMuDist extends ProbDistInitializeByChain<Double0D> {
		private Double0D postM;
		private Double0D postS;
		private NormalDist normalDist;

		public PostMuDist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double0D getPostM() {
			return this.postM;
		}

		protected Double0D getPostS() {
			return this.postS;
		}

		protected void setPostM(Double0D postM) {
			this.postM = postM;
		}

		protected void setPostS(Double0D postS) {
			this.postS = postS;
		}

		protected Double0D getM() {
			return (Double0D) this.fixedParms[0];
		}

		protected Double0D getNu() {
			return (Double0D) this.fixedParms[1];
		}

		protected Double1D getY() {
			return (Double1D) this.chainParms[0];
		}

		protected Double0D getSigma() {
			return (Double0D) this.chainParms[1];
		}

		private Double0D normalVariate() throws ProbDistParmException {
			if (this.normalDist == null) {
				this.normalDist = new NormalDist(this.getPostM(), this
						.getPostS());
				return this.normalDist.variate();
			} else {
				return this.normalDist
						.variate(this.getPostM(), this.getPostS());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(2);
			this.chainParmNames.add("y");
			this.chainParmNames.add("sigma");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("k_l");
			this.fixedParmNames.add("th_l");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					1);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
					2);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			// p' = s^-2 * (N + n)
			// m' = (sum(xi - yi) + nm)/(N + n)
			int N = this.getSamplerData().size();
			Double0D p = this.getSigma().recip().pow(2).mult(this.getNu().plus(N));
			Double0D v = this.getNu().plus(N).recip();
			Double0D m = this.getSamplerData().getCol(0).minus(
					this.getY()).sum().plus(this.getM().mult(this.getNu())).mult(v);
			this.setPostS(p.sqrt().recip());
			this.setPostM(m);
		}

		@Override
		protected Double0D genVariate() throws ProbDistParmException {
			return this.normalVariate();
		}

		@Override
		protected double getDensity(Double0D pt) throws ProbDistParmException {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PostSigmaDist extends ProbDistInitializeByChain<Double0D> {
		private Double0D postK;
		private Double0D postTheta;
		private GammaDist gammaDist;

		public PostSigmaDist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double0D getPostK() {
			return this.postK;
		}

		protected Double0D getPostTheta() {
			return this.postTheta;
		}

		protected void setPostK(Double0D postK) {
			this.postK = postK;
		}

		protected void setPostTheta(Double0D postTheta) {
			this.postTheta = postTheta;
		}

		protected Double0D getK() {
			return (Double0D) this.fixedParms[0];
		}

		protected Double0D getTheta() {
			return (Double0D) this.fixedParms[1];
		}

		protected Double0D getM() {
			return (Double0D) this.fixedParms[2];
		}

		protected Double0D getNu() {
			return (Double0D) this.fixedParms[3];
		}

		protected Double1D getY() {
			return (Double1D) this.chainParms[0];
		}

		private Double0D gammaVariate() throws ProbDistParmException {
			if (this.gammaDist == null) {
				this.gammaDist = new GammaDist(this.getPostK(), this
						.getPostTheta());
				return this.gammaDist.variate();
			} else {
				return this.gammaDist.variate(this.getPostK(), this
						.getPostTheta());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(1);
			this.chainParmNames.add("y");
			this.chainParmNames.add("mu");
			this.fixedParmNames = new ArrayList<String>(4);
			this.fixedParmNames.add("k_s");
			this.fixedParmNames.add("th_s");
			this.fixedParmNames.add("m_m");
			this.fixedParmNames.add("n_m");

			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					2);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
					4);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			int N = this.getSamplerData().size();
			Double0D k = this.getK().plus((N + 1) / 2);
			this.setPostK(k);
			Double1D xmy = this.getSamplerData().getCol(0).minus(this.getY());
			Double0D sxmy2 = xmy.outer(xmy).diag().sum().sum();
			Double0D muTerms = sxmy2.plus(this.getNu().mult(this.getM().pow(2)));
			Double0D ctsTerm = xmy.sum().plus(this.getNu().mult(this.getM())).pow(2).divide(this.getNu().plus(N));
			Double0D th = this.getTheta().plus(muTerms.minus(ctsTerm).mult(0.5));
			this.setPostTheta(th);
		}

		@Override
		protected Double0D genVariate() throws ProbDistParmException {
			return this.gammaVariate().recip().sqrt();
		}

		@Override
		protected double getDensity(Double0D pt) throws ProbDistParmException {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	public ExGaussian() {
		super();
	}

	public ExGaussian(Map<String, Numeric> hypers,
			Map<String, Numeric> init)
			throws ProbDistParmException {
		super(hypers, 1);
		// Set up parameters
		this.params = new HashMap<String, RandomVar<? extends Numeric>>();
		// lambda
		ProbDist<Double0D> priorLambda = new GammaDist((Double0D) this
				.getHyper("k_l"), (Double0D) this.getHyper("th_l"));
		PostLambdaDist postLambda = new PostLambdaDist((Double0D) this
				.getHyper("k_l"), (Double0D) this.getHyper("th_l"));
		Double0D defaultLambda = (Double0D) init.get("lambda");
		RandomVar<Double0D> rvLambda = new RandomVar<Double0D>("lambda",
				priorLambda, postLambda, defaultLambda);
		this.params.put("lambda", rvLambda);
		// mu
		PriorMuDist priorMu = new PriorMuDist((Double0D) this.getHyper("m_m"),
				(Double0D) this.getHyper("n_m"));
		PostMuDist postMu = new PostMuDist((Double0D) this.getHyper("m_m"),
				(Double0D) this.getHyper("n_m"));
		Double0D defaultMu = (Double0D) init.get("mu");
		RandomVar<Double0D> rvMu = new RandomVar<Double0D>("mu", priorMu,
				postMu, defaultMu);
		this.params.put("mu", rvMu);
		// y
		PriorYDist priorY = new PriorYDist();
		PostYDist postY = new PostYDist();
		Double1D defaultY = (Double1D) init.get("y");
		RandomVar<Double1D> rvY = new RandomVar<Double1D>("y", priorY, postY,
				defaultY);
		this.params.put("y", rvY);
		// sigma
		ProbDist<Double0D> priorSigma = new GammaDist((Double0D) this
				.getHyper("k_s"), (Double0D) this.getHyper("th_s"));
		PostSigmaDist postSigma = new PostSigmaDist((Double0D) this
				.getHyper("k_s"), (Double0D) this.getHyper("th_s"),
				(Double0D) this.getHyper("m_m"), (Double0D) this
						.getHyper("n_m"));
		Double0D defaultSigma = (Double0D) init.get("sigma");
		RandomVar<Double0D> rvSigma = new RandomVar<Double0D>("sigma",
				priorSigma, postSigma, defaultSigma);
		this.params.put("sigma", rvSigma);
	}

	@Override
	public ChainLink getInitialLink() {
		return new ChainLink(this.getParam("lambda"), this.getParam("mu"), this
				.getParam("sigma"), this.getParam("y"));
	}

}
