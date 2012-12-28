package org.jqgibbs.models;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.jqgibbs.ChainLink;
import org.jqgibbs.Model;
import org.jqgibbs.SamplerData;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Double3D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Integer1D;
import org.jqgibbs.mathstat.ListSequence;
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.mathstat.RandomVar;
import org.jqgibbs.mathstat.probdist.IIDDist;
import org.jqgibbs.mathstat.probdist.InverseWishartDist;
import org.jqgibbs.mathstat.probdist.MVNormalDist;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistMC;
import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;
import org.jqgibbs.mathstat.probdist.SeqInverseWishartDist;
import org.jqgibbs.mathstat.probdist.SeqMVNormalDist;

public class ToyModel extends Model {

	public int getDims() {
		return this.dims;
	}

	public ToyModel(Map<String, Numeric<? extends Numeric<?>>> hypers,
			Map<String, Numeric<? extends Numeric<?>>> init, int dims)
			throws ProbDistParmException {
		super(hypers, dims);
		// Set up parameters
		this.params = new HashMap<String, RandomVar<? extends Numeric<?>>>();
		// Mu
		ProbDist<Double1D> priorMu = new MVNormalDist((Double1D) this
				.getHyper("m0"), (Double2D) this.getHyper("s0"));
		ProbDistMC<Double1D> postMu = new ProbDistMC<Double1D>(
				(Double1D) this.getHyper("m0"), (Double2D) this.getHyper("s0")) {
			private MVNormalDist mvnDist;
			private Double1D sg0InvM0;
			private Double2D sg0Inv;
			private Double1D postMu;
			private Double2D postSg;

			@Override
			protected void installParmChecks() {
				// Names
				this.chainParmNames = new ArrayList<String>(1);
				this.chainParmNames.add("sg");
				this.fixedParmNames = new ArrayList<String>(2);
				this.fixedParmNames.add("m0");
				this.fixedParmNames.add("s0");
				// Checks
				this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
				this.chainParmCheck.add(null);
				this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
				this.fixedParmCheck.add(null);
				this.fixedParmCheck.add(null);
				// Classes
				this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(1);
				this.chainParmClasses.add(Double2D.class);
				this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(2);
				this.fixedParmClasses.add(Double1D.class);
				this.fixedParmClasses.add(Double2D.class);
			}

			private Double1D mvnVariate() throws ProbDistParmException {
				if (this.mvnDist == null) {
					this.mvnDist = new MVNormalDist(this.getPostMu(),
							this.getPostSg());
					return this.mvnDist.variate();
				} else {
					return this.mvnDist.variate(this.getPostMu(), this.getPostSg());
				}
			}			
			
			protected Double2D getSg() {
				return (Double2D) this.chainParms[0];
			}

			protected Double1D getM0() {
				return (Double1D) this.fixedParms[0];
			}

			protected Double2D getS0() {
				return (Double2D) this.fixedParms[1];
			}

			protected Double1D getS0InvM0() {
				if (this.sg0InvM0 == null) {
					this.sg0InvM0 = this.getS0Inv().mult(this.getM0());
				}
				return this.sg0InvM0;
			}

			protected Double2D getS0Inv() {
				if (this.sg0Inv == null) {
					this.sg0Inv = this.getS0().inverse();
				}
				return this.sg0Inv;
			}

			@Override
			protected void setUpFromChainParms() {
				Double1D sumX = this.getSamplerData().sum();
				int n = this.getSamplerData().size();
				Double2D currSgInv = this.getSg().inverse();
				this.setPostSg(this.getS0Inv().plus(currSgInv.mult(n)).inverse());
				this.setPostMu(this.getPostSg().mult(this.getS0InvM0().plus(
							currSgInv.mult(sumX))));
			}

			private Double1D getPostMu() {
				return this.postMu;
			}
			
			private void setPostMu(Double1D mu) {
				this.postMu = mu;
			}

			private Double2D getPostSg() {
				return this.postSg;
			}

			private void setPostSg(Double2D sg) {
				this.postSg = sg;
			}

			@Override
			protected Double1D genVariate() {
				try {
					// FIXME - create
					return this.mvnVariate();
				} catch (ProbDistParmException e) {
					// FIXME
					return null;
				}
			}

			@Override
			protected double getDensity(Double1D pt) {
				throw new UnsupportedOperationException(
						"Too lazy, come back later");
			}
		};
		Double1D defaultMu = (Double1D) init.get("mu"); // FIXME
		RandomVar<Double1D> rvMu = new RandomVar<Double1D>("mu", priorMu,
				postMu, defaultMu);
		this.params.put("mu", rvMu);
		// Sg
		ProbDist<Double2D> priorSg = new InverseWishartDist((Double2D) this.getHyper("p0"),
						(Integer0D) this.getHyper("k0"));
		ProbDistMC<Double2D> postSg = new ProbDistMC<Double2D>(
				(Double2D) this.getHyper("p0"),
				(Integer0D) this.getHyper("k0")) {
			private InverseWishartDist iwDist;
			private Double2D postPsi;
			private Integer0D postK;

			private Double2D iwVariate() throws ProbDistParmException {
				if (this.iwDist == null) {
					this.iwDist = new InverseWishartDist(this.getPostPsi(),
							this.getPostK());
					return this.iwDist.variate();
				} else {
					return this.iwDist.variate(this.getPostPsi(), this.getPostK());
				}
			}

			private Integer0D getPostK() {
				return this.postK;
			}

			private Double2D getPostPsi() {
				return this.postPsi;
			}
			
			@Override
			protected void installParmChecks() {
				// Names
				this.chainParmNames = new ArrayList<String>(1);
				this.chainParmNames.add("mu");
				this.fixedParmNames = new ArrayList<String>(2);
				this.fixedParmNames.add("p0");
				this.fixedParmNames.add("k0");
				// Checks
				this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(1);
				this.chainParmCheck.add(null);
				this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
				this.fixedParmCheck.add(null);
				this.fixedParmCheck.add(null);
				// Classes
				this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
						1);
				this.chainParmClasses.add(Double1D.class);
				this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
						2);
				this.fixedParmClasses.add(Double2D.class);
				this.fixedParmClasses.add(Integer0D.class);
			}

			protected Double1D getMu() {
				return (Double1D) this.chainParms[0];
			}

			protected Double2D getP0() {
				return (Double2D) this.fixedParms[0];
			}

			protected Integer0D getK0() {
				return (Integer0D) this.fixedParms[1];
			}

			@Override
			protected void setUpFromChainParms() {
				int n = this.getSamplerData().size();
				this.setPostK(new Integer0D(this.getK0().value() + n));
				Double2D s = this.getSamplerData().getNumericValue().minus(this.getMu()).sumOuter();
				this.setPostPsi(this.getP0().plus(s));					
			}

			private void setPostPsi(Double2D psi) {
				this.postPsi = psi;
			}

			private void setPostK(Integer0D K) {
				this.postK = K;
			}

			@Override
			protected Double2D genVariate() throws ProbDistParmException {
				return this.iwVariate();
			}

			@Override
			protected double getDensity(Double2D pt)
					throws ProbDistParmException {
				throw new UnsupportedOperationException(
						"Too lazy, come back later");
			}
		};
		Double2D defaultSg = (Double2D) init.get("sg");
		RandomVar<Double2D> rvSg = new RandomVar<Double2D>("sg", priorSg,
				postSg, defaultSg);
		this.params.put("sg", rvSg);

	}

	@Override
	public ChainLink getInitialLink() {
		return new ChainLink(this.getParam("mu"), this.getParam("sg"));
	}
}
