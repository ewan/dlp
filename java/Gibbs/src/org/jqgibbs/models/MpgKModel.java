package org.jqgibbs.models;

// FIXME - incorporate per-cat nu into Mpg model

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.jqgibbs.ChainLink;
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
import org.jqgibbs.mathstat.probdist.CategoricalDist;
import org.jqgibbs.mathstat.probdist.GammaDist;
import org.jqgibbs.mathstat.probdist.InverseWishartDist;
import org.jqgibbs.mathstat.probdist.MVNormalDist;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistInitializeByChain;
import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;

public class MpgKModel extends MpgFModel {
	
	class PostMuHDist extends PostMuDist {
		public PostMuHDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}	

		protected Double1D getNuBs() {
			return (Double1D) this.chainParms[5];
		}

		@Override
		protected void installParmChecks() {
			super.installParmChecks();
			this.chainParmNames.set(5, "nubs");
			this.chainParmClasses.set(5, Double1D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPostMus(new ListSequence<Double1D>());
			this.setPostSgs(new ListSequence<Double2D>());
			Integer1D active = this.getZ().items();
			Integer0D one = new Integer0D(1);
			for (Integer0D k : active) {
				Integer1D zK = this.getZ().which(k);
				Integer1D zKB1 = zK.intersect(this.getB().which(one));

				Double2D xK = this.getSamplerData().getAll(zK.value());
				Double2D xKB1 = this.getSamplerData().getAll(
						zKB1.value());

				Integer0D nK = new Integer0D(xK.size());
				Integer0D nKB1 = new Integer0D(xKB1.size());

				Double0D postNuB = this.getNuBs().get(k.value()).plus(nKB1);
				Double0D postNuM = this.getNuM().plus(nK).minus(nKB1.pow(2).divide(postNuB));
				this.getPostSgs().add(
						this.getSgs().get(k.value()).mult(postNuM.recip()));

				Double1D mm = xK.sum().plus(this.getMu0().mult(this.getNuM()));
				Double1D mb = xKB1.sum().plus(
						this.getBe0().mult(this.getNuBs().get(k.value()))).mult(nKB1.divide(postNuB));
				this.getPostMus().add(mm.minus(mb).mult(postNuM.recip()));
			}
		}
	}	

	class PostBeHDist extends PostBeDist {
		public PostBeHDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}
		
		protected Double1D getNuBs() {
			return (Double1D) this.chainParms[4];
		}

		@Override
		protected void installParmChecks() {
			super.installParmChecks();
			this.chainParmNames.set(4, "nubs");
			this.chainParmClasses.set(4, Double1D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPostMus(new ListSequence<Double1D>());
			this.setPostSgs(new ListSequence<Double2D>());
			Integer1D active = this.getZ().items();
			Integer0D one = new Integer0D(1);
			for (Integer0D k : active) {
				Integer1D zK = this.getZ().which(k);
				Integer1D zKB1 = zK.intersect(this.getB().which(one));

				Double1D mu = this.getMus().get(k.value());
				Double2D xKB1 = this.getSamplerData().getAll(zKB1.value()).minus(mu);
				Integer0D nKB1 = new Integer0D(xKB1.size());
				Double1D sKB1 = xKB1.sum();

				Double0D postNuB = this.getNuBs().get(k.value()).plus(nKB1);
				this.getPostSgs().add(this.getSgs().get(k.value()).divide(postNuB));
				//FIXME
				this.getPostMus().add(sKB1.plus(this.getBe0().mult(this.getNuBs().get(k.value()))).mult(postNuB.recip()));
			}
		}

	}
	
	class PostSgHDist extends PostSgDist {
		public PostSgHDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double1D getNuBs() {
			return (Double1D) this.chainParms[2];
		}

		@Override
		protected void installParmChecks() {
			super.installParmChecks();
			this.chainParmNames.set(2, "nubs");
			this.chainParmClasses.set(2, Double1D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			Double0D nuM = this.getNuM();

			this.setPostPs(new ListSequence<Double2D>());
			this.setPostKs(new ListSequence<Integer0D>());
			Integer1D active = this.getZ().items();
			Integer0D one = new Integer0D(1);
			for (Integer0D k : active) {
				Integer1D zK = this.getZ().which(k);
				Integer1D zKB1 = zK.intersect(this.getB().which(one));

				Double2D xK = this.getSamplerData().getAll(zK.value());
				Double2D xKB1 = this.getSamplerData().getAll(
						zKB1.value());

				Integer0D nK = new Integer0D(xK.size());
				Integer0D nKB1 = new Integer0D(xKB1.size());

				Double1D sK = xK.sum();
				Double1D sKB1 = xKB1.sum();

				Double0D postNuB = this.getNuBs().get(k.value()).plus(nKB1);
				Double0D postNuM = nuM.plus(nK).minus(nKB1.pow(2).divide(postNuB));

				this.getPostKs().add(this.getK0().plus(nK));

				Double2D s = xK.sumOuter();
				Double2D m0 = this.getMu0().outer(this.getMu0()).mult(nuM);
				Double2D b0 = this.getBe0().outer(this.getBe0()).mult(this.getNuBs().get(k.value()));
				Double1D m = sK.plus(this.getMu0().mult(nuM));
				Double2D mm = m.outer(m).divide(nuM.plus(nK));
				Double1D n = m.mult(nKB1.divide(nuM.plus(nK))).minus(sKB1).plus(
						this.getBe0().mult(this.getNuBs().get(k.value())));
				Double2D nn = n.outer(n).mult(nuM.plus(nK).divide(postNuM.mult(postNuB)));

				Double2D H = this.getP0().plus(s).plus(m0).plus(b0).minus(mm).minus(nn);
				this.getPostPs().add(H);
			}
		}
	}	
	
	class CPHDist extends CPDist {
		public CPHDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}
		
		protected ProbDist<Double0D> baseNuB;

		public void setBaseNuB(ProbDist<Double0D> baseNuB) {
			this.baseNuB = baseNuB;
		}

		public ProbDist<Double0D> getBaseNuB() {
			return this.baseNuB;
		}

		protected Double1D getNuBs() {
			return (Double1D) this.chainParms[8];
		}
		
		protected Double0D getNuBA() {
			return (Double0D) this.fixedParms[2];
		}
		
		protected Double0D getNuBB() {
			return (Double0D) this.fixedParms[3];
		}

		@Override
		protected void installParmChecks() {
			super.installParmChecks();
			this.fixedParmNames = new ArrayList<String>(4);
			this.fixedParmNames.add("p0");
			this.fixedParmNames.add("k0");
			this.fixedParmNames.add("nuba");
			this.fixedParmNames.add("nubb");
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.fixedParmCheck.add(null);
			this.chainParmNames.set(8, "nubs");
			this.chainParmClasses.set(8, Double1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(4);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Integer0D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected Integer1D genVariate() throws ProbDistParmException {
			// FIXME - This sampling scheme destructively modifies
			// mu, sg, be - HACK
			// Now also with nub!
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
				// Sample a new category if necessary, remove the current
				// category
				// from consideration
				Integer1D active = Z.items();
				int maxActive = active.get(active.size() - 1).value();
				Integer0D zedI = Z.get(i);
				Z.set(i, unused);
				Integer1D zZedI = Z.which(zedI);
				int nZZedI = zZedI.size();
				if (nZZedI > 0) { // Sample new category
					Integer0D zedNew = Z.minNotIn(0, maxActive);
					this.resetPostMvn(zedNew);
					Double2D sgNew = this.getBaseSg().variate(this.getP0(),
							this.getK0());
					Double2D sgMu = sgNew.mult(1 / this.getNuM().value());
					Double1D muNew = this.getBaseMu().variate(
							this.getMu0(), sgMu);
					Double0D nuBNew = this.getBaseNuB().variate(this.getNuBA(), this.getNuBB());
					Double2D sgBe = sgNew.mult(1 / nuBNew.value());
					Double1D beNew = this.getBaseBe().variate(this.getBe0(), sgBe);
					// Add it
					this.getSgs().set(zedNew.value(), sgNew);
					this.getMus().set(zedNew.value(), muNew);
					this.getBes().set(zedNew.value(), beNew);
					this.getNuBs().set(zedNew.value(), nuBNew);
					active.add(zedNew);
				}
				// Select a category for this point
				double[] logP = new double[active.size()];
				double maxLogP = -Double.MAX_VALUE; // FIXME??
				int ai = 0;
				for (Integer0D k : active) {
					Double1D xIA = this.getSamplerData().get(i);
					Integer1D zK = Z.which(k);
					int nZK = zK.size();
					this.setUpPostMvn(k);
					Double1D be = this.getBes().get(k.value());
					if (this.getB().get(i).value() == 1) {
						// Get current data point						
						xIA = xIA.minus(be);
					}
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
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}	
	
	class PostBe0HDist extends PostBe0Dist {
		public PostBe0HDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}	

		@Override
		protected void installParmChecks() {
			super.installParmChecks();
			this.chainParmNames.set(2, "nubs");
			this.chainParmClasses.set(2, Double1D.class);
		}
		
		protected Double1D getNuBs() {
			return (Double1D) this.chainParms[2];
		}

		@Override
		protected void setUpFromChainParms() {
			Double2D si = this.getS0Inv();
			Double1D sim = this.getS0InvB0();
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				Double2D pK = this.getSgs().get(k.value()).inverse();
				si = si.plus(pK.mult(this.getNuBs().get(k.value())));
				sim = sim.plus(pK.mult(this.getNuBs().get(k.value())).mult(
						this.getBes().get(k.value()))); // FIXME??
			}
			Double2D s = si.inverse();
			this.setPostSg(s);
			this.setPostMu(s.mult(sim));
		}
	}
	
	class PriorNuBsDist extends ProbDistInitializeByChain<Double1D> {
		private HashMap<Integer0D,GammaDist> gammaDists; // FIXME - again, this is the wrong way to do this
		
		public PriorNuBsDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}		
		
		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(5);
			this.chainParmNames.add("sg");
			this.chainParmNames.add("be");
			this.chainParmNames.add("be0");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("numa");
			this.fixedParmNames.add("numb");
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
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		private Double1D gammaVariates() throws ProbDistParmException {
			Double1D gs = new Double1D(new double[this.getBes().size()]); // FIXME - whoops
			for (Integer0D k : this.gammaDists.keySet()) {
				gs.set(k.value(), this.gammaDists.get(k).variate());
			}
			return gs;
		}

		protected Double3D getSgs() {
			return (Double3D) this.chainParms[0];
		}

		protected Double2D getBes() {
			return (Double2D) this.chainParms[1];
		}

		protected Double1D getBe0() {
			return (Double1D) this.chainParms[2];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[3];
		}

		protected Integer1D getB() {
			return (Integer1D) this.chainParms[4];
		}

		protected Double0D getNuBA() {
			return (Double0D) this.fixedParms[0];
		}

		protected Double0D getNuBB() {
			return (Double0D) this.fixedParms[1];
		}

		@Override
		protected void setUpFromChainParms() {
			this.gammaDists = new HashMap<Integer0D,GammaDist>();
			
			Double0D a0 = this.getNuBA();
			Double0D b0 = this.getNuBB();
			
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				try {
					this.gammaDists.put(k, new GammaDist(a0, b0));
				} catch (ProbDistParmException e) {
					// TODO Auto-generated catch block
					e.printStackTrace(); // FIXME - god...
				}
			}
		}
		@Override
		protected Double1D genVariate() throws ProbDistParmException {
			return this.gammaVariates();
		}

		@Override
		protected double getDensity(Double1D pt) {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}
	
	class PostNuBsDist extends ProbDistInitializeByChain<Double1D> {
		private HashMap<Integer0D,GammaDist> gammaDists; // FIXME - again, this is the wrong way to do this
		private Double0D postNuBA;
		private Double0D postNuBB;

		public PostNuBsDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}
		
		private void setUpPostGamma(Integer0D k) throws ProbDistParmException {
			if (this.gammaDists == null) {
				this.gammaDists = new HashMap<Integer0D,GammaDist>();
			}
			if (!this.gammaDists.containsKey(k)) {
				this.gammaDists.put(k, new GammaDist(this.getPostNuBA(), this.getPostNuBB()));
			}
			// else do nothing
		}		
		
		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(5);
			this.chainParmNames.add("sg");
			this.chainParmNames.add("be");
			this.chainParmNames.add("be0");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("numa");
			this.fixedParmNames.add("numb");
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
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		private Double1D gammaVariates() throws ProbDistParmException {
			Double1D gs = new Double1D(new double[this.getBes().size()]); // FIXME - whoops
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				gs.set(k.value(), this.gammaDists.get(k).variate());
			}
			return gs;
		}

		protected Double3D getSgs() {
			return (Double3D) this.chainParms[0];
		}

		protected Double2D getBes() {
			return (Double2D) this.chainParms[1];
		}

		protected Double1D getBe0() {
			return (Double1D) this.chainParms[2];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[3];
		}

		protected Integer1D getB() {
			return (Integer1D) this.chainParms[4];
		}

		protected Double0D getNuBA() {
			return (Double0D) this.fixedParms[0];
		}

		protected Double0D getNuBB() {
			return (Double0D) this.fixedParms[1];
		}

		@Override
		protected void setUpFromChainParms() {
			Double0D a0 = this.getNuBA();
			Double0D b0 = this.getNuBB();

			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				Double1D be = this.getBes().get(k.value()).minus(this.getBe0());
				Double0D a = a0.plus(0.5*((double) MpgKModel.this.getDims()));
				Double0D b = b0.plus(0.5*be.mult(this.getSgs().get(k.value()).inverse()).mult(be));
				this.setPostNuBA(a);
				this.setPostNuBB(b);
				try {
					this.setUpPostGamma(k);
				} catch (ProbDistParmException e) {
					// TODO Auto-generated catch block
					e.printStackTrace(); // FIXME FIXME FIXME
				}
			}
		}

		private Double0D getPostNuBA() {
			return this.postNuBA;
		}

		private void setPostNuBA(Double0D nuBA) {
			this.postNuBA = nuBA;
		}

		private Double0D getPostNuBB() {
			return this.postNuBB;
		}

		private void setPostNuBB(Double0D nuBB) {
			this.postNuBB = nuBB;
		}

		@Override
		protected Double1D genVariate() throws ProbDistParmException {
			return this.gammaVariates();
		}

		@Override
		protected double getDensity(Double1D pt) {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}	
	
	public MpgKModel(Map<String, Numeric<? extends Numeric<?>>> hypers,
			Map<String, Numeric<? extends Numeric<?>>> init, int dims)
			throws ProbDistParmException {
		super(hypers, init, dims);

		// Mu
		PriorMuDist priorMu = new PriorMuDist();
		PostMuHDist postMu = new PostMuHDist();
		Double2D defaultMu = (Double2D) init.get("mu"); // FIXME
		RandomVar<Double2D> rvMu = new RandomVar<Double2D>("mu", priorMu,
				postMu, defaultMu);
		this.params.put("mu", rvMu);
		
		// Be
		PriorBeDist priorBe = new PriorBeDist();
		PostBeHDist postBe = new PostBeHDist();
		Double2D defaultBe = (Double2D) init.get("be"); // FIXME
		RandomVar<Double2D> rvBe = new RandomVar<Double2D>("be", priorBe,
				postBe, defaultBe);
		this.params.put("be", rvBe);		
		
		// Sg
		PriorSgDist priorSg = new PriorSgDist((Double2D) this.getHyper("p0"), (Integer0D) this.getHyper("k0"));
		PostSgHDist postSg = new PostSgHDist((Double2D) this.getHyper("p0"), (Integer0D) this.getHyper("k0"));
		Double3D defaultSg = (Double3D) init.get("sg"); // FIXME
		RandomVar<Double3D> rvSg = new RandomVar<Double3D>("sg", priorSg,
				postSg, defaultSg);
		this.params.put("sg", rvSg);
		
		// Z
		ProbDist<Integer1D> priorZ = new Deterministic<Integer1D>(init.get("Z"));
		CPHDist postZ = new CPHDist((Double2D) this.getHyper("p0"),
				(Integer0D) this.getHyper("k0"), (Double0D) this.getHyper("nuba"),
				(Double0D) this.getHyper("nubb")); // FIXME -- Hack!
		postZ.setBaseSg(new InverseWishartDist()); // FIXME -- Eek! Hack!
		postZ.setBaseNuB(new GammaDist());
		postZ.setBaseMu(new MVNormalDist()); // FIXME -- Eek! Hack!
		postZ.setBaseBe(new MVNormalDist()); // FIXME -- Eek! Hack!
		Integer1D defaultZ = (Integer1D) init.get("Z");
		RandomVar<Integer1D> rvZ = new RandomVar<Integer1D>("Z", priorZ, postZ,
				defaultZ);
		this.params.put("Z", rvZ); // FIXME
		
		
		// be0
		ProbDist<Double1D> priorBe0 = new MVNormalDist((Double1D) this
				.getHyper("b0"), (Double2D) this.getHyper("sb0"));
		Double1D defaultBe0 = (Double1D) init.get("be0"); // FIXME
		PostBe0HDist postBe0 = new PostBe0HDist((Double1D) this.getHyper("b0"), (Double2D) this.getHyper("sb0"));
		RandomVar<Double1D> rvBe0 = new RandomVar<Double1D>("be0", priorBe0,
				postBe0, defaultBe0);
		this.params.put("be0", rvBe0); // FIXME	
		
		// nubs
		PriorNuBsDist priorNuB = new PriorNuBsDist((Double0D) this
				.getHyper("nuba"), (Double0D) this.getHyper("nubb"));
		Double1D defaultNuB = (Double1D) init.get("nubs"); // FIXME
		PostNuBsDist postNuB = new PostNuBsDist((Double0D) this.getHyper("nuba"), (Double0D) this.getHyper("nubb"));
		RandomVar<Double1D> rvNuB = new RandomVar<Double1D>("nubs", priorNuB,
				postNuB, defaultNuB);
		this.params.put("nubs", rvNuB); // FIXME

	}

	@Override
	public ChainLink getInitialLink() {
		return new ChainLink(this.getParam("Z"), this.getParam("B"),
				this.getParam("be"), this.getParam("mu"), this.getParam("sg"), this
				.getParam("num"), this.getParam("nub"), this.getParam("nubs"), this.getParam("mu0"),
				this.getParam("be0"), this.getParam("x"), this.getParam("al"));
	}
}
