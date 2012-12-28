package org.jqgibbs.models;

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
import org.jqgibbs.mathstat.probdist.BetaDist;
import org.jqgibbs.mathstat.probdist.CategoricalDist;
import org.jqgibbs.mathstat.probdist.GammaDist;
import org.jqgibbs.mathstat.probdist.IIDDist;
import org.jqgibbs.mathstat.probdist.InverseWishartDist;
import org.jqgibbs.mathstat.probdist.MVNormalDist;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistMC;
import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;
import org.jqgibbs.models.MpgModel.PostBDist;
import org.jqgibbs.models.MpgModel.PostNuBDist;

public class MpgJModel extends MpgFModel {
	
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

				SamplerData xK = this.getSamplerData().getAll(zK.value());
				SamplerData xKB1 = this.getSamplerData().getAll(
						zKB1.value());

				Integer0D nK = new Integer0D(xK.size());
				Integer0D nKB1 = new Integer0D(xKB1.size());

				Double1D sK = xK.sum();
				Double1D sKB1 = xKB1.sum();

				Double0D postNuB = this.getNuBs().get(k.value()).plus(nKB1);
				Double0D postNuM = nuM.plus(nK).minus(nKB1.pow(2).divide(postNuB));

				this.getPostKs().add(this.getK0().plus(nK));

				Double2D s = xK.getNumericValue().sumOuter();
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

	class PriorLDist extends ProbDistMC<Integer1D> {
		private IIDDist<Integer1D,Integer0D> catDists;
//		
		public PriorLDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private Integer1D catVariates() throws ProbDistParmException {
			if (this.catDists == null) {
				CategoricalDist catDist = new CategoricalDist(this.getPi());
				this.catDists = new IIDDist<Integer1D,Integer0D>(catDist, this.getBes().size());
			}
//			else {
//				return this.catDists.variate(this.getPostProb());
//			}
			return this.catDists.variate();
		}

		private Double0D getPi() {
			return (Double0D) this.chainParms[0];
		}

		private Integer1D getZ() {
			return (Integer1D) this.chainParms[1];
		}

		private Double2D getBes() {
			return (Double2D) this.chainParms[2];
		}
		
		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(3);
			this.chainParmNames.add("pi");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("be");
			this.fixedParmNames = new ArrayList<String>(0);
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null); // FIXME
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(3);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(0);
		}

		@Override
		protected void setUpFromChainParms() {
			return;
		}

		@Override
		protected Integer1D genVariate() throws ProbDistParmException {
			Integer1D l = new Integer1D(new int[this.getBes().size()]);
			Integer1D samp = this.catVariates();
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				l.set(k.value(), samp.get(k.value()));
			}
			return l;
		}

		@Override
		protected double getDensity(Integer1D pt) {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}
	
	class PostNuBDist extends ProbDistMC<Double0D> {
		private GammaDist gammaDist;
		private Double0D postNuBA;
		private Double0D postNuBB;

		public PostNuBDist(Numeric<?>... fixed) throws ProbDistParmException {
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

		private Double0D gammaVariate() throws ProbDistParmException {
			if (this.gammaDist == null) {
				this.gammaDist = new GammaDist(this.getPostNuBA(), this
						.getPostNuBB());
				return this.gammaDist.variate();
			} else {
				return this.gammaDist.variate(this.getPostNuBA(), this
						.getPostNuBB());
			}
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
			Double0D a = this.getNuBA();
			Double0D b = this.getNuBB();
			Integer1D active = this.getZ().items();
			Integer0D one = new Integer0D(1);
			for (Integer0D k : active) {
//				Integer1D zK = this.getZ().which(k);
//				Integer1D zKB1 = zK.intersect(this.getB().which(one));
//				if (zKB1.size() > 0) {
					Double1D be = this.getBes().get(k.value()).minus(
							this.getBe0());
					a = a.plus(0.5*((double) MpgModel.this.getDims()));
					b = b.plus(0.5*be.mult(this.getSgs().get(k.value()).inverse()).mult(be));
//				}
			}
			this.setPostNuBA(a);
			this.setPostNuBB(b);
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
		protected Double0D genVariate() throws ProbDistParmException {
			return this.gammaVariate();
		}

		@Override
		protected double getDensity(Double0D pt) {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}
	
	class PostPiDist extends ProbDistMC<Double0D> {
		private BetaDist betaDist;
		private PriorLDist lDist;
		private Double0D postPiA;
		private Double0D postPiB;

		public PostPiDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}		
		
		@Override
		protected void installParmChecks() {
			// Names
//			this.chainParmNames = new ArrayList<String>(5);
//			this.chainParmNames.add("sg");
//			this.chainParmNames.add("be");
//			this.chainParmNames.add("be0");
//			this.chainParmNames.add("Z");
//			this.chainParmNames.add("B");
//			this.fixedParmNames = new ArrayList<String>(2);
//			this.fixedParmNames.add("numa");
//			this.fixedParmNames.add("numb");
			// Checks
//			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
//			this.chainParmCheck.add(null);
//			this.chainParmCheck.add(null);
//			this.chainParmCheck.add(null);
//			this.chainParmCheck.add(null);
//			this.chainParmCheck.add(null);
//			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(2);
//			this.fixedParmCheck.add(null);
//			this.fixedParmCheck.add(null);
			// Classes
//			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
//					5);
//			this.chainParmClasses.add(Double3D.class);
//			this.chainParmClasses.add(Double2D.class);
//			this.chainParmClasses.add(Double1D.class);
//			this.chainParmClasses.add(Integer1D.class);
//			this.chainParmClasses.add(Integer1D.class);
//			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
//					2);
//			this.fixedParmClasses.add(Double0D.class);
//			this.fixedParmClasses.add(Double0D.class);
		}

		private Double0D betaVariate() throws ProbDistParmException {
			if (this.betaDist == null) {
				this.betaDist = new BetaDist(this.getPostPiA(), this
						.getPostPiB());
				return this.betaDist.variate();
			} else {
				return this.betaDist.variate(this.getPostPiA(), this
						.getPostPiB());
			}
		}
		
		private Integer1D lVariate() throws ProbDistParmException {
			if (this.lDist == null) {
				this.lDist = new PriorLDist(this.getChainLink(), this.getSamplerData());
				return this.lDist.variate();
			} else {
				return this.lDist.variate(this.getChainLink(), this.getSamplerData());
			}
		}		
//
//		protected Double3D getSgs() {
//			return (Double3D) this.chainParms[0];
//		}
//
//		protected Double2D getBes() {
//			return (Double2D) this.chainParms[1];
//		}
//
//		protected Double1D getBe0() {
//			return (Double1D) this.chainParms[2];
//		}
//
//		protected Integer1D getZ() {
//			return (Integer1D) this.chainParms[3];
//		}
//
//		protected Integer1D getB() {
//			return (Integer1D) this.chainParms[4];
//		}
//
//		protected Double0D getNuBA() {
//			return (Double0D) this.fixedParms[0];
//		}
//
//		protected Double0D getNuBB() {
//			return (Double0D) this.fixedParms[1];
//		}

		@Override
		protected void setUpFromChainParms() {
			Integer1D l;
			try {
				l = this.lVariate();
			} catch (ProbDistParmException e) {
				// TODO Auto-generated catch block
				e.printStackTrace(); // FIXME
			}
//			Double0D a = this.getNuBA();
//			Double0D b = this.getNuBB();
//			Integer1D active = this.getZ().items();
//			Integer0D one = new Integer0D(1);
//			for (Integer0D k : active) {
////				Integer1D zK = this.getZ().which(k);
////				Integer1D zKB1 = zK.intersect(this.getB().which(one));
////				if (zKB1.size() > 0) {
//					Double1D be = this.getBes().get(k.value()).minus(
//							this.getBe0());
//					a = a.plus(0.5*((double) MpgModel.this.getDims()));
//					b = b.plus(0.5*be.mult(this.getSgs().get(k.value()).inverse()).mult(be));
////				}
//			}
//			this.setPostNuBA(a);
//			this.setPostNuBB(b);
		}

		private Double0D getPostPiA() {
			return this.postPiA;
		}

		private void setPostPiA(Double0D piA) {
			this.postPiA = piA;
		}

		private Double0D getPostPiB() {
			return this.postPiB;
		}

		private void setPostPiB(Double0D piB) {
			this.postPiB = piB;
		}

		@Override
		protected Double0D genVariate() throws ProbDistParmException {
			return this.betaVariate();
		}

		@Override
		protected double getDensity(Double0D pt) {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}	
	
	public MpgJModel(Map<String, Numeric<? extends Numeric<?>>> hypers,
			Map<String, Numeric<? extends Numeric<?>>> init, int dims)
			throws ProbDistParmException {
		super(hypers, init, dims);

		// Mu
		PriorMuJDist priorMu = new PriorMuJDist();
		PostMuJDist postMu = new PostMuJDist();
		Double2D defaultMu = (Double2D) init.get("mu"); // FIXME
		RandomVar<Double2D> rvMu = new RandomVar<Double2D>("mu", priorMu,
				postMu, defaultMu);
		this.params.put("mu", rvMu);
		
		// Be
		PriorBeJDist priorBe = new PriorBeJDist();
		PostBeJDist postBe = new PostBeJDist();
		Double2D defaultBe = (Double2D) init.get("be"); // FIXME
		RandomVar<Double2D> rvBe = new RandomVar<Double2D>("be", priorBe,
				postBe, defaultBe);
		this.params.put("be", rvBe);		
		
		// Sg
		PriorSgJDist priorSg = new PriorSgJDist((Double2D) this.getHyper("p0"), (Integer0D) this.getHyper("k0"));
		PostSgJDist postSg = new PostSgJDist((Double2D) this.getHyper("p0"), (Integer0D) this.getHyper("k0"));
		Double3D defaultSg = (Double3D) init.get("sg"); // FIXME
		RandomVar<Double3D> rvSg = new RandomVar<Double3D>("sg", priorSg,
				postSg, defaultSg);
		this.params.put("sg", rvSg);
		
		// Z
		ProbDist<Integer1D> priorZ = new Deterministic<Integer1D>(init.get("Z"));
		CPJDist postZ = new CPJDist((Double2D) this.getHyper("p0"),
				(Integer0D) this.getHyper("k0"), (Double0D) this.getHyper("nuba"),
				(Double0D) this.getHyper("nubb")); // FIXME -- Hack!
		postZ.setBaseL(new CategoricalDist());
		postZ.setBaseSg(new InverseWishartDist()); // FIXME -- Eek! Hack!
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
		PostBe0JDist postBe0 = new PostBe0JDist((Double1D) this.getHyper("b0"), (Double2D) this.getHyper("sb0"));
		RandomVar<Double1D> rvBe0 = new RandomVar<Double1D>("be0", priorBe0,
				postBe0, defaultBe0);
		this.params.put("be0", rvBe0); // FIXME		

		// nub
		ProbDist<Double0D> priorNuB = new GammaDist((Double0D) this
				.getHyper("nuba"), (Double0D) this.getHyper("nubb"));
		Double0D defaultNuB = (Double0D) init.get("nub"); // FIXME
		PostNuBJDist postNuB = new PostNuBJDist((Double0D) this.getHyper("nuba"), (Double0D) this.getHyper("nubb"));
		RandomVar<Double0D> rvNuB = new RandomVar<Double0D>("nub", priorNuB,
				postNuB, defaultNuB);
		this.params.put("nub", rvNuB); // FIXME		
		
		// pi
		ProbDist<Double0D> priorPi = new BetaDist(this.getHyper("lshape1"), this.getHyper("lshape2"));
		Double0D defaultPi = (Double0D) init.get("pi"); // FIXME
		PostPiDist postPi = new PostPiDist(this.getHyper("lshape1"), this.getHyper("lshape2"));		
		RandomVar<Double0D> rvPi = new RandomVar<Double0D>("pi", priorPi, postPi, defaultPi);
	}

	@Override
	public ChainLink getInitialLink() {
		return new ChainLink(this.getParam("Z"), this.getParam("B"),
				this.getParam("be"), this.getParam("mu"), this.getParam("sg"), this
				.getParam("num"), this.getParam("nub"), this.getParam("mu0"),
				this.getParam("be0"), this.getParam("L"), this.getParam("x"),
				this.getParam("al"));
	}
}
