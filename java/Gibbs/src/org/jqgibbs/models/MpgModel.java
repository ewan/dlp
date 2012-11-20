package org.jqgibbs.models;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
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
import org.jqgibbs.mathstat.probdist.BetaDist;
import org.jqgibbs.mathstat.probdist.CategoricalDistInitializeByK;
import org.jqgibbs.mathstat.probdist.CategoricalDistInitializeByP;
import org.jqgibbs.mathstat.probdist.GammaDist;
import org.jqgibbs.mathstat.probdist.IIDDist;
import org.jqgibbs.mathstat.probdist.InverseWishartDist;
import org.jqgibbs.mathstat.probdist.MVNormalDist;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistInitializeByChain;
import org.jqgibbs.mathstat.probdist.ProbDistInitializeDirectly;
import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;
import org.jqgibbs.mathstat.probdist.SeqCategoricalDist;
import org.jqgibbs.mathstat.probdist.SeqInverseWishartDist;
import org.jqgibbs.mathstat.probdist.SeqMVNormalDist;

import com.sun.jdi.Value;

public class MpgModel extends Model {

	public int getDims() {
		return this.dims;
	}
	
	class PriorMuDist extends ProbDistInitializeByChain<Double2D> {
		private ListSequence<Double1D> mu0Seq;
		private ListSequence<Double2D> priorSgs;

		private SeqMVNormalDist mvnDists;
		
		public PriorMuDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private void setPriorSgs(ListSequence<Double2D> priorSgs) {
			this.priorSgs = priorSgs;
		}

		private ListSequence<Double1D> getMu0Seq() {
			if (this.mu0Seq == null) {
				this.mu0Seq = new ListSequence<Double1D>();
				Integer1D active = this.getZ().items();
				for (int i = 0; i < active.size(); i++) {
					this.mu0Seq.add(this.getMu0());
				}
			}
			return this.mu0Seq;
		}

		private ListSequence<Double2D> getPriorSgs() {
			return this.priorSgs;
		}

		private Double3D getSgs() {
			return (Double3D) this.chainParms[0];
		}

		private Double1D getMu0() {
			return (Double1D) this.chainParms[1];
		}

		private Double0D getNuM() {
			return (Double0D) this.chainParms[2];
		}

		private Integer1D getZ() {
			return (Integer1D) this.chainParms[3];
		}

		private Double2D mvnVariates() throws ProbDistParmException {
			if (this.mvnDists == null) {
				this.mvnDists = new SeqMVNormalDist(this.getMu0Seq(), this
						.getPriorSgs());
				return this.mvnDists.variate();
			} else {
				return this.mvnDists.variate(this.getMu0Seq(), this
						.getPriorSgs());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(4);
			this.chainParmNames.add("sg");
			this.chainParmNames.add("mu0");
			this.chainParmNames.add("num");
			this.chainParmNames.add("Z");
			this.fixedParmNames = new ArrayList<String>(0);
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					4);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					0);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPriorSgs(new ListSequence<Double2D>());
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				Double2D sg = this.getSgs().get(k.value()).mult(
						1 / this.getNuM().value());
				this.getPriorSgs().add(sg);
			}
		}

		@Override
		protected Double2D genVariate() throws ProbDistParmException {
			return this.mvnVariates();
		}

		@Override
		protected double getDensity(Double2D pt) {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}	
	
	class PostMuDist extends ProbDistInitializeByChain<Double2D> {
		private ListSequence<Double1D> postMus;
		private ListSequence<Double2D> postSgs;
		private SeqMVNormalDist mvnDists;
		
		public PostMuDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}	

		private Double2D mvnVariates() throws ProbDistParmException {
			if (this.mvnDists == null) {
				this.mvnDists = new SeqMVNormalDist(this.getPostMus(), this
						.getPostSgs());
				return this.mvnDists.variate();
			} else {
				return this.mvnDists.variate(this.getPostMus(), this
						.getPostSgs());
			}
		}

		protected ListSequence<Double1D> getPostMus() {
			return this.postMus;
		}

		protected ListSequence<Double2D> getPostSgs() {
			return this.postSgs;
		}

		protected void setPostMus(ListSequence<Double1D> postMus) {
			this.postMus = postMus;
		}

		protected void setPostSgs(ListSequence<Double2D> postSgs) {
			this.postSgs = postSgs;
		}

		protected Double3D getSgs() {
			return (Double3D) this.chainParms[0];
		}

		protected Double1D getMu0() {
			return (Double1D) this.chainParms[1];
		}

		protected Double0D getNuM() {
			return (Double0D) this.chainParms[2];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[3];
		}

		protected Integer1D getB() {
			return (Integer1D) this.chainParms[4];
		}

		protected Double0D getNuB() {
			return (Double0D) this.chainParms[5];
		}

		protected Double1D getBe0() {
			return (Double1D) this.chainParms[6];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(7);
			this.chainParmNames.add("sg");
			this.chainParmNames.add("mu0");
			this.chainParmNames.add("num");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.chainParmNames.add("nub");
			this.chainParmNames.add("be0");
			this.fixedParmNames = new ArrayList<String>(0);
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(7);
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					7);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(0);
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

				Double0D postNuB = this.getNuB().plus(nKB1);
				Double0D postNuM = this.getNuM().plus(nK).minus(nKB1.pow(2).divide(postNuB));
				this.getPostSgs().add(
						this.getSgs().get(k.value()).mult(postNuM.recip()));

				Double1D mm = xK.sum().plus(this.getMu0().mult(this.getNuM()));
				Double1D mb = xKB1.sum().plus(
						this.getBe0().mult(this.getNuB())).mult(nKB1.divide(postNuB));
				this.getPostMus().add(mm.minus(mb).mult(postNuM.recip()));
			}
		}

		@Override
		protected Double2D genVariate() throws ProbDistParmException {
			double[][] v = new double[this.getSgs().size()][MpgModel.this.getDims()];
			Double2D activeUpdates = this.mvnVariates();
			int i=0;
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], MpgModel.this.getDims());
				i++;
			}
			return new Double2D(v);
		}

		@Override
		protected double getDensity(Double2D pt) {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}
	
	class PriorBeDist extends ProbDistInitializeByChain<Double2D> {
		private ListSequence<Double1D> be0Seq;
		private ListSequence<Double2D> priorSgs;
		
		public PriorBeDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private SeqMVNormalDist mvnDists;

		private void setPriorSgs(ListSequence<Double2D> priorSgs) {
			this.priorSgs = priorSgs;
		}

		private ListSequence<Double1D> getBe0Seq() {
			if (this.be0Seq == null) {
				this.be0Seq = new ListSequence<Double1D>();
				Integer1D active = this.getZ().items();
				for (int i = 0; i < active.size(); i++) {
					this.be0Seq.add(this.getBe0());
				}
			}
			return this.be0Seq;
		}

		private ListSequence<Double2D> getPriorSgs() {
			return this.priorSgs;
		}

		private Double3D getSgs() {
			return (Double3D) this.chainParms[0];
		}

		private Double1D getBe0() {
			return (Double1D) this.chainParms[1];
		}

		private Double0D getNuB() {
			return (Double0D) this.chainParms[2];
		}

		private Integer1D getZ() {
			return (Integer1D) this.chainParms[3];
		}

		private Double2D mvnVariates() throws ProbDistParmException {
			if (this.mvnDists == null) {
				this.mvnDists = new SeqMVNormalDist(this.getBe0Seq(), this
						.getPriorSgs());
				return this.mvnDists.variate();
			} else {
				return this.mvnDists.variate(this.getBe0Seq(), this
						.getPriorSgs());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(4);
			this.chainParmNames.add("sg");
			this.chainParmNames.add("be0");
			this.chainParmNames.add("nub");
			this.chainParmNames.add("Z");
			this.fixedParmNames = new ArrayList<String>(0);
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(4);
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					4);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					0);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPriorSgs(new ListSequence<Double2D>());
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				Double2D sg = this.getSgs().get(k.value()).mult(
						1 / this.getNuB().value());
				this.getPriorSgs().add(sg);
			}
		}

		@Override
		protected Double2D genVariate() throws ProbDistParmException {
			return this.mvnVariates();
		}

		@Override
		protected double getDensity(Double2D pt) {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}

	class PostBeDist extends ProbDistInitializeByChain<Double2D> {
		private ListSequence<Double1D> postMus;
		private ListSequence<Double2D> postSgs;
		private SeqMVNormalDist mvnDists;

		public PostBeDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}
		
		private Double2D mvnVariates() throws ProbDistParmException {
			if (this.mvnDists == null) {
				this.mvnDists = new SeqMVNormalDist(this.getPostMus(), this
						.getPostSgs());
				return this.mvnDists.variate();
			} else {
				return this.mvnDists.variate(this.getPostMus(), this
						.getPostSgs());
			}
		}

		protected ListSequence<Double1D> getPostMus() {
			return this.postMus;
		}

		protected ListSequence<Double2D> getPostSgs() {
			return this.postSgs;
		}

		protected void setPostMus(ListSequence<Double1D> postMus) {
			this.postMus = postMus;
		}

		protected void setPostSgs(ListSequence<Double2D> postSgs) {
			this.postSgs = postSgs;
		}

		protected Double3D getSgs() {
			return (Double3D) this.chainParms[0];
		}

		protected Double2D getMus() {
			return (Double2D) this.chainParms[1];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[2];
		}

		protected Integer1D getB() {
			return (Integer1D) this.chainParms[3];
		}

		protected Double0D getNuB() {
			return (Double0D) this.chainParms[4];
		}

		protected Double1D getBe0() {
			return (Double1D) this.chainParms[5];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(6);
			this.chainParmNames.add("sg");
			this.chainParmNames.add("mu");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.chainParmNames.add("nub");
			this.chainParmNames.add("be0");
			this.fixedParmNames = new ArrayList<String>(0);
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(6);
			this.chainParmCheck.add(null); // FIXME
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.chainParmCheck.add(null);
			this.fixedParmCheck = new ArrayList<ProbDistParmCheck[]>(0);
			// Classes
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					6);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					0);
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

				Double0D postNuB = this.getNuB().plus(nKB1);
				this.getPostSgs().add(this.getSgs().get(k.value()).divide(postNuB));
				//FIXME
				this.getPostMus().add(sKB1.plus(this.getBe0().mult(this.getNuB())).mult(postNuB.recip()));
			}
		}

		@Override
		protected Double2D genVariate() throws ProbDistParmException {
			double[][] v = new double[this.getSgs().size()][MpgModel.this.getDims()];
			Double2D activeUpdates = this.mvnVariates();
			int i=0;
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], MpgModel.this.getDims());
				i++;
			}
			return new Double2D(v);
		}

		@Override
		protected double getDensity(Double2D pt) {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}	
	
	class PriorSgDist extends ProbDistInitializeByChain<Double3D> {
		private SeqInverseWishartDist iwDists;
		private ListSequence<Double2D> p0Seq;
		private ListSequence<Integer0D> k0Seq;

		public PriorSgDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}
		
		private ListSequence<Double2D> getP0Seq() {
			if (this.p0Seq == null) {
				this.p0Seq = new ListSequence<Double2D>();
				Integer1D active = this.getZ().items();
				for (int i = 0; i < active.size(); i++) {
					this.p0Seq.add(this.getP0());
				}
			}
			return this.p0Seq;
		}

		private ListSequence<Integer0D> getK0Seq() {
			if (this.k0Seq == null) {
				this.k0Seq = new ListSequence<Integer0D>();
				Integer1D active = this.getZ().items();
				for (int i = 0; i < active.size(); i++) {
					this.k0Seq.add(this.getK0());
				}
			}
			return this.k0Seq;
		}

		private Double2D getP0() {
			return (Double2D) this.fixedParms[0];
		}

		private Integer0D getK0() {
			return (Integer0D) this.fixedParms[1];
		}

		private Integer1D getZ() {
			return (Integer1D) this.chainParms[0];
		}

		private Double3D iwVariates() throws ProbDistParmException {
			if (this.iwDists == null) {
				this.iwDists = new SeqInverseWishartDist(this.getP0Seq(),
						this.getK0Seq());
				return this.iwDists.variate();
			} else {
				return this.iwDists.variate(this.getP0Seq(), this
						.getK0Seq());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(1);
			this.chainParmNames.add("Z");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("p0");
			this.fixedParmNames.add("k0");
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
			this.fixedParmClasses.add(Integer0D.class);
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
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}	

	class PostSgDist extends ProbDistInitializeByChain<Double3D> {
		private ListSequence<Double2D> postPs;
		private ListSequence<Integer0D> postKs;
		private SeqInverseWishartDist iwDists;

		public PostSgDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}
		
		protected ListSequence<Double2D> getPostPs() {
			return this.postPs;
		}

		protected ListSequence<Integer0D> getPostKs() {
			return this.postKs;
		}

		protected void setPostPs(ListSequence<Double2D> postPs) {
			this.postPs = postPs;
		}

		protected void setPostKs(ListSequence<Integer0D> postKs) {
			this.postKs = postKs;
		}

		protected Double2D getP0() {
			return (Double2D) this.fixedParms[0];
		}

		protected Integer0D getK0() {
			return (Integer0D) this.fixedParms[1];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[0];
		}

		protected Integer1D getB() {
			return (Integer1D) this.chainParms[1];
		}

		protected Double0D getNuB() {
			return (Double0D) this.chainParms[2];
		}

		protected Double0D getNuM() {
			return (Double0D) this.chainParms[3];
		}

		protected Double1D getMu0() {
			return (Double1D) this.chainParms[4];
		}

		protected Double1D getBe0() {
			return (Double1D) this.chainParms[5];
		}
		
		protected Double3D getSgs() {
			return (Double3D) this.chainParms[6];
		}

		private Double3D iwVariates() throws ProbDistParmException {
			if (this.iwDists == null) {
				this.iwDists = new SeqInverseWishartDist(this.getPostPs(),
						this.getPostKs());
				return this.iwDists.variate();
			} else {
				return this.iwDists.variate(this.getPostPs(), this
						.getPostKs());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(7);
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.chainParmNames.add("nub");
			this.chainParmNames.add("num");
			this.chainParmNames.add("mu0");
			this.chainParmNames.add("be0");
			this.chainParmNames.add("sg");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("p0");
			this.fixedParmNames.add("k0");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(7);
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
					7);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Double3D.class);				
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Integer0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			Double0D nuM = this.getNuM();
			Double0D nuB = this.getNuB();

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

				Double0D postNuB = nuB.plus(nKB1);
				Double0D postNuM = nuM.plus(nK).minus(nKB1.pow(2).divide(postNuB));

				this.getPostKs().add(this.getK0().plus(nK));

				Double2D s = xK.transpose().mult(xK);
				Double2D m0 = this.getMu0().outer(this.getMu0()).mult(nuM);
				Double2D b0 = this.getBe0().outer(this.getBe0()).mult(nuB);
				Double1D m = sK.plus(this.getMu0().mult(nuM));
				Double2D mm = m.outer(m).divide(nuM.plus(nK));
				Double1D n = m.mult(nKB1.divide(nuM.plus(nK))).minus(sKB1).plus(
						this.getBe0().mult(nuB));
				Double2D nn = n.outer(n).mult(nuM.plus(nK).divide(postNuM.mult(postNuB)));

				Double2D H = this.getP0().plus(s).plus(m0).plus(b0).minus(mm).minus(nn);
				this.getPostPs().add(H);
			}
		}

		@Override
		protected Double3D genVariate() throws ProbDistParmException {
			int d = MpgModel.this.getDims();
			double[][][] v = new double[this.getSgs().size()][d][d];
			Double3D activeUpdates = this.iwVariates();
			int i=0;
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], MpgModel.this.getDims());
				i++;
			}
			return new Double3D(v);	
		}

		@Override
		protected double getDensity(Double3D pt)
				throws ProbDistParmException {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}	
	
	class CPDist extends ProbDistInitializeByChain<Integer1D> {
		public CPDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double1D postProb;
		private Double1D postMu;
		private Double2D postSg;
		private CategoricalDistInitializeByP catDist;
		private MVNormalDist mvnDist;
		private HashMap<Integer0D,MVNormalDist> mvnDists; // Shouldn't you be using (and hacking) SeqProbDist for this? FIXME

		private ProbDist<Double2D> baseSg;
		private ProbDist<Double1D> baseMu;
		private ProbDist<Double1D> baseBe;

		public void setBaseSg(ProbDist<Double2D> baseSg) {
			this.baseSg = baseSg;
		}

		public ProbDist<Double2D> getBaseSg() {
			return this.baseSg;
		}

		public void setBaseMu(ProbDist<Double1D> baseMu) {
			this.baseMu = baseMu;
		}

		public ProbDist<Double1D> getBaseMu() {
			return this.baseMu;
		}

		public void setBaseBe(ProbDist<Double1D> baseBe) {
			this.baseBe = baseBe;
		}

		public ProbDist<Double1D> getBaseBe() {
			return this.baseBe;
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
				this.mvnDists = new HashMap<Integer0D,MVNormalDist>();
			}
			if (!this.mvnDists.containsKey(k)) {
				this.mvnDists.put(k, new MVNormalDist(this.getMus().get(k.value()), this.getSgs().get(k.value())));
			}
			// else do nothing
		}
		
		protected void resetPostMvn(Integer0D k) {
			if (this.mvnDists != null && this.mvnDists.containsKey(k)) {
				this.mvnDists.remove(k);
			}
		}
		
//		private Double1D getPostMu() {
//			return this.postMu;
//		}
//
//		private void setPostMu(Double1D mu) {
//			this.postMu = mu;
//		}
//
//		private Double2D getPostSg() {
//			return this.postSg;
//		}
//
//		private void setPostSg(Double2D sg) {
//			this.postSg = sg;
//		}

		private Double1D getPostProb() {
			return this.postProb;
		}

		protected void setPostProb(Double1D prob) {
			this.postProb = prob;
		}

		protected Double3D getSgs() {
			return (Double3D) this.chainParms[0];
		}

		protected Double2D getMus() {
			return (Double2D) this.chainParms[1];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[2];
		}

		protected Integer1D getB() {
			return (Integer1D) this.chainParms[3];
		}

		protected Double2D getBes() {
			return (Double2D) this.chainParms[4];
		}

		protected Double1D getMu0() {
			return (Double1D) this.chainParms[5];
		}

		protected Double1D getBe0() {
			return (Double1D) this.chainParms[6];
		}

		protected Double0D getNuM() {
			return (Double0D) this.chainParms[7];
		}

		private Double0D getNuB() {
			return (Double0D) this.chainParms[8];
		}

		protected Double0D getAl() {
			return (Double0D) this.chainParms[9];
		}

		protected Double2D getP0() {
			return (Double2D) this.fixedParms[0];
		}

		protected Integer0D getK0() {
			return (Integer0D) this.fixedParms[1];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(10);
			this.chainParmNames.add("sg");
			this.chainParmNames.add("mu");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.chainParmNames.add("be");
			this.chainParmNames.add("mu0");
			this.chainParmNames.add("be0");
			this.chainParmNames.add("num");
			this.chainParmNames.add("nub");
			this.chainParmNames.add("al");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("p0");
			this.fixedParmNames.add("k0");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(10);
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
					10);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Double1D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Integer0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			this.mvnDists = null;
		}

		@Override
		protected Integer1D genVariate() throws ProbDistParmException {
			// FIXME - This sampling scheme destructively modifies
			// mu, sg, be - HACK
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
					Double2D sgNew = this.getBaseSg().variateFast(this.getP0(),
							this.getK0());
					Double2D sgMu = sgNew.mult(1 / this.getNuM().value());
					Double1D muNew = this.getBaseMu().variateFast(
							this.getMu0(), sgMu);
					Double2D sgBe = sgNew.mult(1 / this.getNuB().value());
					Double1D beNew = this.getBaseBe().variateFast(
							this.getBe0(), sgBe);
					// Add it
					this.getSgs().set(zedNew.value(), sgNew);
					this.getMus().set(zedNew.value(), muNew);
					this.getBes().set(zedNew.value(), beNew);
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
					//this.setPostMu(this.getMus().get(k.value()));
					//this.setPostSg(this.getSgs().get(k.value()));
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
	
	class PostBe0Dist extends ProbDistInitializeByChain<Double1D> {
		protected MVNormalDist mvnDist;
		protected Double1D sg0InvB0;
		protected Double2D sg0Inv;
		protected Double1D postMu;
		protected Double2D postSg;
		
		public PostBe0Dist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}			

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(5);
			this.chainParmNames.add("sg");
			this.chainParmNames.add("be");
			this.chainParmNames.add("nub");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("b0");
			this.fixedParmNames.add("sb0");
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
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double1D.class);
			this.fixedParmClasses.add(Double2D.class);
		}

		protected Double1D mvnVariate() throws ProbDistParmException {
			if (this.mvnDist == null) {
				this.mvnDist = new MVNormalDist(this.getPostMu(), this
						.getPostSg());
				return this.mvnDist.variate();
			} else {
				return this.mvnDist.variate(this.getPostMu(), this
						.getPostSg());
			}
		}

		protected Double3D getSgs() {
			return (Double3D) this.chainParms[0];
		}

		protected Double2D getBes() {
			return (Double2D) this.chainParms[1];
		}

		protected Double0D getNuB() {
			return (Double0D) this.chainParms[2];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[3];
		}

		protected Integer1D getB() {
			return (Integer1D) this.chainParms[4];
		}

		protected Double1D getB0() {
			return (Double1D) this.fixedParms[0];
		}

		protected Double2D getSB0() {
			return (Double2D) this.fixedParms[1];
		}

		protected Double1D getS0InvB0() {
			if (this.sg0InvB0 == null) {
				this.sg0InvB0 = this.getS0Inv().mult(this.getB0());
			}
			return this.sg0InvB0;
		}

		protected Double2D getS0Inv() {
			if (this.sg0Inv == null) {
				this.sg0Inv = this.getSB0().inverse();
			}
			return this.sg0Inv;
		}

		@Override
		protected void setUpFromChainParms() {
			Double2D si = this.getS0Inv();
			Double1D sim = this.getS0InvB0();
			Integer1D active = this.getZ().items();
			Integer0D one = new Integer0D(1);
			for (Integer0D k : active) {
//				Integer1D zK = this.getZ().which(k);
//				Integer1D zKB1 = zK.intersect(this.getB().which(one));
//				if (zKB1.size() > 0) {
					Double2D pK = this.getSgs().get(k.value()).inverse();
					si = si.plus(pK.mult(this.getNuB()));
					sim = sim.plus(pK.mult(this.getNuB()).mult(
							this.getBes().get(k.value()))); // FIXME??
//				}
			}
			Double2D s = si.inverse();
			this.setPostSg(s);
			this.setPostMu(s.mult(sim));
		}

		protected Double1D getPostMu() {
			return this.postMu;
		}

		protected void setPostMu(Double1D mu) {
			this.postMu = mu;
		}

		protected Double2D getPostSg() {
			return this.postSg;
		}

		protected void setPostSg(Double2D sg) {
			this.postSg = sg;
		}

		@Override
		protected Double1D genVariate() throws ProbDistParmException {
			return this.mvnVariate();
		}

		@Override
		protected double getDensity(Double1D pt) {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}
	
	class PostNuBDist extends ProbDistInitializeByChain<Double0D> {
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
	
	class PostBDist extends ProbDistInitializeByChain<Integer1D> {
		private Double1D postProb;
		private Double1D postMu;
		private Double2D postSg;
		private CategoricalDistInitializeByP catDist;
		private MVNormalDist mvnDist;
		private HashMap<Integer0D,MVNormalDist> mvnDists;
		
		public PostBDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private Integer0D catVariate() throws ProbDistParmException {
			if (this.catDist == null) {
				this.catDist = new CategoricalDistInitializeByP(this
						.getPostProb());
				return this.catDist.variate();
			} else {
				return this.catDist.variate(this.getPostProb());
			}
		}

		private double mvnLogDensity(Integer0D k, Double1D pt)
				throws ProbDistParmException {
			return this.mvnDists.get(k).logDensity(pt);
		}
		
		private void setUpPostMvn(Integer0D k) throws ProbDistParmException {
			if (this.mvnDists == null) {
				this.mvnDists = new HashMap<Integer0D,MVNormalDist>();
			}
			if (!this.mvnDists.containsKey(k)) {
				this.mvnDists.put(k, new MVNormalDist(this.getMus().get(k.value()), this.getSgs().get(k.value())));
			}
			// else do nothing
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

		private Double1D getPostProb() {
			return this.postProb;
		}

		private void setPostProb(Double1D prob) {
			this.postProb = prob;
		}

		private Double3D getSgs() {
			return (Double3D) this.chainParms[0];
		}

		private Double2D getMus() {
			return (Double2D) this.chainParms[1];
		}

		private Integer1D getZ() {
			return (Integer1D) this.chainParms[2];
		}

		private Integer1D getB() {
			return (Integer1D) this.chainParms[3];
		}

		private Double2D getBes() {
			return (Double2D) this.chainParms[4];
		}

		private Double0D getBshape1() {
			return (Double0D) this.fixedParms[0];
		}

		private Double0D getBshape2() {
			return (Double0D) this.fixedParms[1];
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(5);
			this.chainParmNames.add("sg");
			this.chainParmNames.add("mu");
			this.chainParmNames.add("Z");
			this.chainParmNames.add("B");
			this.chainParmNames.add("be");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("bshape1");
			this.fixedParmNames.add("bshape2");
			// Checks
			this.chainParmCheck = new ArrayList<ProbDistParmCheck[]>(5);
			this.chainParmCheck.add(null); // FIXME
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
			this.chainParmClasses.add(Integer1D.class);
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

		@Override
		protected Integer1D genVariate() throws ProbDistParmException {
			int N = this.getSamplerData().size();
			Integer1D B = null;
			try {
				B = (Integer1D) this.getB().clone();
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace(); // FIXME
			}
			Integer0D unused = new Integer0D(-1);
			for (int i = 0; i < N; i++) {
				B.set(i, unused);

				Integer0D one = new Integer0D(1);
				Integer0D zedI = this.getZ().get(i);
				Integer1D zZedI = this.getZ().which(zedI);
				Integer1D zZedIB1 = zZedI.intersect(this.getB().which(one));

				int nZZedI = zZedI.size();
				int nZZedIB1 = zZedIB1.size();
				int nZZedIB0 = nZZedI - nZZedIB1;

				//this.setPostMu(this.getMus().get(zedI.value()));
				//this.setPostSg(this.getSgs().get(zedI.value()));
				this.setUpPostMvn(zedI);
				Double1D be = this.getBes().get(zedI.value());

				Double1D xIA = this.getSamplerData().get(i);
				Double1D xIAS = xIA.minus(be);
				double[] lpx = new double[2];
				lpx[0] = this.mvnLogDensity(zedI, xIA)
						+ Math.log(nZZedIB0 + this.getBshape1().value());
				lpx[1] = this.mvnLogDensity(zedI, xIAS)
						+ Math.log(nZZedIB1 + this.getBshape2().value());
				Double1D p = new Double1D(lpx).minus(
						Math.max(lpx[0], lpx[1])).exp();
				double sum = p.sum().value();
				this.setPostProb(p.mult(1/sum));
				
				B.set(i, this.catVariate());
			}
			return B;
		}

		@Override
		protected double getDensity(Integer1D pt) {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}


	public MpgModel(Map<String, Numeric<? extends Numeric<?>>> hypers,
			Map<String, Numeric<? extends Numeric<?>>> init, int dims)
			throws ProbDistParmException {
		super(hypers, dims);
		// Set up parameters
		this.params = new HashMap<String, RandomVar<? extends Numeric<?>>>();
		// Mu
		PriorMuDist priorMu = new PriorMuDist();
		PostMuDist postMu = new PostMuDist();
		Double2D defaultMu = (Double2D) init.get("mu"); // FIXME
		RandomVar<Double2D> rvMu = new RandomVar<Double2D>("mu", priorMu,
				postMu, defaultMu);
		this.params.put("mu", rvMu);

		// Sg
		PriorSgDist priorSg = new PriorSgDist((Double2D) this.getHyper("p0"), (Integer0D) this.getHyper("k0"));
		PostSgDist postSg = new PostSgDist((Double2D) this.getHyper("p0"), (Integer0D) this.getHyper("k0"));
		Double3D defaultSg = (Double3D) init.get("sg"); // FIXME
		RandomVar<Double3D> rvSg = new RandomVar<Double3D>("sg", priorSg,
				postSg, defaultSg);
		this.params.put("sg", rvSg);

		// Be
		Double2D defaultBe = (Double2D) init.get("be"); // FIXME
		PriorBeDist priorBe = new PriorBeDist();
		PostBeDist postBe = new PostBeDist();
		RandomVar<Double2D> rvBe = new RandomVar<Double2D>("be", priorBe,
				postBe, defaultBe);
		this.params.put("be", rvBe);

		// Z
		ProbDist<Integer1D> priorZ = new Deterministic<Integer1D>(init.get("Z"));
		CPDist postZ = new CPDist((Double2D) this.getHyper("p0"),
				(Integer0D) this.getHyper("k0")); // FIXME -- Hack!
		postZ.setBaseSg(new InverseWishartDist()); // FIXME -- Eek! Hack!
		postZ.setBaseMu(new MVNormalDist()); // FIXME -- Eek! Hack!
		postZ.setBaseBe(new MVNormalDist()); // FIXME -- Eek! Hack!
		Integer1D defaultZ = (Integer1D) init.get("Z");
		RandomVar<Integer1D> rvZ = new RandomVar<Integer1D>("Z", priorZ, postZ,
				defaultZ);
		this.params.put("Z", rvZ); // FIXME

		// B
		ProbDist<Integer1D> priorB = new Deterministic<Integer1D>(init.get("B"));
		PostBDist postB = new PostBDist(this.getHyper("bshape1"), this.getHyper("bshape2"));		
		Integer1D defaultB = (Integer1D) init.get("B"); // FIXME
		RandomVar<Integer1D> rvB = new RandomVar<Integer1D>("B", priorB, postB,
				defaultB);
		this.params.put("B", rvB); // FIXME
		// num
		ProbDist<Double0D> priorNuM = new GammaDist((Double0D) this
				.getHyper("numa"), (Double0D) this.getHyper("numb"));
		ProbDistInitializeByChain<Double0D> postNuM = new ProbDistInitializeByChain<Double0D>(
				(Double0D) this.getHyper("numa"), (Double0D) this
						.getHyper("numb")) {
			private GammaDist gammaDist;
			private Double0D postNuMA;
			private Double0D postNuMB;

			@Override
			protected void installParmChecks() {
				// Names
				this.chainParmNames = new ArrayList<String>(4);
				this.chainParmNames.add("sg");
				this.chainParmNames.add("mu");
				this.chainParmNames.add("Z");
				this.chainParmNames.add("mu0");
				this.fixedParmNames = new ArrayList<String>(2);
				this.fixedParmNames.add("numa");
				this.fixedParmNames.add("numb");
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
				this.chainParmClasses.add(Double3D.class);
				this.chainParmClasses.add(Double2D.class);
				this.chainParmClasses.add(Integer1D.class);
				this.chainParmClasses.add(Double1D.class);
				this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
						2);
				this.fixedParmClasses.add(Double0D.class);
				this.fixedParmClasses.add(Double0D.class);
			}

			private Double0D gammaVariate() throws ProbDistParmException {
				if (this.gammaDist == null) {
					this.gammaDist = new GammaDist(this.getPostNuMA(), this
							.getPostNuMB());
					return this.gammaDist.variate();
				} else {
					return this.gammaDist.variate(this.getPostNuMA(), this
							.getPostNuMB());
				}
			}

			protected Double3D getSgs() {
				return (Double3D) this.chainParms[0];
			}

			protected Double2D getMus() {
				return (Double2D) this.chainParms[1];
			}

			protected Integer1D getZ() {
				return (Integer1D) this.chainParms[2];
			}

			protected Double1D getMu0() {
				return (Double1D) this.chainParms[3];
			}

			protected Double0D getNuMA() {
				return (Double0D) this.fixedParms[0];
			}

			protected Double0D getNuMB() {
				return (Double0D) this.fixedParms[1];
			}

			@Override
			protected void setUpFromChainParms() {
				Double0D a = this.getNuMA();
				Double0D b = this.getNuMB();
				Integer1D active = this.getZ().items();
				for (Integer0D k : active) {
					Double1D m = this.getMus().get(k.value()).minus(
							this.getMu0());
					a = a.plus(0.5*((double) MpgModel.this.getDims()));
					b = b.plus(0.5*m.mult(this.getSgs().get(k.value()).inverse()).mult(m));
				}
				this.setPostNuMA(a);
				this.setPostNuMB(b);
			}

			private Double0D getPostNuMA() {
				return this.postNuMA;
			}

			private void setPostNuMA(Double0D nuMA) {
				this.postNuMA = nuMA;
			}

			private Double0D getPostNuMB() {
				return this.postNuMB;
			}

			private void setPostNuMB(Double0D nuMB) {
				this.postNuMB = nuMB;
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
		};
		Double0D defaultNuM = (Double0D) init.get("num"); // FIXME
		RandomVar<Double0D> rvNuM = new RandomVar<Double0D>("num", priorNuM,
				postNuM, defaultNuM);
		this.params.put("num", rvNuM); // FIXME

		// nub
		ProbDist<Double0D> priorNuB = new GammaDist((Double0D) this
				.getHyper("nuba"), (Double0D) this.getHyper("nubb"));
		Double0D defaultNuB = (Double0D) init.get("nub"); // FIXME
		PostNuBDist postNuB = new PostNuBDist((Double0D) this.getHyper("nuba"), (Double0D) this.getHyper("nubb"));
		RandomVar<Double0D> rvNuB = new RandomVar<Double0D>("nub", priorNuB,
				postNuB, defaultNuB);
		this.params.put("nub", rvNuB); // FIXME

		// mu0
		ProbDist<Double1D> priorMu0 = new MVNormalDist((Double1D) this
				.getHyper("m0"), (Double2D) this.getHyper("sm0"));
		ProbDistInitializeByChain<Double1D> postMu0 = new ProbDistInitializeByChain<Double1D>(
				(Double1D) this.getHyper("m0"), (Double2D) this.getHyper("sm0")) {
			private MVNormalDist mvnDist;
			private Double1D sg0InvM0;
			private Double2D sg0Inv;
			private Double1D postMu;
			private Double2D postSg;

			@Override
			protected void installParmChecks() {
				// Names
				this.chainParmNames = new ArrayList<String>(4);
				this.chainParmNames.add("sg");
				this.chainParmNames.add("mu");
				this.chainParmNames.add("Z");
				this.chainParmNames.add("num");
				this.fixedParmNames = new ArrayList<String>(2);
				this.fixedParmNames.add("m0");
				this.fixedParmNames.add("sm0");
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
				this.chainParmClasses.add(Double3D.class);
				this.chainParmClasses.add(Double2D.class);
				this.chainParmClasses.add(Integer1D.class);
				this.chainParmClasses.add(Double0D.class);
				this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
						2);
				this.fixedParmClasses.add(Double1D.class);
				this.fixedParmClasses.add(Double2D.class);
			}

			private Double1D mvnVariate() throws ProbDistParmException {
				if (this.mvnDist == null) {
					this.mvnDist = new MVNormalDist(this.getPostMu(), this
							.getPostSg());
					return this.mvnDist.variate();
				} else {
					return this.mvnDist.variate(this.getPostMu(), this
							.getPostSg());
				}
			}

			protected Double3D getSgs() {
				return (Double3D) this.chainParms[0];
			}

			protected Double2D getMus() {
				return (Double2D) this.chainParms[1];
			}

			protected Integer1D getZ() {
				return (Integer1D) this.chainParms[2];
			}

			protected Double0D getNuM() {
				return (Double0D) this.chainParms[3];
			}

			protected Double1D getM0() {
				return (Double1D) this.fixedParms[0];
			}

			protected Double2D getSM0() {
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
					this.sg0Inv = this.getSM0().inverse();
				}
				return this.sg0Inv;
			}

			@Override
			protected void setUpFromChainParms() {
				Double2D si = this.getS0Inv();
				Double1D sim = this.getS0InvM0();
				Integer1D active = this.getZ().items();
				for (Integer0D k : active) {
					Double2D pK = this.getSgs().get(k.value()).inverse();
					si = si.plus(pK.mult(this.getNuM()));
					sim = sim.plus(pK.mult(this.getNuM()).mult(
							this.getMus().get(k.value()))); // FIXME??
				}
				Double2D s = si.inverse();
				this.setPostSg(s);
				this.setPostMu(s.mult(sim));
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
			protected Double1D genVariate() throws ProbDistParmException {
				return this.mvnVariate();
			}

			@Override
			protected double getDensity(Double1D pt) {
				throw new UnsupportedOperationException(
						"Too lazy, come back later");
			}
		};
		Double1D defaultMu0 = (Double1D) init.get("mu0"); // FIXME
		RandomVar<Double1D> rvMu0 = new RandomVar<Double1D>("mu0", priorMu0,
				postMu0, defaultMu0);
		this.params.put("mu0", rvMu0);

		// be0
		ProbDist<Double1D> priorBe0 = new MVNormalDist((Double1D) this
				.getHyper("b0"), (Double2D) this.getHyper("sb0"));
		Double1D defaultBe0 = (Double1D) init.get("be0"); // FIXME
		PostBe0Dist postBe0 = new PostBe0Dist((Double1D) this.getHyper("b0"), (Double2D) this.getHyper("sb0"));
		RandomVar<Double1D> rvBe0 = new RandomVar<Double1D>("be0", priorBe0,
				postBe0, defaultBe0);
		this.params.put("be0", rvBe0); // FIXME

		// x
		ProbDist<Double0D> priorX = new BetaDist(
				(Double0D) this.getHyper("xa"), (Double0D) this.getHyper("xb"));
		ProbDistInitializeByChain<Double0D> postX = new ProbDistInitializeByChain<Double0D>() {
			private BetaDist betaDist;
			private Double0D postXA;
			private Double0D postXB;

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
					this.betaDist = new BetaDist(this.getPostXA(), this
							.getPostXB());
					return this.betaDist.variate();
				} else {
					return this.betaDist.variate(this.getPostXA(), this
							.getPostXB());
				}
			}

			protected Double0D getAl() {
				return (Double0D) this.chainParms[0];
			}

			@Override
			protected void setUpFromChainParms() {
				this.setPostXA(this.getAl().plus(1));
				this.setPostXB(new Double0D((double) this.getSamplerData()
						.size()));
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
				throw new UnsupportedOperationException(
						"Too lazy, come back later");
			}
		};
		Double0D defaultX = (Double0D) init.get("x"); // FIXME
		RandomVar<Double0D> rvX = new RandomVar<Double0D>("x", priorX, postX,
				defaultX);
		this.params.put("x", rvX); // FIXME

		// al
		ProbDist<Double0D> priorAl = new GammaDist((Double0D) this
				.getHyper("ala"), (Double0D) this.getHyper("alb"));
		ProbDistInitializeByChain<Double0D> postAl = new ProbDistInitializeByChain<Double0D>(
				(Double0D) this.getHyper("ala"), (Double0D) this
						.getHyper("alb")) {
			private GammaDist gammaDist;
			private CategoricalDistInitializeByP catDist;
			private Double0D postAlA;
			private Double0D postAlB;
			private Double1D postMix;

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
				throw new UnsupportedOperationException(
						"Too lazy, come back later");
			}
		};
		Double0D defaultAl = (Double0D) init.get("al"); // FIXME
		RandomVar<Double0D> rvAl = new RandomVar<Double0D>("al", priorAl,
				postAl, defaultAl);
		this.params.put("al", rvAl); // FIXME
	}

	@Override
	public ChainLink getInitialLink() {
		return new ChainLink(this.getParam("Z"), this.getParam("B"),
				this.getParam("be"), this.getParam("mu"), this.getParam("sg"), this
				.getParam("num"), this.getParam("nub"), this.getParam("mu0"),
				this.getParam("be0"), this.getParam("x"), this.getParam("al"));
	}
}
