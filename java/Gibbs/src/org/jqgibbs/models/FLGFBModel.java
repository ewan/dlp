package org.jqgibbs.models;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.jqgibbs.ChainLink;
import org.jqgibbs.Model;
import org.jqgibbs.mathstat.AbstractSequence;
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

public class FLGFBModel extends Model {

	class PriorOmegaDist extends ProbDistMC<Double2D> {
		private InverseWishartDist iwDist; // FIXME - Shouldn't this use the IID?


		public PriorOmegaDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private Double2D getPhi() {
			return (Double2D) this.fixedParms[0];
		}

		private Integer0D getLambda() {
			return (Integer0D) this.fixedParms[1];
		}

		private Integer1D getZ() {
			return (Integer1D) this.chainParms[0];
		}

		private Double2D iwVariate() throws ProbDistParmException {
			if (this.iwDist == null) {
				this.iwDist = new InverseWishartDist(this.getPhi(), this.getLambda());
				return this.iwDist.variate();
			} else {
				return this.iwDist.variate(this.getPhi(), this.getLambda());
			}
		}

		@Override
		protected void installParmChecks() {
			// Names
			this.chainParmNames = new ArrayList<String>(1);
			this.chainParmNames.add("Z");
			this.fixedParmNames = new ArrayList<String>(2);
			this.fixedParmNames.add("Phi");
			this.fixedParmNames.add("lambda");
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
		protected Double2D genVariate() throws ProbDistParmException {
			return this.iwVariate();
		}

		@Override
		protected double getDensity(Double2D pt) {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}	

	class PostOmegaDist extends ProbDistMC<Double2D> {
		private Double2D postPhi;
		private Integer0D postLambda;
		private InverseWishartDist iwDist;

		public PostOmegaDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}
		
		protected Double2D getPostPhi() {
			return this.postPhi;
		}

		protected Integer0D getPostLambda() {
			return this.postLambda;
		}

		protected void setPostPhi(Double2D postPhi) {
			this.postPhi = postPhi;
		}

		protected void setPostLambda(Integer0D postLambda) {
			this.postLambda = postLambda;
		}

		protected Double2D getPhi() {
			return (Double2D) this.fixedParms[0];
		}

		protected Integer0D getLambda() {
			return (Integer0D) this.fixedParms[1];
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
				this.iwDist = new InverseWishartDist(this.getPostPhi(),
						this.getPostLambda());
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
			this.fixedParmClasses.add(Integer0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPostPhi(this.getPhi());
			Integer1D active = this.getZ().items();
			int d = FLGFBModel.this.dims;
			this.setPostLambda(this.getLambda().plus(d*active.size()));			
			for (Integer0D k : active) {
				Double2D AM = this.getA().get(k.value()).minus(this.getM());
				Double2D sgInv = this.getSg().get(k.value()).inverse();
				this.setPostPhi(this.getPostPhi().plus(AM.mult(sgInv).mult(AM.transpose())));
			}
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
	}	

	class PriorSgDist extends ProbDistMC<Double3D> {
		private SeqInverseWishartDist iwDists; // FIXME - Shouldn't this use the IID?
		private ListSequence<Double2D> repPsi;
		private ListSequence<Integer0D> repKappa;

		public PriorSgDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}
		
		private ListSequence<Double2D> getRepPsi() {
			if (this.repPsi == null) {
				this.repPsi = new ListSequence<Double2D>();
				Integer1D active = this.getZ().items();
				for (int i = 0; i<active.size(); i++) {
					this.repPsi.add(this.getPsi());
				}
			}
			return this.repPsi;
		}

		private ListSequence<Integer0D> getRepKappa() {
			if (this.repKappa == null) {
				this.repKappa = new ListSequence<Integer0D>();
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

		private Integer0D getKappa() {
			return (Integer0D) this.fixedParms[1];
		}

		private Integer1D getZ() {
			return (Integer1D) this.chainParms[0];
		}

		private Double3D iwVariates() throws ProbDistParmException {
			if (this.iwDists == null) {
				this.iwDists = new SeqInverseWishartDist(this.getRepPsi(), this.getRepKappa());
				return this.iwDists.variate();
			} else {
				return this.iwDists.variate(this.getRepKappa(), this.getRepKappa());
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
	
	class PostSgDist extends ProbDistMC<Double3D> {
		private ListSequence<Double2D> postPsi;
		private ListSequence<Integer0D> postKappa;
		private SeqInverseWishartDist iwDists;

		public PostSgDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}
		
		protected ListSequence<Double2D> getPostPsi() {
			return this.postPsi;
		}

		protected ListSequence<Integer0D> getPostKappa() {
			return this.postKappa;
		}

		protected void setPostPsi(ListSequence<Double2D> postPsi) {
			this.postPsi = postPsi;
		}

		protected void setPostKappa(ListSequence<Integer0D> postKappa) {
			this.postKappa = postKappa;
		}

		protected Double2D getPsi() {
			return (Double2D) this.fixedParms[0];
		}

		protected Integer0D getKappa() {
			return (Integer0D) this.fixedParms[1];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[0];
		}

		protected Integer2D getB() {
			return (Integer2D) this.chainParms[1];
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
			this.chainParmNames = new ArrayList<String>(5);
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
			this.chainParmClasses.add(Integer2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);				
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(
					2);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Integer0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPostPsi(new ListSequence<Double2D>());
			this.setPostKappa(new ListSequence<Integer0D>());
			Integer1D active = this.getZ().items();
			int h = this.getA().size();
			for (Integer0D k : active) {
				Integer1D zK = this.getZ().which(k);
				Double2D xK = this.getSamplerData().getAll(zK.value());				
				Integer0D nK = new Integer0D(xK.size());
				
				this.getPostKappa().add(this.getKappa().plus(nK).plus(h));
				
				Integer2D B = this.getB().getAll(zK.value());	
				Double2D BTB = B.transpose().mult(B);
				Double2D omegaInv = this.getOmega().inverse();
				Double2D MTOmegaInv = this.getM().transpose().mult(omegaInv);
				Double2D MTOmegaInvM = MTOmegaInv.mult(this.getM());
				Double2D XTB = B.transpose().mult(xK).transpose();
				Double2D Mp = XTB.plus(MTOmegaInv);
				Double2D postOmega = omegaInv.plus(BTB).inverse();
				Double2D MpPOmegaMpT = Mp.mult(postOmega).mult(Mp.transpose()); 
				Double2D XTX = xK.transpose().mult(xK);
				this.getPostPsi().add(this.getPsi().plus(XTX).plus(MTOmegaInvM).minus(MpPOmegaMpT));
			}
		}

		@Override
		protected Double3D genVariate() throws ProbDistParmException {
			int d = FLGFBModel.this.dims;
			double[][][] v = new double[this.getA().size()][d][d];
			Double3D activeUpdates = this.iwVariates();
			int i=0;
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], d);
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
	
	class PriorADist extends ProbDistMC<Double3D> {
		private SeqMatrixNormalDist matnDists;
		
		public PriorADist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}
		
		private Double3D getRepOmega() {
			Double2D Omega = (Double2D) this.chainParms[0];
			Integer1D active = this.getZ().items();
			Double3D repOmega = new Double3D(active.size(),this.getB().numCols(),this.getB().numCols());				
			for (int i = 0; i < active.size(); i++) {
				repOmega.add(Omega);
			}
			return repOmega;			
		}		

		private Double3D getActiveSg() {
			Double3D sg =  (Double3D) this.chainParms[1];
			Integer1D active = this.getZ().items();	
			return sg.getAll(active.value());
		}

		private Double3D getRepM() {
			Double2D M = (Double2D) this.chainParms[2];
			Integer1D active = this.getZ().items();
			Double3D repM = new Double3D(active.size(),this.getB().numCols(),FLGFBModel.this.dims);				
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
				this.matnDists = new SeqMatrixNormalDist(this.getRepM(), this.getActiveSg(), this.getRepOmega());
				return this.matnDists.variate();
			} else {
				return this.matnDists.variate(this.getRepM(), this.getActiveSg(), this.getRepOmega());
			}
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(5);
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
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}	

	class PostADist extends ProbDistMC<Double3D> {
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
				this.matnDists = new SeqMatrixNormalDist(this.getPostM(), sgActive, this.getPostOmega());
				// FIXME - Will this work if all arguments are not ListSequences?
				return this.matnDists.variate();
			} else {
				return this.matnDists.variate(this.getPostM(), sgActive, this.getPostOmega());
			}
		}

		protected ListSequence<Double2D> getPostM() {
			return this.postMs;
		}
		
		protected ListSequence<Double2D> getPostOmega() {
			return this.postOmegas;
		}

		protected void setPostM(ListSequence<Double2D> postMs) {
			this.postMs = postMs;
		}

		protected void setPostOmega(ListSequence<Double2D> postOmegas) {
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

		protected Integer2D getB() {
			return (Integer2D) this.chainParms[4];
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(5);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(0);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPostOmega(new ListSequence<Double2D>());
			this.setPostM(new ListSequence<Double2D>());			
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				Integer1D zK = this.getZ().which(k);
				Double2D xK = this.getSamplerData().getAll(zK.value());				
				Integer2D B = this.getB().getAll(zK.value());
				Double2D BTB = B.transpose().mult(B);
				Double2D omegaInv = this.getOmega().inverse();
				Double2D postOmegaInv = omegaInv.plus(BTB);
				Double2D postOmega = postOmegaInv.inverse();
				this.getPostOmega().add(postOmega);
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
			int d = FLGFBModel.this.dims;
			double[][][] v = new double[K][h][d];
			Double3D activeUpdates = this.matnVariates();
			int i=0;
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], h); // FIXME?? Check dims	
				i++;
			}
			return new Double3D(v);
		}

		@Override
		protected double getDensity(Double3D pt) {
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}	
	
	class CPDist extends ProbDistMC<Integer1D> {
		public CPDist(Numeric<?>... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double1D postProb;
		private CategoricalDist catDist;
		private HashMap<Integer0D,MVNormalDist> mvnDists; // Shouldn't you be using (and hacking) SeqProbDist for this? FIXME

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
				this.mvnDists = new HashMap<Integer0D,MVNormalDist>();
			}
			if (!this.mvnDists.containsKey(k)) {
				double[] zeroD = new double[FLGFBModel.this.dims];
				Double1D zero = new Double1D(zeroD);
				this.mvnDists.put(k, new MVNormalDist(zero, this.getSg().get(k.value())));
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

		protected Integer2D getB() {
			return (Integer2D) this.chainParms[4];
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

		protected Integer0D getKappa() {
			return (Integer0D) this.fixedParms[1];
		}
		
		protected Double2D getPhi() {
			return (Double2D) this.fixedParms[2];
		}

		protected Integer0D getLambda() {
			return (Integer0D) this.fixedParms[3];
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(7);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(4);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Integer0D.class);
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
					Double2D sgNew = this.getBaseSg().variate(this.getPsi(), this.getKappa());
					Double2D aNew = this.getBaseA().variate(this.getM(), sgNew, this.getOmega());
					// Add it
					this.getSg().set(zedNew.value(), sgNew);
					this.getA().set(zedNew.value(), aNew);
					active.add(zedNew);
				}
				// Select a category for this point
				double[] logP = new double[active.size()]; // FIXME - Probably you could speed this up by a lot
				double maxLogP = -Double.MAX_VALUE; // FIXME??
				int ai = 0;
				for (Integer0D k : active) {
					Integer1D zK = Z.which(k);
					Integer1D b = this.getB().get(i);
					Double1D aTB = b.mult(this.getA().get(k.value()));
					Double1D xIA = this.getSamplerData().get(i).minus(aTB);					
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
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}	

	class PostMDist extends ProbDistMC<Double2D> {
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

		protected Integer2D getB() {
			return (Integer2D) this.chainParms[1];
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
			this.chainParmClasses.add(Integer2D.class);
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
			for (Integer0D k : active) {
				Double2D sgInv = this.getSg().get(k.value()).inverse();
				Double2D omegaInv = this.getOmega().inverse();
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
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
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
				this.catDist = new CategoricalDist(this
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric<?>>>(1);
			this.chainParmClasses.add(Integer2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric<?>>>(0);
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
			throw new UnsupportedOperationException(
					"Too lazy, come back later");
		}
	}	
	
	public FLGFBModel(Map<String, Numeric<? extends Numeric<?>>> hypers,
			Map<String, Numeric<? extends Numeric<?>>> init, int dims)
			throws ProbDistParmException {
		super(hypers, dims);
		// Set up parameters
		this.params = new HashMap<String, RandomVar<? extends Numeric<?>>>();

		// Omega
		PriorOmegaDist priorOmega = new PriorOmegaDist((Double2D) this.getHyper("Phi"), (Integer0D) this.getHyper("lambda"));
		PostOmegaDist postOmega = new PostOmegaDist((Double2D) this.getHyper("Phi"), (Integer0D) this.getHyper("lambda"));
		Double2D defaultOmega = (Double2D) init.get("Omega"); // FIXME
		RandomVar<Double2D> rvOmega = new RandomVar<Double2D>("Omega", priorOmega, postOmega, defaultOmega);
		this.params.put("Omega", rvOmega);		

		// Sg
		PriorSgDist priorSg = new PriorSgDist((Double2D) this.getHyper("Psi"), (Integer0D) this.getHyper("kappa"));
		PostSgDist postSg = new PostSgDist((Double2D) this.getHyper("Psi"), (Integer0D) this.getHyper("kappa"));
		Double3D defaultSg = (Double3D) init.get("Sg"); // FIXME
		RandomVar<Double3D> rvSg = new RandomVar<Double3D>("Sg", priorSg, postSg, defaultSg);
		this.params.put("Sg", rvSg);
		
		// A
		PriorADist priorA = new PriorADist();
		PostADist postA = new PostADist();
		Double3D defaultA = (Double3D) init.get("A"); // FIXME
		RandomVar<Double3D> rvA = new RandomVar<Double3D>("A", priorA, postA, defaultA);
		this.params.put("A", rvA);

		// Z
		ProbDist<Integer1D> priorZ = new Deterministic<Integer1D>(init.get("Z"));
		CPDist postZ = new CPDist((Double2D) this.getHyper("Psi"),
				(Integer0D) this.getHyper("kappa"), (Double2D) this.getHyper("Phi"),
				(Integer0D) this.getHyper("lambda")); // FIXME -- Hack!
		postZ.setBaseSg(new InverseWishartDist()); // FIXME -- Eek! Hack!
		postZ.setBaseA(new MatrixNormalDist()); // FIXME -- Eek! Hack!
		Integer1D defaultZ = (Integer1D) init.get("Z");
		RandomVar<Integer1D> rvZ = new RandomVar<Integer1D>("Z", priorZ, postZ,
				defaultZ);
		this.params.put("Z", rvZ); // FIXME

		// B
		ProbDist<Integer2D> priorB = new Deterministic<Integer2D>(init.get("B"));
		PostBDist postB = new PostBDist();		
		Integer2D defaultB = (Integer2D) init.get("B"); // FIXME
		RandomVar<Integer2D> rvB = new RandomVar<Integer2D>("B", priorB, postB,
				defaultB);
		this.params.put("B", rvB); // FIXME
		
		// M
		Double2D defaultM = (Double2D) init.get("M"); // FIXME		
		ProbDist<Double2D> priorM = new MatrixNormalDist((Double2D) this
				.getHyper("W"), (Double2D) this.getHyper("S"), Double2D.ident(defaultM.numRows()));
		PostMDist postM = new PostMDist((Double2D) this
				.getHyper("W"), (Double2D) this.getHyper("S"));
		RandomVar<Double2D> rvM = new RandomVar<Double2D>("M", priorM, postM, defaultM);
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
		ProbDistMC<Double0D> postAl = new PostAlDist((Double0D) this
				.getHyper("ala"), (Double0D) this.getHyper("alb"));		
		Double0D defaultAl = (Double0D) init.get("al"); // FIXME
		RandomVar<Double0D> rvAl = new RandomVar<Double0D>("al", priorAl,
				postAl, defaultAl);
		this.params.put("al", rvAl); // FIXME
	}	

	@Override
	public ChainLink getInitialLink() {
		return new ChainLink(this.getParam("Z"), this.getParam("B"),
				this.getParam("A"), this.getParam("Sg"), this.getParam("Omega"),
				this.getParam("M"), this.getParam("x"), this.getParam("al"));
	}	
}
