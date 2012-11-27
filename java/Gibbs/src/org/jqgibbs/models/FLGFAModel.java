package org.jqgibbs.models;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import cern.colt.matrix.linalg.LUDecomposition;
import cern.jet.stat.Gamma;

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
import org.jqgibbs.mathstat.List;
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.mathstat.RandomVar;
import org.jqgibbs.mathstat.UndirectedGraph;
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
import org.jqgibbs.mathstat.probdist.SeqMVNormalDist;
import org.jqgibbs.mathstat.probdist.SeqMatrixNormalDist;

public class FLGFAModel extends Model {

	class PriorOmegaDist extends ProbDistInitializeByChain<Double3D> {
		private SeqInverseWishartDist iwDists; // FIXME - Shouldn't this use the
		// IID?
		private List<Double2D> repPhi;
		private List<Double0D> repLambda;

		public PriorOmegaDist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private List<Double2D> getRepPhi() {
			if (this.repPhi == null) {
				this.repPhi = new ArrayList<Double2D>();
				Integer1D active = this.getZ().items();
				for (int i = 0; i < active.size(); i++) {
					this.repPhi.add(this.getPhi());
				}
			}
			return this.repPhi;
		}

		private List<Double0D> getRepLambda() {
			if (this.repLambda == null) {
				this.repLambda = new List<Double0D>();
				Integer1D active = this.getZ().items();
				for (int i = 0; i < active.size(); i++) {
					this.repLambda.add(this.getLambda());
				}
			}
			return this.repLambda;
		}

		private Double2D getPhi() {
			return (Double2D) this.fixedParms[0];
		}

		private Double0D getLambda() {
			return (Double0D) this.fixedParms[1];
		}

		private Integer1D getZ() {
			return (Integer1D) this.chainParms[0];
		}

		private Double3D iwVariates() throws ProbDistParmException {
			if (this.iwDists == null) {
				this.iwDists = new SeqInverseWishartDist(this.getRepPhi(), this
						.getRepLambda());
				return this.iwDists.variate();
			} else {
				return this.iwDists.variate(this.getRepPhi(), this
						.getRepLambda());
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					1);
			this.chainParmClasses.add(Integer1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
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

	class PostOmegaDist extends ProbDistInitializeByChain<Double3D> {
		private List<Double2D> postPhis;
		private List<Double0D> postLambdas;
		private SeqInverseWishartDist iwDists;

		public PostOmegaDist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected List<Double2D> getPostPhis() {
			return this.postPhis;
		}

		protected List<Double0D> getPostLambdas() {
			return this.postLambdas;
		}

		protected void setPostPhis(List<Double2D> postPhis) {
			this.postPhis = postPhis;
		}

		protected void setPostLambdas(List<Double0D> postLambdas) {
			this.postLambdas = postLambdas;
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

		private Double3D iwVariates() throws ProbDistParmException {
			if (this.iwDists == null) {
				this.iwDists = new SeqInverseWishartDist(this.getPostPhis(),
						this.getPostLambdas());
				return this.iwDists.variate();
			} else {
				return this.iwDists.variate(this.getPostPhis(), this
						.getPostLambdas());
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					4);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
					2);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPostPhis(new List<Double2D>());
			this.setPostLambdas(new List<Double0D>());
			Integer1D active = this.getZ().items();
			int d = FLGFAModel.this.dims;
			for (Integer0D k : active) {
				Double2D AM = this.getA().get(k.value()).minus(this.getM());
				Double2D sgInv = this.getSg().get(k.value()).inverse();
				this.getPostLambdas().add(this.getLambda().plus(d));
				this.getPostPhis()
						.add(
								this.getPhi().plus(
										AM.mult(sgInv).mult(AM.transpose())));
			}
		}

		@Override
		protected Double3D genVariate() throws ProbDistParmException {
			int h = this.getM().size();
			double[][][] v = new double[this.getA().size()][h][h];
			Double3D activeUpdates = this.iwVariates();
			int i = 0;
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				v[k.value()] = Arrays.copyOf(activeUpdates.value()[i], h);
				i++;
			}
			return new Double3D(v);
		}

		@Override
		protected double getDensity(Double3D pt) throws ProbDistParmException {
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PriorSgDist extends ProbDistInitializeByChain<Double3D> {
		private SeqInverseWishartDist iwDists; // FIXME - Shouldn't this use the
		// IID?
		private List<Double2D> repPsi;
		private List<Double0D> repKappa;

		public PriorSgDist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private List<Double2D> getRepPsi() {
			if (this.repPsi == null) {
				this.repPsi = new List<Double2D>();
				Integer1D active = this.getZ().items();
				for (int i = 0; i < active.size(); i++) {
					this.repPsi.add(this.getPsi());
				}
			}
			return this.repPsi;
		}

		private List<Double0D> getRepKappa() {
			if (this.repKappa == null) {
				this.repKappa = new List<Double0D>();
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					1);
			this.chainParmClasses.add(Integer1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
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
		private List<Double2D> postPsi;
		private List<Double0D> postKappa;
		private SeqInverseWishartDist iwDists;

		public PostSgDist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
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
			return (Double2D) this.fixedParms[0];
		}

		protected Double0D getKappa() {
			return (Double0D) this.fixedParms[1];
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[0];
		}

		protected Integer2D getB() {
			return (Integer2D) this.chainParms[1];
		}

		protected Double3D getOmega() {
			return (Double3D) this.chainParms[2];
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					5);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
					2);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPostPsi(new List<Double2D>());
			this.setPostKappa(new List<Double0D>());
			Integer1D active = this.getZ().items();
			int h = this.getA().size();
			for (Integer0D k : active) {
				Integer1D zK = this.getZ().which(k);
				Double2D xK = this.getSamplerData().getAll(zK.value());
				Integer0D nK = new Integer0D(xK.size());

//				this.getPostKappa().add(this.getKappa().plus(nK).plus(h));
				this.getPostKappa().add(this.getKappa().plus(nK));

				Integer2D B = this.getB().getAll(zK.value());
				Double2D BTB = B.transpose().mult(B);
				Double2D omegaInv = this.getOmega().get(k.value()).inverse();
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
			int d = FLGFAModel.this.dims;
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

		public PriorADist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private Double3D getActiveOmega() {
			Double3D omega = (Double3D) this.chainParms[0];
			Integer1D active = this.getZ().items();
			return omega.getAll(active.value());
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
					FLGFAModel.this.dims);
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
				this.matnDists = new SeqMatrixNormalDist(this.getRepM(), this
						.getActiveSg(), this.getActiveOmega());
				return this.matnDists.variate();
			} else {
				return this.matnDists.variate(this.getRepM(), this
						.getActiveSg(), this.getActiveOmega());
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					5);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
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
		private List<Double2D> postMs;
		private List<Double2D> postSgs;
		private List<Double2D> postOmegas;
		private SeqMatrixNormalDist matnDists;

		public PostADist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
		}

		private Double3D matnVariates() throws ProbDistParmException {
			Integer1D active = this.getZ().items();
			Double3D sgActive = this.getSg().getAll(active.value());
			if (this.matnDists == null) {
				this.matnDists = new SeqMatrixNormalDist(this.getPostM(),
						sgActive, this.getPostOmega());
				// FIXME - Will this work if all arguments are not
				// Lists?
				return this.matnDists.variate();
			} else {
				return this.matnDists.variate(this.getPostM(), sgActive, this
						.getPostOmega());
			}
		}

		protected List<Double2D> getPostM() {
			return this.postMs;
		}

		protected List<Double2D> getPostOmega() {
			return this.postOmegas;
		}

		protected void setPostM(List<Double2D> postMs) {
			this.postMs = postMs;
		}

		protected void setPostOmega(List<Double2D> postOmegas) {
			this.postOmegas = postOmegas;
		}

		protected Double3D getOmega() {
			return (Double3D) this.chainParms[0];
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					5);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
					0);
		}

		@Override
		protected void setUpFromChainParms() {
			this.setPostOmega(new List<Double2D>());
			this.setPostM(new List<Double2D>());
			Integer1D active = this.getZ().items();
			for (Integer0D k : active) {
				Integer1D zK = this.getZ().which(k);
				Double2D xK = this.getSamplerData().getAll(zK.value());
				Integer2D B = this.getB().getAll(zK.value());
				Double2D BTB = B.transpose().mult(B);
				Double2D omegaInv = this.getOmega().get(k.value()).inverse();
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
			int d = FLGFAModel.this.dims;
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
					for (int m=0; m<h; m++) {
						for (int n=0; n<d; n++) {
							if (Double.isInfinite(v[k.value()][m][n]) || Double.isNaN(v[k.value()][m][n])) {
								System.err.println("Warning: bad sample for A: " + activeUpdates.get(i));
								System.err.println("Sg: " + this.getSg().get(k.value()));
								System.err.println("Omega: " + this.getOmega().get(k.value()));
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

	// Neal's (1998) Algorithm 8
	class CPDist extends ProbDistInitializeByChain<Integer1D> {
		public CPDist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Double1D postProb;
		private CategoricalDistInitializeByP catDist;
		private HashMap<Integer0D, MVNormalDist> mvnDists; // Shouldn't you be
		// using (and
		// hacking)
		// SeqProbDist for
		// this? FIXME

		private ProbDist<Double2D> baseOmega;
		private ProbDist<Double2D> baseSg;
		private ProbDist<Double2D> baseA;

		public ProbDist<Double2D> getBaseOmega() {
			return this.baseOmega;
		}

		public void setBaseOmega(ProbDist<Double2D> baseOmega) {
			this.baseOmega = baseOmega;
		}

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
				double[] zeroD = new double[FLGFAModel.this.dims];
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

		protected Double3D getOmega() {
			return (Double3D) this.chainParms[0];
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

		protected Double0D getKappa() {
			return (Double0D) this.fixedParms[1];
		}

		protected Double2D getPhi() {
			return (Double2D) this.fixedParms[2];
		}

		protected Double0D getLambda() {
			return (Double0D) this.fixedParms[3];
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					7);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.chainParmClasses.add(Double2D.class);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
					4);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
			this.fixedParmClasses.add(Double2D.class);
			this.fixedParmClasses.add(Double0D.class);
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
					// Double2D omegaNew =
					// this.getBaseOmega().variateFast(this.getPhi(),
					// this.getLambda());
					// Double2D sgNew =
					// this.getBaseSg().variateFast(this.getPsi(),
					// this.getKappa());
					Double2D omegaNew = null;
					Double2D sgNew = null;
					Double2D kronSg = null;
					boolean ksOk = false;
					while (!ksOk) {
						omegaNew = this.getBaseOmega().variateFast();
						sgNew = this.getBaseSg().variateFast();
						kronSg = sgNew.kron(omegaNew);
						ksOk = kronSg.isWellConditioned();
						if (!ksOk) {
							System.err.println("Warning: threw out IW samples that had ill-conditioned Kronecker product");
						}						
					} 
					Double2D aNew = this.getBaseA().variateFast(this.getM(),
							sgNew, omegaNew);
					// Add it
					this.getOmega().set(zedNew.value(), omegaNew);
					this.getSg().set(zedNew.value(), sgNew);
					this.getA().set(zedNew.value(), aNew);
					active.add(zedNew);
				}
				// Select a category for this point
				double[] logP = new double[active.size()]; // FIXME - Probably
				// you could speed
				// this up by a lot
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
			throw new UnsupportedOperationException("Too lazy, come back later");
		}
	}

	class PostMDist extends ProbDistInitializeByChain<Double2D> {
		protected MVNormalDist mvnDist;
		protected Double1D vecW;
		protected Double2D SInv;
		protected Double1D postVecMu;
		protected Double2D postKronSg;

		public PostMDist(Numeric... fixed) throws ProbDistParmException {
			super(fixed);
		}

		protected Integer1D getZ() {
			return (Integer1D) this.chainParms[0];
		}

		protected Integer2D getB() {
			return (Integer2D) this.chainParms[1];
		}

		protected Double3D getOmega() {
			return (Double3D) this.chainParms[2];
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					5);
			this.chainParmClasses.add(Integer1D.class);
			this.chainParmClasses.add(Integer2D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.chainParmClasses.add(Double3D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
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
				Double2D omegaInv = this.getOmega().get(k.value()).inverse();
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

		public PostXDist(Numeric... fixed) throws ProbDistParmException {
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					1);
			this.chainParmClasses.add(Double0D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
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

		public PostAlDist(Numeric... fixed) throws ProbDistParmException {
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					2);
			this.chainParmClasses.add(Double0D.class);
			this.chainParmClasses.add(Integer1D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
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

	class PostBDist extends ProbDistInitializeByChain<Integer2D> {
		public PostBDist(Numeric... fixed) throws ProbDistParmException {
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
			this.chainParmClasses = new ArrayList<Class<? extends Numeric>>(
					1);
			this.chainParmClasses.add(Integer2D.class);
			this.fixedParmClasses = new ArrayList<Class<? extends Numeric>>(
					0);
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
			// M
			RandomVar<Double2D> rvM = (RandomVar<Double2D>) cli.get("M");
			Double2D M = rvM.getNumericValue();
			apost += rvM.getPrior().logDensity(M);
			// Hyperparameters
			RandomVar<Double0D> rvAl = (RandomVar<Double0D>) cli.get("al");
			Double0D al = rvAl.getNumericValue();
			apost += rvAl.getPrior().logDensity(al);
			// Model parameters and likelihood
			RandomVar<Double3D> rvOmega = (RandomVar<Double3D>) cli
					.get("Omega");
			RandomVar<Double3D> rvSg = (RandomVar<Double3D>) cli.get("Sg");
			RandomVar<Double3D> rvA = (RandomVar<Double3D>) cli.get("A");
			// B
			RandomVar<Integer2D> rvB = (RandomVar<Integer2D>) cli.get("B");
			Integer2D B = rvB.getNumericValue();
			// skip
			// Z
			RandomVar<Integer1D> rvZ = (RandomVar<Integer1D>) cli.get("Z");
			Integer1D Z = rvZ.getNumericValue();
			for (int j = 0; j < Z.size(); j++) {
				int[] zIndices = new int[j];
				for (int k=0; k<j; k++) {
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
				Double2D Omega = ((Double3D) rvOmega.getNumericValue()).get(z
						.value());
				Double2D Sg = ((Double3D) rvSg.getNumericValue())
						.get(z.value());
				Double2D A = ((Double3D) rvA.getNumericValue()).get(z.value());
				CPDist postZ = (CPDist) rvZ.getPosterior();
				apost += postZ.getBaseOmega().logDensity(Omega);
				apost += postZ.getBaseSg().logDensity(Sg);
				apost += postZ.getBaseA().logDensity(A, M, Sg, Omega);
				// Likelihood
				Integer1D whichz = Z.which(z);
				Integer2D db = B.getAll(whichz.value());
				Double2D dz = d.getAll(whichz.value());
				for (int j = 0; j < dz.size(); j++) {
					Double1D pt = dz.get(j);
					Integer1D b = db.get(j);
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
	// List<Integer2D> BSeq = new List<Integer2D>();
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
	// Numeric>> vars = new HashMap<String, AbstractSequence<? extends
	// AbstractSequence<?, ?>, ? extends Numeric>>();
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
	// ((List<Double3D>) vars.get("A"))
	// .add((Double3D) next.get("A").getNumericValue()
	// .clone());
	// ((List<Double3D>) vars.get("Sg"))
	// .add((Double3D) next.get("Sg").getNumericValue()
	// .clone());
	// ((List<Double3D>) vars.get("Omega"))
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

	public FLGFAModel() {
		super();
	}

	public FLGFAModel(Map<String, Numeric> hypers,
			Map<String, Numeric> init, int dims)
			throws ProbDistParmException {
		super(hypers, dims);
		
		// Set up parameters
		this.params = new HashMap<String, RandomVar<? extends Numeric>>();

		
		// Omega
		PriorOmegaDist priorOmega = new PriorOmegaDist((Double2D) this
				.getHyper("Phi"), (Double0D) this.getHyper("lambda"));
		PostOmegaDist postOmega = new PostOmegaDist((Double2D) this
				.getHyper("Phi"), (Double0D) this.getHyper("lambda"));
		Double3D defaultOmega = (Double3D) init.get("Omega"); // FIXME
		RandomVar<Double3D> rvOmega = new RandomVar<Double3D>("Omega",
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
		CPDist postZ = new CPDist((Double2D) this.getHyper("Psi"),
				(Double0D) this.getHyper("kappa"), (Double2D) this
						.getHyper("Phi"), (Double0D) this.getHyper("lambda")); // FIXME
		// --
		// Hack!
		postZ.setBaseOmega(new InverseWishartDist(this.getHyper("Phi"), this
				.getHyper("lambda"))); // FIXME -- Eek! Hack!
		postZ.setBaseSg(new InverseWishartDist(this.getHyper("Psi"), this
				.getHyper("kappa"))); // FIXME -- Eek! Hack!
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
		try {
		ProbDist<Double2D> priorM = new MatrixNormalDist((Double2D) this
				.getHyper("W"), (Double2D) this.getHyper("S"), Double2D
				.ident(defaultM.numRows()));
		PostMDist postM = new PostMDist((Double2D) this.getHyper("W"),
				(Double2D) this.getHyper("S"));
		RandomVar<Double2D> rvM = new RandomVar<Double2D>("M", priorM, postM,
				defaultM);
		this.params.put("M", rvM);
		
		} catch (ProbDistParmException e) {
			System.err.println(e.getMessage());
			throw e;
		}

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
