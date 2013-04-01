package org.jqgibbs.models;

import java.util.HashMap;
import java.util.Map;

import org.jqgibbs.Chain;
import org.jqgibbs.ChainLink;
import org.jqgibbs.Flattenable;
import org.jqgibbs.Model;
import org.jqgibbs.RandomVar;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Double3D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Integer1D;
import org.jqgibbs.mathstat.probdist.BetaDist;
import org.jqgibbs.mathstat.probdist.CategoricalDist;
import org.jqgibbs.mathstat.probdist.GammaDist;
import org.jqgibbs.mathstat.probdist.InverseWishartDist;
import org.jqgibbs.mathstat.probdist.MVNormalDist;
import org.jqgibbs.mathstat.probdist.MatrixNormalDist;
import org.jqgibbs.mathstat.probdist.ProbDistMC;

import cern.colt.matrix.impl.DenseDoubleMatrix3D;
import cern.jet.stat.Gamma;

public class MLM_sample_params extends Model {

	private static double lΓ(double x, int p) {
		double r = 0.25*p*(p-1)*Math.log(Math.PI);
		for (int i=0; i<p; i++) {
			r += Gamma.logGamma(x+i/2);
		}
		return r;
	}
	
	public static Map<String,Object> extend_hypers(Map<String,Flattenable> hypers, Double2D Y, Double2D X) {
		Map<String,Object> extended_hypers = new HashMap<String,Object>();
		
		Double0D κ = (Double0D) hypers.get("kappa");
		extended_hypers.put("κ", κ);
		extended_hypers.put("κ_1", κ.plus(1));
		
		Double2D Ψ = (Double2D) hypers.get("Psi");
		extended_hypers.put("Ψ", Ψ);
		extended_hypers.put("ldet_Ψ", 0.5*κ.value()*Math.log(Ψ.det()));
		
		Double2D Φ = (Double2D) hypers.get("Phi");
		Double0D λ = (Double0D) hypers.get("lambda");
		extended_hypers.put("Φ", Φ);
		extended_hypers.put("λ", λ);
		
		int h = X.numCols();
		extended_hypers.put("h", h);
		
		Double2D Sⁿ = ((Double2D) hypers.get("S")).inverse();
		Double1D vecW = ((Double2D) hypers.get("W")).colVec();
		Double2D Sⁿ_kr_Ih = Sⁿ.kron(Double2D.ident(h));
		Double1D Sⁿ_kr_Ih_vW = Sⁿ_kr_Ih.mult(vecW);
		extended_hypers.put("Sⁿ_kr_Ih", Sⁿ_kr_Ih);
		extended_hypers.put("Sⁿ_kr_Ih_vW", Sⁿ_kr_Ih_vW);
		
		int p = Y.numCols();
		extended_hypers.put("p", p);
		extended_hypers.put("lΓκ", lΓ(κ.value(), p));
		extended_hypers.put("lΓκ_1", lΓ(κ.value()+1, p));
		
		double lπ = 0.5*p*Math.log(Math.PI);
		double l2π = 0.5*p*Math.log(2*Math.PI);
		extended_hypers.put("lπ", lπ);
		extended_hypers.put("l2π", l2π);
	
		extended_hypers.put("Y", Y);
		extended_hypers.put("X", X);
	
		int N = Y.numRows();
		extended_hypers.put("N", N);
		
		Double2D[] xxʹ = new Double2D[N];
		Double2D[] yyʹ = new Double2D[N];
		Double2D[] xyʹ = new Double2D[N];

		extended_hypers.put("xxʹ", xxʹ);
		extended_hypers.put("xyʹ", xyʹ);
		extended_hypers.put("yyʹ", yyʹ);
		
		extended_hypers.put("alpha_a", (Double0D) hypers.get("alpha_a"));
		extended_hypers.put("alpha_b", (Double0D) hypers.get("alpha_b"));
		
		return extended_hypers;
	}	

	class PostOmegaDist extends ProbDistMC<Double2D> {
		// Fixed parameters
		Double2D Φ = (Double2D) MLM_sample_params.this.extended_hypers.get("Φ");
		Double0D λ = (Double0D) MLM_sample_params.this.extended_hypers.get("λ");
		int p = (Integer) MLM_sample_params.this.extended_hypers.get("p");
		int N = (Integer) MLM_sample_params.this.extended_hypers.get("N");

		// Chain parameters
		Integer1D R;
		Double3D Σ;
		Double3D A;
		Double2D A0;

		private InverseWishartDist iwDist = new InverseWishartDist();
		
		public PostOmegaDist() { }

		@Override
		public void setMCState(ChainLink l) {
			this.R = ((Integer1D) l.get("z").getNumericValue()).items();
			this.Σ = (Double3D) l.get("Sigma").getNumericValue();
			this.A = (Double3D) l.get("A").getNumericValue();
			this.A0 = (Double2D) l.get("A0").getNumericValue();
			this.initialized = true;
		}

		@Override
		protected Double2D genVariate() {
			Double2D Φᵤ = Φ.copy(); 
			int rj = 0;
			for (rj=0; rj<R.size(); rj++) {
				int r = R.value()[rj];
				Double2D A_A0 = A.get(r).minus(A0);
				Φᵤ.plusEquals(A_A0.transposeMult(Σ.get(r).inverse()).mult(A_A0));
			}
			
			iwDist.setParms(Φᵤ, λ.plus(p*N));
			return(iwDist.variate());
		}
	}

	class PostSigmaDist extends ProbDistMC<Double3D> {
		// Fixed parameters
		private Double0D κ = (Double0D) MLM_sample_params.this.extended_hypers.get("κ");
		private int p = (Integer) MLM_sample_params.this.extended_hypers.get("p");
		private Double2D Ψ = (Double2D) MLM_sample_params.this.extended_hypers.get("Ψ");
		private Double2D Y = (Double2D) MLM_sample_params.this.extended_hypers.get("Y");
		private Double2D X = (Double2D) MLM_sample_params.this.extended_hypers.get("X");

		// Chain parameters
		private Double2D Ωⁿ;
		private Double2D A0;
		private Double2D ΩⁿA0;
		private Double2D Ψ_A0ʹΩⁿA0;
		private Integer1D z;
		private Integer1D R;
		private int max_r;
		
		private InverseWishartDist iwDist = new InverseWishartDist();

		public PostSigmaDist() { }

		@Override
		public void setMCState(ChainLink l) {
			this.Ωⁿ = ((Double2D) l.get("Omega").getNumericValue()).inverse();
			this.A0 = (Double2D) l.get("A0").getNumericValue();
			this.ΩⁿA0 = Ωⁿ.mult(A0);
			this.Ψ_A0ʹΩⁿA0 = Ψ.plus(A0.transposeMult(ΩⁿA0)); 
			this.z = (Integer1D) l.get("z").getNumericValue();
			this.R = z.items();
			this.max_r = R.value()[R.size()-1];

			this.initialized = true;
		}

		@Override
		protected Double3D genVariate() {
			DenseDoubleMatrix3D Σ = new DenseDoubleMatrix3D(max_r+1,p,p);
			int rj = 0;
			for (rj=0; rj<R.size(); rj++) {
				int r = R.value()[rj];
				int[] zᵣ = z.which(r).value();
				int Nᵣ = zᵣ.length; 
				
				Double2D Xᵣ = X.getAll(zᵣ);
				Double2D Yᵣ = Y.getAll(zᵣ);
				
				Double2D XᵣʹXᵣ = Xᵣ.transposeMult(Xᵣ); // FIXME: here - extended chain
				Double2D XᵣʹYᵣ = Xᵣ.transposeMult(Yᵣ); // FIXME: here - extended chain
				Double2D YᵣʹYᵣ = Yᵣ.transposeMult(Yᵣ); // FIXME: here - extended chain
				
				Double2D Ωⁿ_XᵣʹXᵣ = Ωⁿ.plus(XᵣʹXᵣ);
				Double2D ΩⁿA0_XᵣʹYᵣ = ΩⁿA0.plus(XᵣʹYᵣ);
				Double2D Ψᵣᵤ = Ψ_A0ʹΩⁿA0.plus(YᵣʹYᵣ).minus(ΩⁿA0_XᵣʹYᵣ.transposeMult(Ωⁿ_XᵣʹXᵣ.inverse().mult(ΩⁿA0_XᵣʹYᵣ)));
				
				iwDist.setParms(Ψᵣᵤ, κ.plus(Nᵣ));
				Double2D Σᵣ = iwDist.variate();
				
				int slice = r; // FIXME: this is what you would change for "compression"
				for (int row=0; row<p; row++) {
					for (int col=0; col<p; col++) {
						Σ.setQuick(slice, row, col, Σᵣ.getDm().getQuick(row, col));
					}
				}
			}
								
			return new Double3D(Σ);			
		}
	}

	class PostADist extends ProbDistMC<Double3D> {
		// Fixed parameters
		private int p = (Integer) MLM_sample_params.this.extended_hypers.get("p");
		private int h = (Integer) MLM_sample_params.this.extended_hypers.get("h");
		private Double2D Y = ((Double2D) MLM_sample_params.this.extended_hypers.get("Y"));
		private Double2D X = ((Double2D) MLM_sample_params.this.extended_hypers.get("X"));

		// Chain parameters
		private Integer1D z;
		private Double2D Ωⁿ;
		private Double3D Σ;
		private Double2D A0;
		private Double2D ΩⁿA0;
		private Integer1D R;
		private int max_r;
		
		private MatrixNormalDist mnDist = new MatrixNormalDist();

		public PostADist() { }

		@Override
		public void setMCState(ChainLink l) {
			this.A0 = (Double2D) l.get("A0").getNumericValue();
			this.Ωⁿ = ((Double2D) l.get("Omega").getNumericValue()).inverse();
			this.ΩⁿA0 = Ωⁿ.mult(A0);
			this.Σ = (Double3D) l.get("Sigma").getNumericValue();
			this.z = (Integer1D) l.get("z").getNumericValue();
			this.R = z.items();
			this.max_r = R.value()[R.size()-1];
			this.initialized = true;
		}

		@Override
		protected Double3D genVariate() {
			// FIXME - Remember to add the compressèdness and extra layer of crud you need to
			// support this compressèdness
			DenseDoubleMatrix3D A = new DenseDoubleMatrix3D(max_r+1,h,p);
			int rj = 0;
			for (rj=0; rj<R.size(); rj++) {
				int r = R.value()[rj];
				int[] zᵣ = z.which(r).value();
				Double2D Xᵣ = X.getAll(zᵣ);
				Double2D Yᵣ = Y.getAll(zᵣ);
				
				Double2D XᵣʹXᵣ = Xᵣ.transposeMult(Xᵣ); // FIXME: here - extended chain
				Double2D Ωᵣᵤ = Ωⁿ.plus(XᵣʹXᵣ).inverse(); 
				
				Double2D XᵣʹYᵣ = Xᵣ.transposeMult(Yᵣ); // FIXME: here - extended chain
				Double2D Α0ᵣᵤ = Ωᵣᵤ.mult(XᵣʹYᵣ.plus(ΩⁿA0));
				
				mnDist.setParms(Α0ᵣᵤ, Σ.get(r), Ωᵣᵤ);
				Double2D Aᵣ = mnDist.variate();
				
				int slice = r; // FIXME: this is what you would change for "compression"
				for (int row=0; row<h; row++) {
					for (int col=0; col<p; col++) {
						A.setQuick(slice, row, col, Aᵣ.getDm().getQuick(row, col));
					}
				}
			}
								
			return new Double3D(A);
		}
	}

	class CRPDist extends ProbDistMC<Integer1D> {
		// Fixed parameters
		private Double0D κ_1 = ((Double0D) MLM_sample_params.this.extended_hypers.get("κ_1"));
		private Double2D Ψ = ((Double2D) MLM_sample_params.this.extended_hypers.get("Ψ"));
		
		double ldet_Ψ = (Double) (MLM_sample_params.this.extended_hypers.get("ldet_Ψ"));
		private double lπ = (Double) MLM_sample_params.this.extended_hypers.get("lπ");
		private double l2π = (Double) MLM_sample_params.this.extended_hypers.get("l2π");
		private int N = (Integer) MLM_sample_params.this.extended_hypers.get("N");
		private int p = (Integer) MLM_sample_params.this.extended_hypers.get("p");
		private Double2D X = ((Double2D) MLM_sample_params.this.extended_hypers.get("X"));
		private Double2D Y = ((Double2D) MLM_sample_params.this.extended_hypers.get("Y"));
		private Double2D[] xxʹ = (Double2D[]) MLM_sample_params.this.extended_hypers.get("xxʹ");
		private Double2D[] xyʹ = (Double2D[]) MLM_sample_params.this.extended_hypers.get("xyʹ");
		private Double2D[] yyʹ = (Double2D[]) MLM_sample_params.this.extended_hypers.get("yyʹ");
		private double lΓκ = (Double) MLM_sample_params.this.extended_hypers.get("lΓκ");
		private double lΓκ_1 = (Double) MLM_sample_params.this.extended_hypers.get("lΓκ_1");
		
		private Map<Integer,Integer> Nᵣ = new HashMap<Integer,Integer>();
		
		// Chain parameters
		private Double2D A0;
		private Double2D Ωⁿ;
		private Double2D ΩⁿA0;
		private Double2D Ψ_A0ʹΩⁿA0;
		private Double3D Σ;
		private Double3D A;
		private Integer1D z;
		private Double0D α;
		private double[] prior_lil = new double[N];

		private Double1D postProb;
		private CategoricalDist catDist = new CategoricalDist();
		private InverseWishartDist iwDist = new InverseWishartDist();
		private MatrixNormalDist mnDist = new MatrixNormalDist();

		public CRPDist() { }

		@Override
		public void setMCState(ChainLink l) {
			this.A0 = (Double2D) l.get("A0").getNumericValue();
			this.Ωⁿ = ((Double2D) l.get("Omega").getNumericValue()).inverse();
			this.ΩⁿA0 = Ωⁿ.mult(A0);
			this.Ψ_A0ʹΩⁿA0 = Ψ.plus(A0.transposeMult(ΩⁿA0)); 
			this.Σ = (Double3D) l.get("Sigma").getNumericValue();
			this.A = (Double3D) l.get("A").getNumericValue();
			this.z = (Integer1D) l.get("z").getNumericValue();
			this.α = (Double0D) l.get("alpha").getNumericValue();
			
			double ldet_Ωⁿ = 0.5*p*Math.log(Ωⁿ.det());	
			for (int i=0; i<N; i++) {
				Double1D xᵢ = X.get(i);
				Double1D yᵢ = Y.get(i);
				
				xxʹ[i] = xᵢ.outer(xᵢ);
				xyʹ[i] = xᵢ.outer(yᵢ);
				yyʹ[i] = yᵢ.outer(yᵢ);
				
				Double2D Ωᵤ = Ωⁿ.plus(xxʹ[i]).inverse();
				Double2D ΩⁿA0_xᵢyᵢʹ = ΩⁿA0.plus(xyʹ[i]);
				Double2D Ψᵤ = Ψ_A0ʹΩⁿA0.plus(yyʹ[i]).minus(ΩⁿA0_xᵢyᵢʹ.transposeMult(Ωᵤ).mult(ΩⁿA0_xᵢyᵢʹ));
				
				double ldet_Ωᵤ = 0.5*p*Math.log(Ωᵤ.det());
				double ldet_Ψᵤ = 0.5*κ_1.value()*Math.log(Ψᵤ.det());				
				
				prior_lil[i] = lΓκ_1 + ldet_Ψ + ldet_Ωⁿ + ldet_Ωᵤ - lπ - lΓκ - ldet_Ψᵤ;
			}		
			
			this.initialized = true;
		}

		protected Integer0D catVariate() {
			this.catDist.setParms(this.postProb);
			return this.catDist.variate();
		}


		@Override
		protected Integer1D genVariate() {
			double lα = Math.log(α.value());
			Integer1D z = new Integer1D(this.z.value().clone());
			
			for (int i = 0; i < N; i++) {
				Integer1D R = z.items();
				int zi_old_value = z.value()[i];
				z.set(i,-1);
				
				int Nᵣᵢ;
				if (Nᵣ.containsKey(zi_old_value)) {
					Nᵣᵢ = Nᵣ.get(zi_old_value)-1;
				} else {
					int[] zᵣ = z.which(zi_old_value).value();
					Nᵣᵢ = zᵣ.length;
					if (Nᵣᵢ > 0) {
						Nᵣ.put(zi_old_value,Nᵣᵢ+1);
					} else {
						Nᵣ.put(zi_old_value,1);
					}
				}
				
				// Existing categories
				double[] logP = new double[R.size() + 1];	
				
				double maxLogP = -Double.MAX_VALUE;
				int rj = 0;
				for (rj=0; rj<R.size(); rj++) {
					int r = R.value()[rj];
					
					Double1D Aᵣʹxᵢ = A.get(r).transposeMult(X.get(i));
					Double1D yᵢ_Aᵣʹxᵢ = Y.get(i).minus(Aᵣʹxᵢ);
					double llik_e_term = -0.5*yᵢ_Aᵣʹxᵢ.mult(Σ.get(r).inverse().mult(yᵢ_Aᵣʹxᵢ));
					double llik_nc_term = -l2π - 0.5*Math.log(Σ.get(r).det());
		
					if (r == zi_old_value) {
						if (Nᵣᵢ == 0) {
							logP[rj] = -Double.MAX_VALUE;
						} else {
							logP[rj] = Math.log(Nᵣᵢ) + llik_nc_term + llik_e_term;
						}
					} else {
						if (!Nᵣ.containsKey(r)) {
							int[] zᵣ = z.which(r).value();
							Nᵣ.put(r, zᵣ.length);
						}
						logP[rj] = Math.log(Nᵣ.get(r)) + llik_nc_term + llik_e_term;
					}
					if (logP[rj] == Double.POSITIVE_INFINITY) {
						logP[rj] = Double.MAX_VALUE;
					}				
					if (logP[rj] > maxLogP) {
						maxLogP = logP[rj];
					}
				}
				// New category
				logP[rj] = lα + prior_lil[i];
				
				// Select a category for this point
				int newc_rjindex = rj; // last index
				this.postProb = (new Double1D(logP)).minus(maxLogP).exp();
				int zi_new_rjindex = this.catVariate().value();
				int zi_new_value;
				if (zi_new_rjindex == newc_rjindex) {
					int maxActive = R.value()[R.size()-1];
					zi_new_value = z.minNotIn(0, maxActive);
				} else {
					zi_new_value = R.value()[zi_new_rjindex];
				}
				
				// Sample new category if necessary
				if (zi_new_rjindex == newc_rjindex) {
					Double2D ΩⁿA0_xᵢyᵢʹ = ΩⁿA0.plus(xyʹ[i]);
					Double2D Ωᵤ = Ωⁿ.plus(xxʹ[i]).inverse();
					Double2D Ωᵤ_ΩⁿA0_xᵢyᵢʹ = Ωᵤ.mult(ΩⁿA0_xᵢyᵢʹ);
					Double2D Ψᵤ = Ψ_A0ʹΩⁿA0.plus(yyʹ[i]).minus(ΩⁿA0_xᵢyᵢʹ.transposeMult(Ωᵤ_ΩⁿA0_xᵢyᵢʹ));
					iwDist.setParms(Ψᵤ, κ_1);
					Double2D Σᵢ = iwDist.variate();
					mnDist.setParms(Ωᵤ_ΩⁿA0_xᵢyᵢʹ, Σᵢ, Ωᵤ);
					Double2D Aᵢ = mnDist.variate();
					A.set(zi_new_value, Aᵢ);
					Σ.set(zi_new_value, Σᵢ);
				}
				
				// Update counts
				if (zi_new_value != zi_old_value) {
					if (zi_new_rjindex == newc_rjindex) {
						Nᵣ.put(zi_new_value, 1);
					} else {
						Nᵣ.put(zi_new_value,Nᵣ.get(zi_new_value)+1);
					}
					if (Nᵣᵢ == 0) {
						Nᵣ.remove(zi_old_value);
					} else {
						Nᵣ.put(zi_old_value,Nᵣᵢ);
					}
				}
				z.set(i, zi_new_value);
			}
			
			return z;	
		}
	}

	class PostA0Dist extends ProbDistMC<Double2D> {
		// Fixed parameters
		private Double2D Sⁿ_kr_Ih = (Double2D) MLM_sample_params.this.extended_hypers.get("Sⁿ_kr_Ih");
		private Double1D Sⁿ_kr_Ih_vW = (Double1D) MLM_sample_params.this.extended_hypers.get("Sⁿ_kr_Ih_vW");
		private int h = (Integer) MLM_sample_params.this.extended_hypers.get("h");
		private MVNormalDist mvnDist = new MVNormalDist();

		// Chain parameters
		private Integer1D z;
		private Double2D Ωⁿ;
		private Double3D Σ;
		private Double3D A;

		public PostA0Dist() { }

		@Override
		public void setMCState(ChainLink l) {
			this.z = (Integer1D) l.get("z").getNumericValue();
			this.Ωⁿ = ((Double2D) l.get("Omega").getNumericValue()).inverse();
			this.Σ = (Double3D) l.get("Sigma").getNumericValue();
			this.A = (Double3D) l.get("A").getNumericValue();
			
			int[] R = z.items().value();
			Double2D precision = Sⁿ_kr_Ih.copy();
			Double1D unscaled_centre = Sⁿ_kr_Ih_vW.copy();
			for (int r : R) {
				Double2D Σⁿ = Σ.get(r).inverse();
				Double2D Σⁿ_kr_Ωⁿ = Σⁿ.kron(Ωⁿ);
				precision.plusEquals(Σⁿ_kr_Ωⁿ);
				unscaled_centre.plusEquals(Σⁿ_kr_Ωⁿ.mult(A.get(r).colVec()));
			}
			Double2D covariance = precision.inverse();
			mvnDist.setParms(covariance.mult(unscaled_centre), covariance);	
			
			this.initialized = true;
		}

		@Override
		protected Double2D genVariate() {
			return this.mvnDist.variate().toDouble2D(h);
		}
	}

	class PostAlphaDist extends ProbDistMC<Double0D> {
		private BetaDist betaDist = new BetaDist();
		private GammaDist gammaDist = new GammaDist();
		private CategoricalDist catDist = new CategoricalDist();

		// Fixed parameters
		private Double0D alpha_a = (Double0D) MLM_sample_params.this.extended_hypers.get("alpha_a");
		private Double0D alpha_b = (Double0D) MLM_sample_params.this.extended_hypers.get("alpha_b");
		private Double0D N = new Double0D((Integer) MLM_sample_params.this.extended_hypers.get("N")); 

		// Chain parameters
		private Double0D α;
		private Integer1D z;

		public PostAlphaDist() { }

		@Override
		public void setMCState(ChainLink l) {
			this.α = (Double0D) l.get("alpha").getNumericValue();
			this.z = (Integer1D) l.get("z").getNumericValue();
			this.initialized = true;
		}

		@Override
		protected Double0D genVariate() {
			betaDist.setParms(α.plus(1), N);	
			double x = betaDist.variate().value();
			double log_x = Math.log(x);
		
			int K = z.items().size();
			
			double p1 = alpha_a.value() + K - 1;
			double p2 = N.value() * (alpha_b.value() - log_x);
			catDist.setParms(new Double1D(p1, p2));
			int i = catDist.variate().value();

			Double0D post_a;
			if (i == 0) {
				post_a = alpha_a.plus(K);
			} else {
				post_a = alpha_a.plus(K-1);
			}
			this.gammaDist.setParms(post_a, alpha_b.plus(-log_x));			
			return gammaDist.variate();
		}
	}	
	
	@SuppressWarnings("unchecked")
	public static ChainLink pointEstimate(Chain c, Map<String, Object> extended_hypers) {
		double max = Double.NEGATIVE_INFINITY;
		int argmax = c.size() - 1;
		
		GammaDist gammaDist = new GammaDist((Double0D) extended_hypers.get("alpha_a"), (Double0D) extended_hypers.get("alpha_b"));	
		InverseWishartDist iwDist_Ω = new InverseWishartDist((Double2D) extended_hypers.get("Φ"), (Double0D) extended_hypers.get("lambda"));	
		int h = (Integer) extended_hypers.get("h");
		MatrixNormalDist mnDist = new MatrixNormalDist((Double2D) extended_hypers.get("W"), (Double2D) extended_hypers.get("S"), Double2D.ident(h));	
		
		for (int m = 0; m < c.size(); m++) {
			double posterior = 0;
			ChainLink cl_m = c.get(m);
			
			// α
			Double0D α = ((RandomVar<Double0D>) cl_m.get("alpha")).getNumericValue();
			posterior += gammaDist.logDensity(α);
						
			// Omega
			Double2D Ω = ((RandomVar<Double2D>) cl_m.get("Omega")).getNumericValue();
			posterior += iwDist_Ω.logDensity(Ω);
			
			// A0
			Double2D A0 = ((RandomVar<Double2D>) cl_m.get("A0")).getNumericValue();
			posterior += mnDist.logDensity(A0);
	
			// Model parameters and likelihood
//			Double3D Σ = ((RandomVar<Double3D>) cl_m.get("Σ")).getNumericValue();
//			Double3D A = ((RandomVar<Double3D>) cl_m.get("A")).getNumericValue();
//			double lα = Math.log(α.value());
//			int[] z = ((RandomVar<Integer1D>) cl_m.get("z")).getNumericValue().value();
//			Map<Integer,Object> cat_changed = new HashMap<Integer,Object>();
//			for (int i = 0; i < N; i++) {
////				if (m == 0) {
////					log_cumul_lil[i] = new HashMap<Integer,Double>();
////					log_cumul_lil[i] = new HashMap<Integer,Double>();
////				}
////				
//				int rᵢ = z[i];
//				
//				if (i > 0) {
//					int[] zᵢ_ints = new int[i];
//					System.arraycopy(z,0,zᵢ_ints,0,i); // FIXME - yeah, this is dumb - fix it later
//					Integer1D zᵢ = new Integer1D(zᵢ_ints);
//					Integer1D Rᵢ = zᵢ.items();
//					int[] zᵣᵢ = zᵢ.which(rᵢ).value();
//					int Nᵣᵢ = zᵣᵢ.length;
//		
//					if (cat_changed.containsKey(rᵢ)) {
//						// Category has changed since last sample
//						// The likelihood integral for this point needs
//						// to be recomputed
//						if (Nᵣᵢ > 0) {
//							// Category still has points before i
//							Double2D Xᵣᵢ = X.getAll(zᵣᵢ);
//							Double2D Yᵣᵢ = Y.getAll(zᵣᵢ);
//							Double2D XᵣᵢʹXᵣᵢ = Xᵣᵢ.transposeMult(Xᵣᵢ);
//							Double2D YᵣᵢʹYᵣᵢ = Yᵣᵢ.transposeMult(Yᵣᵢ); 
//							Double2D XᵣᵢʹYᵣᵢ = Xᵣᵢ.transposeMult(Yᵣᵢ);
//							
//							double κ_Nᵣᵢ = κ+Nᵣᵢ;
//							double κ_Nᵣᵢ_1 = κ_Nᵣᵢ+1;
//							
//							Double2D Ωⁿ_XᵣᵢʹXᵣᵢ = Ωⁿ.plus(XᵣᵢʹXᵣᵢ);
//							Double2D Ψ_A0ʹΩⁿA0_YᵣᵢʹYᵣᵢ = Ψ_A0ʹΩⁿA0.plus(YᵣᵢʹYᵣᵢ);
//							Double2D ΩⁿA0_XᵣᵢʹYᵣᵢ = ΩⁿA0.plus(XᵣᵢʹYᵣᵢ);
//							Double2D ΩⁿA0_XᵣᵢʹYᵣᵢ_xᵢyᵢʹ = ΩⁿA0_XᵣᵢʹYᵣᵢ.plus(xyʹ[i]);
//							
//							Double2D Ωᵥ = Ωⁿ_XᵣᵢʹXᵣᵢ.inverse();
//							Double2D Ωᵤᵥ = Ωⁿ_XᵣᵢʹXᵣᵢ.plus(xxʹ[i]).inverse();
//							Double2D Ψᵥ = Ψ_A0ʹΩⁿA0_YᵣᵢʹYᵣᵢ.minus(ΩⁿA0_XᵣᵢʹYᵣᵢ.transposeMult(Ωᵥ.mult(ΩⁿA0_XᵣᵢʹYᵣᵢ)));
//							Double2D Ψᵤᵥ = Ψ_A0ʹΩⁿA0_YᵣᵢʹYᵣᵢ.plus(yyʹ[i]).minus(ΩⁿA0_XᵣᵢʹYᵣᵢ_xᵢyᵢʹ.transposeMult(Ωᵤᵥ.mult(ΩⁿA0_XᵣᵢʹYᵣᵢ_xᵢyᵢʹ)));
//							
//							double ldet_Ωᵥ = 0.5*p*Math.log(Ωᵥ.det());
//							double ldet_Ωᵤᵥ = 0.5*p*Math.log(Ωᵤᵥ.det());
//							double ldet_Ψᵥ = 0.5*κ_Nᵣᵢ*Math.log(Ψᵥ.det());
//							double ldet_Ψᵤᵥ = 0.5*κ_Nᵣᵢ_1*Math.log(Ψᵤᵥ.det());	
//	
//							double lΓκ_Nᵣᵢ = lΓ((κ_Nᵣᵢ)/2, p);
//							double lΓκ_Nᵣᵢ_1 = lΓ((κ_Nᵣᵢ_1)/2, p);
//							
//							double lil = Math.log(Nᵣᵢ) + lΓκ_Nᵣᵢ_1 + ldet_Ωᵤᵥ + ldet_Ψᵥ - lπ - lΓκ_Nᵣᵢ - ldet_Ωᵥ - ldet_Ψᵤᵥ;
//							
//							log_cumul_lil[i].put(rᵢ, lil);
//							cumul_lil[i].put(rᵢ, Math.exp(lil));						
//						} else {
//							// Category now has no points before i: use the "new category" term
//							double lil = lα + prior_lil[i];
//							log_cumul_lil[i].put(rᵢ, lil);
//							cumul_lil[i].put(rᵢ, Math.exp(lil));
//						}
//					}					
//					// Compute normalizing constant
//					double nc = cumul_lil[i].get(rᵢ);
//					for (int j_r=0; j_r<Rᵢ.size(); j_r++) {
//						int r = Rᵢ.value()[j_r];
//						if (r != rᵢ) {
//							if (cat_changed.containsKey(r)) {
//								// Another category has changed: recompute posterior integrated likelihood
//								int[] zᵣ = zᵢ.which(r).value();
//								int Nᵣ = zᵣ.length;
//					
//								Double2D Xᵣ = X.getAll(zᵣ);
//								Double2D Yᵣ = Y.getAll(zᵣ);
//								Double2D XᵣʹXᵣ = Xᵣ.transposeMult(Xᵣ);
//								Double2D YᵣʹYᵣ = Yᵣ.transposeMult(Yᵣ); 
//								Double2D XᵣʹYᵣ = Xᵣ.transposeMult(Yᵣ); 
//								
//								double κ_Nᵣ = κ+Nᵣ;
//								double κ_Nᵣ_1 = κ_Nᵣ+1;
//								
//								Double2D Ωⁿ_XᵣʹXᵣ = Ωⁿ.plus(XᵣʹXᵣ);
//								Double2D Ψ_A0ʹΩⁿA0_YᵣʹYᵣ = Ψ_A0ʹΩⁿA0.plus(YᵣʹYᵣ);
//								Double2D ΩⁿA0_XᵣʹYᵣ = ΩⁿA0.plus(XᵣʹYᵣ);
//								Double2D ΩⁿA0_XᵣʹYᵣ_xᵢyᵢʹ = ΩⁿA0_XᵣʹYᵣ.plus(xyʹ[i]);
//								
//								Double2D Ωᵥ = Ωⁿ_XᵣʹXᵣ.inverse();
//								Double2D Ωᵤᵥ = Ωⁿ_XᵣʹXᵣ.plus(xxʹ[i]).inverse();
//								Double2D Ψᵥ = Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.minus(ΩⁿA0_XᵣʹYᵣ.transposeMult(Ωᵥ.mult(ΩⁿA0_XᵣʹYᵣ)));
//								Double2D Ψᵤᵥ = Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.plus(yyʹ[i]).minus(ΩⁿA0_XᵣʹYᵣ_xᵢyᵢʹ.transposeMult(Ωᵤᵥ.mult(ΩⁿA0_XᵣʹYᵣ_xᵢyᵢʹ)));
//								
//								double ldet_Ωᵥ = 0.5*p*Math.log(Ωᵥ.det());
//								double ldet_Ωᵤᵥ = 0.5*p*Math.log(Ωᵤᵥ.det());
//								double ldet_Ψᵥ = 0.5*κ_Nᵣ*Math.log(Ψᵥ.det());
//								double ldet_Ψᵤᵥ = 0.5*κ_Nᵣ_1*Math.log(Ψᵤᵥ.det());	
//		
//								double lΓκ_Nᵣ = lΓ((κ_Nᵣ)/2, p);
//								double lΓκ_Nᵣ_1 = lΓ((κ_Nᵣ_1)/2, p);
//								
//								double lil = Math.log(Nᵣ) + lΓκ_Nᵣ_1 + ldet_Ωᵤᵥ + ldet_Ψᵥ - lπ - lΓκ_Nᵣ - ldet_Ωᵥ - ldet_Ψᵤᵥ;
//								log_cumul_lil[i].put(r, lil);
//								cumul_lil[i].put(r, Math.exp(lil));							
//							}
//							nc += cumul_lil[i].get(r);
//						}
//					}
//					if (Nᵣᵢ > 0) {
//						nc += Math.exp(lα + prior_lil[i]);
//					}
//					posterior += log_cumul_lil[i].get(rᵢ) - Math.log(nc);
//				} else { // i == 0: probability 1
//					posterior += 0;
//				}
//
//				if (m == 0) {
//					cat_changed.put(rᵢ,null);
//				} else {
//					if (!cat_changed.containsKey(rᵢ)) {
//						// Up to i-1, category had not changed since last sample...
//						if (rᵢ != previous_z[i]) {
//							// ... but category (and old counterpart) just changed: notify future points
//							cat_changed.put(rᵢ,null);
//							cat_changed.put(previous_z[i],null);
//						}
//					}					
//				}
//			}				
			
			// Set max
			if (posterior > max) {
				argmax = m;
				max = posterior;
			}
			
		}

		return c.get(argmax);			
	}

	public MLM_sample_params() {
		super();
	}
	
	public MLM_sample_params(Map<String, Flattenable> hypers, Map<String, Flattenable> init, Double2D data) {
		super(hypers);
		
		// Set up parameters
		this.params = new HashMap<String, RandomVar<? extends Flattenable>>();
		
		// Extend hyperparameters
		this.extended_hypers = MLM_sample_params.extend_hypers(hypers, data, (Double2D) init.get("X"));	
		
		// Omega
		this.params.put("Omega", new RandomVar<Double2D>("Omega", new PostOmegaDist(), (Double2D) init.get("Omega")));

		// Sg
		this.params.put("Sigma", new RandomVar<Double3D>("Sigma", new PostSigmaDist(), (Double3D) init.get("Sigma")));

		// A
		this.params.put("A", new RandomVar<Double3D>("A", new PostADist(), (Double3D) init.get("A")));

		// z
		this.params.put("z", new RandomVar<Integer1D>("z", new CRPDist(), (Integer1D) init.get("z")));

		// A0
		this.params.put("A0", new RandomVar<Double2D>("A0", new PostA0Dist(), (Double2D) init.get("A0")));

		// alpha
		this.params.put("alpha", new RandomVar<Double0D>("alpha", new PostAlphaDist(), (Double0D) init.get("alpha")));
	}

	@Override
	public ChainLink getInitialLink() {
		ChainLink cl = new ChainLink();
		cl.add(this.params.get("z"));
		cl.add(this.params.get("A"));
		cl.add(this.params.get("Sigma"));
		cl.add(this.params.get("Omega"));
		cl.add(this.params.get("A0"));
		cl.add(this.params.get("alpha"));
		return cl;
	}
}
