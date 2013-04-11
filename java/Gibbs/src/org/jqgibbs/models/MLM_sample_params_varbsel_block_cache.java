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

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix3D;
import cern.jet.stat.Gamma;

public class MLM_sample_params_varbsel_block_cache extends Model {

	private static double lΓ(double x, int p) {
		double r = 0.25*p*(p-1)*Math.log(Math.PI);
		for (int i=0; i<p; i++) {
			r += Gamma.logGamma(x-i/2.0);
		}
		return r;
	}
	
	public static Map<String,Object> extend_hypers(Map<String,Flattenable> hypers, Double2D Y, Double2D X) {
		Map<String,Object> extended_hypers = new HashMap<String,Object>();
		
		Integer z_lag = ((Integer0D) hypers.get("zlag")).value();
		extended_hypers.put("z_lag", z_lag);
		
		Double τ = ((Double0D) hypers.get("tau")).value();
		extended_hypers.put("τ", τ);
		extended_hypers.put("lτ", Math.log(τ));
		extended_hypers.put("lτʹ", Math.log(1-τ));
		
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
		Double2D W = (Double2D) hypers.get("W");
		extended_hypers.put("Sⁿ", Sⁿ);
		extended_hypers.put("W", W);
		extended_hypers.put("ldet_Sⁿ", Math.log(Sⁿ.det()));
		
		int p = Y.numCols();
		extended_hypers.put("p", p);
		extended_hypers.put("lΓκ", lΓ(κ.value()/2, p));
		extended_hypers.put("lΓκ_1", lΓ((κ.value()+1)/2, p));
		
		double lπ = 0.5*p*Math.log(Math.PI);
		double l2π = 0.5*p*Math.log(2*Math.PI);
		extended_hypers.put("lπ", lπ);
		extended_hypers.put("l2π", l2π);
	
		extended_hypers.put("Y", Y);
		extended_hypers.put("X", X);
	
		int N = Y.numRows();
		extended_hypers.put("N", N);
		

		
		Map<Integer1D,Double2D> X_γs = new HashMap<Integer1D,Double2D>();
		Map<Integer1D,Double2D[]> xxʹ_γs = new HashMap<Integer1D,Double2D[]>();
		Map<Integer1D,Double2D[]> xyʹ_γs = new HashMap<Integer1D,Double2D[]>();
		
		int num_γs = (int) Math.pow(2, h-1);
		Integer1D[] γs = new Integer1D[num_γs];
		for (int γ_sum=0; γ_sum<num_γs; γ_sum++) {
			int[] bits = new int[h];
			bits[0] = 1;
			if (γ_sum>0) {
				int max_pos = (int) Math.floor(Math.log(γ_sum)/Math.log(2));
				for (int i=0; i<max_pos+1; i++) {
					bits[(h-1)-i] = ((γ_sum & (1 << i)) != 0) ? 1 : 0;
				}
			}
			γs[γ_sum] = new Integer1D(bits);
		}
		
		for (Integer1D γ : γs) {
			int[] γ_indices = γ.which(1).value();
			Double2D X_γ = X.subMatrix(null, γ_indices);
			X_γs.put(γ, X_γ);
			Double2D[] xxʹ_γ = new Double2D[N];
			Double2D[] xyʹ_γ = new Double2D[N];
			for (int i=0; i<N; i++) {
				Double1D x_γ = X_γ.get(i);
				Double1D y = Y.get(i);
				xxʹ_γ[i] = x_γ.outer(x_γ);
				xyʹ_γ[i] = x_γ.outer(y);
			}
			xxʹ_γs.put(γ, xxʹ_γ);
			xyʹ_γs.put(γ, xyʹ_γ);
		}
		extended_hypers.put("X_γs", X_γs);
		extended_hypers.put("xxʹ_γs", xxʹ_γs);
		extended_hypers.put("xyʹ_γs", xyʹ_γs);
		
		Double2D[] yyʹ = new Double2D[N];
		for (int i=0; i<N; i++) {
			Double1D y = Y.get(i);
			yyʹ[i] = y.outer(y);
		}
		extended_hypers.put("yyʹ", yyʹ);
		
		extended_hypers.put("alpha_a", hypers.get("alpha_a"));
		extended_hypers.put("alpha_b", hypers.get("alpha_b"));
		
		return extended_hypers;
	}	
	
	class TemperatureSchedule extends ProbDistMC<Double0D> {
		private double lT0;
		private double lTf;
		private int deadline;
		
		private int t = 0;
		
		TemperatureSchedule(double T0, double Tf, int deadline) {
			this.lT0 = Math.log(T0);
			this.lTf = Math.log(Tf);
			this.deadline = deadline;
		}
		
		@Override
		public void setMCState(ChainLink l) {};
		
		@Override
		protected Double0D genVariate() {
			double lT;
			if (t < deadline) {
				double scale = (++t)/deadline;
				lT = (1-scale)*lT0 + scale*lTf;
			} else {
				lT = lTf;
			}
			return new Double0D(Math.exp(lT));
		}
	}
	
	class PostOmegaDist extends ProbDistMC<Double2D> {
		// Fixed parameters
		Double2D Φ = (Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("Φ");
		Double0D λ = (Double0D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("λ");
		int p = (Integer) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("p");
		int N = (Integer) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("N");

		// Chain parameters
		Integer1D γ;
		Double2D Φ_γ;
		Double0D λ_γ;
		int p_γ;
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
				Φᵤ.plusEquals(A_A0.mult(Σ.get(r).inverse()).multTranspose(A_A0));
			}
			
			iwDist.setParms(Φᵤ, λ.plus(p*R.size()));
			return(iwDist.variate());
		}
	}

	class PostSigmaDist extends ProbDistMC<Double3D> {
		// Fixed parameters
		private Double0D κ = (Double0D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("κ");
		private int p = (Integer) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("p");
		private Double2D Ψ = (Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("Ψ");
		private Double2D Y = (Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("Y");
		private Double2D X = (Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("X");

		// Chain parameters
		private Double2D X_γ;
		private Double2D Ωⁿ_γ;
		private Double2D A0_γ;
		private Double2D ΩⁿA0_γ;
		private Double2D Ψ_A0ʹΩⁿA0_γ;
		private Integer1D z;
		private Integer1D R;
		private int max_r;
		private int[] γ_indices;
		
		private InverseWishartDist iwDist = new InverseWishartDist();

		public PostSigmaDist() { }

		@Override
		public void setMCState(ChainLink l) {
			this.γ_indices = ((Integer1D) l.get("gamma").getNumericValue()).which(1).value();
			this.X_γ = X.subMatrix(null, γ_indices);
			this.Ωⁿ_γ = ((Double2D) l.get("Omega").getNumericValue()).subMatrix(γ_indices, γ_indices).inverse();
			this.A0_γ = ((Double2D) l.get("A0").getNumericValue()).subMatrix(γ_indices, null);
			this.ΩⁿA0_γ = Ωⁿ_γ.mult(A0_γ);
			this.Ψ_A0ʹΩⁿA0_γ = Ψ.plus(A0_γ.transposeMult(ΩⁿA0_γ)); 
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
				
				Double2D Xᵣ_γ = X_γ.getAll(zᵣ);
				Double2D Yᵣ = Y.getAll(zᵣ);
				
				Double2D XᵣʹXᵣ_γ = Xᵣ_γ.transposeMult(Xᵣ_γ); // FIXME: here - extended chain
				Double2D XᵣʹYᵣ_γ = Xᵣ_γ.transposeMult(Yᵣ); // FIXME: here - extended chain
				Double2D YᵣʹYᵣ = Yᵣ.transposeMult(Yᵣ); // FIXME: here - extended chain
				
				Double2D Ωⁿ_XᵣʹXᵣ_γ = Ωⁿ_γ.plus(XᵣʹXᵣ_γ);
				Double2D ΩⁿA0_XᵣʹYᵣ_γ = ΩⁿA0_γ.plus(XᵣʹYᵣ_γ);
				Double2D Ψᵣᵤ_γ = Ψ_A0ʹΩⁿA0_γ.plus(YᵣʹYᵣ).minus(ΩⁿA0_XᵣʹYᵣ_γ.transposeMult(Ωⁿ_XᵣʹXᵣ_γ.inverse().mult(ΩⁿA0_XᵣʹYᵣ_γ)));
				
				iwDist.setParms(Ψᵣᵤ_γ, κ.plus(Nᵣ));
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
		private int p = (Integer) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("p");
		private int h = (Integer) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("h");
		private Double2D Y = ((Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("Y"));
		private Double2D X = ((Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("X"));

		// Chain parameters
		private Double2D X_γ;
		private Integer1D z;
		private Double2D Ωⁿ_γ;
		private Double3D Σ;
		private Double2D A0_γ;
		private Double2D ΩⁿA0_γ;
		private Integer1D R;
		private int max_r;
		private int[] γ_indices;
		
		private MatrixNormalDist mnDist_γ;

		public PostADist() { }

		@Override
		public void setMCState(ChainLink l) {
			this.γ_indices = ((Integer1D) l.get("gamma").getNumericValue()).which(1).value();
			this.X_γ = X.subMatrix(null,γ_indices);
			this.A0_γ = ((Double2D) l.get("A0").getNumericValue()).subMatrix(γ_indices, null);
			this.Ωⁿ_γ = ((Double2D) l.get("Omega").getNumericValue()).subMatrix(γ_indices,γ_indices).inverse();
			this.ΩⁿA0_γ = Ωⁿ_γ.mult(A0_γ);
			this.mnDist_γ = new MatrixNormalDist();			
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
				Double2D Xᵣ_γ = X_γ.getAll(zᵣ);
				Double2D Yᵣ = Y.getAll(zᵣ);
				
				Double2D XᵣʹXᵣ_γ = Xᵣ_γ.transposeMult(Xᵣ_γ); // FIXME: here - extended chain
				Double2D Ωᵣᵤ_γ = Ωⁿ_γ.plus(XᵣʹXᵣ_γ).inverse(); 
				
				Double2D XᵣʹYᵣ_γ = Xᵣ_γ.transposeMult(Yᵣ); // FIXME: here - extended chain
				Double2D Α0ᵣᵤ_γ = Ωᵣᵤ_γ.mult(XᵣʹYᵣ_γ.plus(ΩⁿA0_γ));
				
				mnDist_γ.setParms(Α0ᵣᵤ_γ, Σ.get(r), Ωᵣᵤ_γ);
				Double2D Aᵣ_γ = mnDist_γ.variate();
				
				int slice = r; // FIXME: this is what you would change for "compression"
				int row_Aᵣ_γ = 0;
				for (int row : γ_indices) {
					for (int col=0; col<p; col++) {
						A.setQuick(slice, row, col, Aᵣ_γ.getDm().getQuick(row_Aᵣ_γ, col));
					}
					row_Aᵣ_γ++;
				}
			}
								
			return new Double3D(A);
		}
	}

	@SuppressWarnings("unchecked")
	class CRPDist extends ProbDistMC<Integer1D> {
		// Fixed parameters
		private double λ = ((Double0D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("λ")).value();
		private Double0D κ_1 = ((Double0D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("κ_1"));
		private Double2D Ψ = ((Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("Ψ"));
		private Double2D Φ = ((Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("Φ"));
		private Double2D W = ((Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("W"));
		Double2D Sⁿ = (Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("Sⁿ");
		
		int z_lag = (Integer) MLM_sample_params_varbsel_block_cached.this.extended_hypers.get("z_lag");
		double ldet_Sⁿ = (Double) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("ldet_Sⁿ");
		double lτ = (Double) (MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("lτ"));
		double lτʹ = (Double) (MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("lτʹ"));
		double ldet_Ψ = (Double) (MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("ldet_Ψ"));
		private double lπ = (Double) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("lπ");
		private double l2π = (Double) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("l2π");
		private int N = (Integer) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("N");
		private int p = (Integer) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("p");
		private int h = (Integer) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("h");
		private Double2D X = ((Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("X"));
		private Double2D Y = ((Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("Y"));
		private Double2D[] yyʹ = (Double2D[]) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("yyʹ");
		private double lΓκ = (Double) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("lΓκ");
		private double lΓκ_1 = (Double) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("lΓκ_1");
		Map<Integer1D,Double2D[]> xxʹ_γs = (Map<Integer1D,Double2D[]>) extended_hypers.get("xxʹ_γs"); // FIXME
		Map<Integer1D,Double2D[]> xyʹ_γs = (Map<Integer1D,Double2D[]>) extended_hypers.get("xyʹ_γs"); // FIXME
		
		// Chain parameters
		private Double2D A0;
		private Double2D A0_W;
		private Double2D Ω;
		private Double3D Σ; // FIXME: Blocked, modified
		private Double3D A; // FIXME: Blocked, modified
		private Integer1D γ; // FIXME: Blocked, modified
		private Integer1D z;
		private Double0D α;
		
		private double Tz;
		private double Tγ;

		private Double1D postProb;
		private CategoricalDist catDist = new CategoricalDist();
		private InverseWishartDist iwDist = new InverseWishartDist();

		public CRPDist() { }

		@Override
		public void setMCState(ChainLink l) {
			this.γ = (Integer1D) l.get("gamma").getNumericValue();
			this.Ω = ((Double2D) l.get("Omega").getNumericValue());
			this.A0 = (Double2D) l.get("A0").getNumericValue();
			this.A0_W = A0.minus(W);
			this.Σ = (Double3D) l.get("Sigma").getNumericValue();
			this.A = (Double3D) l.get("A").getNumericValue();
			this.z = (Integer1D) l.get("z").getNumericValue();
			this.α = (Double0D) l.get("alpha").getNumericValue();
			
			this.Tz = ((Double0D) l.get("Tz").getNumericValue()).value();
			this.Tγ = ((Double0D) l.get("Tγ").getNumericValue()).value();
			
			for (int i=0; i<N; i++) {
				Double1D yᵢ = Y.get(i);
				yyʹ[i] = yᵢ.outer(yᵢ);
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
			for (int varb = 1; varb<h; varb++) {
				Integer1D[] z_γ = new Integer1D[2];
				double[] logP_γz = new double[2];
				logP_γz[0] = lτʹ;
				logP_γz[1] = lτ;
				int γ_val_old = γ.value()[varb];
				for (int γ_val = 0; γ_val<2; γ_val++) {
					γ.set(varb,γ_val);
					int[] γ_indices = γ.which(1).items().value();
					
					Double2D Ω_γ = Ω.subMatrix(γ_indices,γ_indices);
					Double2D Ωⁿ_γ = Ω_γ.inverse();
					double ldet_Ωⁿ_γ = 0.5*p*Math.log(Ωⁿ_γ.det());	
					Double2D A0_γ = A0.subMatrix(γ_indices, null);
					Double2D ΩⁿA0_γ = Ωⁿ_γ.mult(A0_γ);
					Double2D Ψ_A0ʹΩⁿA0_γ = Ψ.plus(A0_γ.transposeMult(ΩⁿA0_γ)); 
					Double2D Φ_γ = Φ.subMatrix(γ_indices,γ_indices);
					Double2D A0_W_γ = A0_W.subMatrix(γ_indices,null);
					Double2D[] xxʹ_γ = xxʹ_γs.get(γ);
					Double2D[] xyʹ_γ = xyʹ_γs.get(γ);
					MatrixNormalDist mnDist = new MatrixNormalDist();
					
					int h_γ = γ_indices.length;
					double λ_γ = λ - (h - h_γ);
					logP_γz[γ_val] = 0.5*λ_γ*Math.log(Φ_γ.det()) - 0.5*(λ_γ+h_γ+1)*Math.log(Ω_γ.det()) - 0.5*h_γ*λ_γ*Math.log(2) - lΓ(λ_γ/2.0,h_γ) - 0.5*Φ_γ.mult(Ω_γ.inverse()).trace().value()
							    - h_γ*l2π + 0.5*h_γ*ldet_Sⁿ - 0.5*Sⁿ.mult(A0_W_γ.transposeMult(A0_W_γ)).trace().value(); // Forget log(1)!
					
					z_γ[γ_val] = new Integer1D(z.value().clone());
					for (int l = 0; l < z_lag; l++) {
					for (int i = 0; i < N; i++) {
						Integer1D R = z_γ[γ_val].items();
						int zi_old_value = z_γ[γ_val].value()[i];
						z_γ[γ_val].set(i,-1);
						
						Double1D xᵢ = X.get(i);
						int Nᵣᵢ = z_γ[γ_val].which(zi_old_value).size();
					
						double[] logP_z = new double[R.size() + 1];	
						double maxLogP = -Double.MAX_VALUE;										
						
						// Existing categories
						int rj = 0;
							for (rj=0; rj<R.size(); rj++) {
								int r = R.value()[rj];
								
								Double1D Aᵣʹxᵢ = A.get(r).sparsifyRows(γ_indices).transposeMult(xᵢ);
								Double1D yᵢ_Aᵣʹxᵢ = Y.get(i).minus(Aᵣʹxᵢ);
								
								double llik_e_term = -0.5*yᵢ_Aᵣʹxᵢ.mult(Σ.get(r).inverse().mult(yᵢ_Aᵣʹxᵢ));
								double llik_nc_term = -l2π - 0.5*Math.log(Σ.get(r).det());
								if (r == zi_old_value) {
									if (Nᵣᵢ == 0) {
										logP_z[rj] = -Double.MAX_VALUE;
									} else {
										logP_z[rj] = Math.log(Nᵣᵢ) + llik_nc_term + llik_e_term;
									}
								} else {
									int Nᵣ = z_γ[γ_val].which(r).size();
									logP_z[rj] = Math.log(Nᵣ) + llik_nc_term + llik_e_term;
								}
								if (logP_z[rj] == Double.POSITIVE_INFINITY) {
									logP_z[rj] = Double.MAX_VALUE;
								}				
								if (logP_z[rj] > maxLogP) {
									maxLogP = logP_z[rj];
								}
							}	
							
							// New category
							Double2D Ωᵤ_γ = Ωⁿ_γ.plus(xxʹ_γ[i]).inverse();
							Double2D ΩⁿA0_xᵢyᵢʹ_γ = ΩⁿA0_γ.plus(xyʹ_γ[i]);
							Double2D Ωᵤ_ΩⁿA0_xᵢyᵢʹ_γ = Ωᵤ_γ.mult(ΩⁿA0_xᵢyᵢʹ_γ);
							Double2D Ψᵤ = Ψ_A0ʹΩⁿA0_γ.plus(yyʹ[i]).minus(ΩⁿA0_xᵢyᵢʹ_γ.transposeMult(Ωᵤ_ΩⁿA0_xᵢyᵢʹ_γ));
							double ldet_Ωᵤ_γ = 0.5*p*Math.log(Ωᵤ_γ.det());
							double ldet_Ψᵤ = 0.5*κ_1.value()*Math.log(Ψᵤ.det());
							double prior_lil = lΓκ_1 + ldet_Ψ + ldet_Ωⁿ_γ + ldet_Ωᵤ_γ - lπ - lΓκ - ldet_Ψᵤ;
							logP_z[rj] = lα + prior_lil;
							if (logP_z[rj] == Double.POSITIVE_INFINITY) {
								logP_z[rj] = Double.MAX_VALUE;
							}				
							if (logP_z[rj] > maxLogP) {
								maxLogP = logP_z[rj];
							}			
							
							// Normalize and do annealing
							Double1D lP_scaled = (new Double1D(logP_z)).minus(maxLogP);
							Double1D prob_unnorm = lP_scaled.exp();
							double p_z_nc = prob_unnorm.sum().value();
							double lp_z_nc = Math.log(p_z_nc);
							double[] logP_z_annealed = new double[logP_z.length];
							int newc_rjindex = rj; // last index
							maxLogP = -Double.MAX_VALUE;										
							for (rj=0; rj<logP_z.length; rj++) {
								if (rj == newc_rjindex || (R.value()[rj] != zi_old_value)) {
									logP_z_annealed[rj] = (1-1/Tz)*lp_z_nc + (1/Tz)*lP_scaled.value()[rj];
								} else {
									logP_z_annealed[rj] = lP_scaled.value()[rj];
								}
								if (logP_z_annealed[rj] > maxLogP) {
									maxLogP = logP_z_annealed[rj];
								}
							}
							
							// Select a category for this point
							catDist.setParms((new Double1D(logP_z_annealed)).minus(maxLogP).exp());
							int zi_new_rjindex = catDist.variate().value();
							int zi_new_value;
							if (zi_new_rjindex == newc_rjindex) {
								int maxActive = R.value()[R.size()-1];
								zi_new_value = z_γ[γ_val].minNotIn(0, maxActive);
							} else {
								zi_new_value = R.value()[zi_new_rjindex];
							}						
							
							// Sample new category if necessary
							if (zi_new_rjindex == newc_rjindex) {
								// Sample Σ
								iwDist.setParms(Ψᵤ, κ_1);
								Double2D Σᵢ = iwDist.variate();
								Σ.set(zi_new_value, Σᵢ);
								// Sample A
								mnDist.setParms(Ωᵤ_ΩⁿA0_xᵢyᵢʹ_γ, Σᵢ, Ωᵤ_γ);
								Double2D Aᵢ_γ = mnDist.variate();
								A.set(zi_new_value, new Double2D(h,p));
								int row_Aᵢ_γ = 0;
								for (int row : γ_indices) {
									for (int col=0; col<p; col++) {
										A.getDm().set(zi_new_value, row, col, Aᵢ_γ.getDm().get(row_Aᵢ_γ, col));
									}
									row_Aᵢ_γ++;
								}
							}
							z_γ[γ_val].set(i, zi_new_value);
							logP_γz[γ_val] += logP_z[zi_new_rjindex];
						}
					}
				}
				// Sample γ
				double maxLogP = logP_γz[0] > logP_γz[1] ? logP_γz[0] : logP_γz[1];
				
				// Normalize and do annealing
				Double1D prob_unnorm = (new Double1D(logP_γz)).minus(maxLogP).exp();
				double p_γz_nc = prob_unnorm.sum().value();
				double lp_γz_nc = Math.log(p_γz_nc);
				double[] logP_γz_annealed = new double[2];
				for (int γ_val=0; γ_val<2; γ_val++) {
					if (γ_val != γ_val_old) {
						logP_γz_annealed[γ_val] = (1-1/Tγ)*lp_γz_nc + (1/Tγ)*logP_γz[γ_val];
					} else {
						logP_γz_annealed[γ_val] = logP_γz[γ_val];
					}
				}
					
				catDist.setParms((new Double1D(logP_γz_annealed)).minus(maxLogP).exp());
				γ.set(varb,catDist.variate().value());
				z = z_γ[γ.value()[varb]];
			}
			return z;
		}
	}

	class PostA0Dist extends ProbDistMC<Double2D> {
		// Fixed parameters
		private Double2D Sⁿ = (Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("Sⁿ");
		private Double2D W = (Double2D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("W");
		private int h = (Integer) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("h");
		private int p = (Integer) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("p");
		private MVNormalDist mvnDist;

		// Chain parameters
		private Integer1D z;
		private Double2D Ωⁿ_γ;
		private Double1D vecW_γ;
		private Double3D Σ;
		private Double3D A;
		private int[] γ_indices;
		private int h_γ;

		public PostA0Dist() { }

		@Override
		public void setMCState(ChainLink l) {
			this.γ_indices = ((Integer1D) l.get("gamma").getNumericValue()).which(1).value();
			this.h_γ = γ_indices.length;
			this.z = (Integer1D) l.get("z").getNumericValue();
			this.Ωⁿ_γ = ((Double2D) l.get("Omega").getNumericValue()).subMatrix(γ_indices,γ_indices).inverse();
			this.vecW_γ = W.subMatrix(γ_indices,null).colVec();
			this.Σ = (Double3D) l.get("Sigma").getNumericValue();
			this.A = (Double3D) l.get("A").getNumericValue();
			this.mvnDist = new MVNormalDist();
			this.initialized = true;
		}

		@Override
		protected Double2D genVariate() {
			int[] R = z.items().value();
			Double2D precision = Sⁿ.kron(Double2D.ident(h_γ));
			Double1D unscaled_centre = precision.mult(vecW_γ);
			for (int r : R) {
				Double2D Σᵣⁿ = Σ.get(r).inverse();
				Double2D Σᵣⁿ_kr_Ωⁿ_γ = Σᵣⁿ.kron(Ωⁿ_γ);
				Double1D vecAᵣ_γ = A.get(r).subMatrix(γ_indices, null).colVec();
				precision.plusEquals(Σᵣⁿ_kr_Ωⁿ_γ);
				unscaled_centre.plusEquals(Σᵣⁿ_kr_Ωⁿ_γ.mult(vecAᵣ_γ));
			}
			Double2D covariance = precision.inverse();
			mvnDist.setParms(covariance.mult(unscaled_centre), covariance);		
			Double2D A0_γ = this.mvnDist.variate().toDouble2D(h_γ);
			
			DoubleMatrix2D A0_m = new DenseDoubleMatrix2D(h,p);
			int row_A0_γ = 0;
			for (int row : γ_indices) {
				for (int col=0; col<p; col++) {
					A0_m.set(row, col, A0_γ.getDm().get(row_A0_γ, col));
				}
				row_A0_γ++;
			}
			return new Double2D(A0_m);
		}
	}

	class PostAlphaDist extends ProbDistMC<Double0D> {
		private BetaDist betaDist = new BetaDist();
		private GammaDist gammaDist = new GammaDist();
		private CategoricalDist catDist = new CategoricalDist();

		// Fixed parameters
		private Double0D alpha_a = (Double0D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("alpha_a");
		private Double0D alpha_b = (Double0D) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("alpha_b");
		private Double0D N = new Double0D((Integer) MLM_sample_params_varbsel_block_cache.this.extended_hypers.get("N")); 

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
		Double2D Φ = (Double2D) extended_hypers.get("Φ");
		Double0D λ = (Double0D) extended_hypers.get("λ");	
		Double2D Ψ = (Double2D) extended_hypers.get("Ψ");
		double ldet_Ψ = (Double) extended_hypers.get("ldet_Ψ");
		double κ = ((Double0D) extended_hypers.get("κ")).value();	
		double κ_1 = ((Double0D) extended_hypers.get("κ_1")).value();	
		double lΓκ = (Double) extended_hypers.get("lΓκ");
		double lΓκ_1 = (Double) extended_hypers.get("lΓκ_1");
		int h = (Integer) extended_hypers.get("h");
		int p = (Integer) extended_hypers.get("p");
		int N = (Integer) extended_hypers.get("N");
		double l2π = (Double) extended_hypers.get("l2π");
		double lπ = (Double) extended_hypers.get("lπ");
		Double2D W = (Double2D) extended_hypers.get("W");
		Double2D Sⁿ = (Double2D) extended_hypers.get("Sⁿ");
		double lτ = (Double) extended_hypers.get("lτ");
		double lτʹ = (Double) extended_hypers.get("lτʹ");
		double ldet_Sⁿ = (Double) extended_hypers.get("ldet_Sⁿ");
		
		Double2D Y = (Double2D) extended_hypers.get("Y");
		Double2D[] yyʹ = (Double2D[]) extended_hypers.get("yyʹ");
		Map<Integer1D,Double2D> X_γs = (Map<Integer1D,Double2D>) extended_hypers.get("X_γs");
		Map<Integer1D,Double2D[]> xxʹ_γs = (Map<Integer1D,Double2D[]>) extended_hypers.get("xxʹ_γs"); // FIXME
		Map<Integer1D,Double2D[]> xyʹ_γs = (Map<Integer1D,Double2D[]>) extended_hypers.get("xyʹ_γs"); // FIXME
		
		for (int m = 0; m < c.size(); m++) {
			double posterior = 0;
			ChainLink cl_m = c.get(m);
			
			// γ
			Integer1D γ = ((RandomVar<Integer1D>) cl_m.get("gamma")).getNumericValue();
			for (int i=1; i<h; i++) {
				posterior += γ.value()[i] == 1 ? lτ : lτʹ;
			}
			int[] γ_indices = γ.which(1).value();
			
			// α
			Double0D α = ((RandomVar<Double0D>) cl_m.get("alpha")).getNumericValue();
			posterior += gammaDist.logDensity(α);
						
			// Omega
			Double2D Ω = ((RandomVar<Double2D>) cl_m.get("Omega")).getNumericValue();
			Double2D Ω_γ = Ω.subMatrix(γ_indices, γ_indices);
			Double2D Ωⁿ_γ = Ω_γ.inverse();
			double ldet_Ω_γ = Math.log(Ω_γ.det());
			
			int h_γ = γ_indices.length;
			
			double λ_γ = λ.value() - (h - h_γ);
			Double2D Φ_γ = Φ.subMatrix(γ_indices, γ_indices);
			double ldet_Φ_γ = Math.log(Φ_γ.det());
			
			posterior += 0.5*λ_γ*ldet_Φ_γ - 0.5*λ_γ*h_γ*Math.log(2) - lΓ(λ_γ/2., h_γ) - 0.5*(λ_γ+h_γ+1)*ldet_Ω_γ - 0.5*Φ_γ.mult(Ωⁿ_γ).trace().value();
			
			// A0
			Double2D A0 = ((RandomVar<Double2D>) cl_m.get("A0")).getNumericValue();
			Double2D A0_γ = A0.subMatrix(γ_indices, null);
			Double2D A0_W_γ = A0_γ.minus(W.subMatrix(γ_indices, null));
			posterior += -h_γ*l2π + h_γ*ldet_Sⁿ - 0.5*Sⁿ.mult(A0_W_γ.transposeMult(A0_W_γ)).trace().value();
	
			// Model parameters and likelihood
			Double3D Σ = ((RandomVar<Double3D>) cl_m.get("Sigma")).getNumericValue();
			Double3D A = ((RandomVar<Double3D>) cl_m.get("A")).getNumericValue();
			double lα = Math.log(α.value());
			int[] z = ((RandomVar<Integer1D>) cl_m.get("z")).getNumericValue().value();
//			Map<Integer,Object> cat_changed = new HashMap<Integer,Object>();
			
			double ldet_Ωⁿ_γ = Math.log(Ωⁿ_γ.det());
			Double2D ΩⁿA0_γ = Ωⁿ_γ.mult(A0_γ);
			Double2D Ψ_A0ʹΩⁿA0_γ = Ψ.mult(A0_γ.transposeMult(ΩⁿA0_γ));
			
			Double2D X_γ = X_γs.get(γ); 
			Double2D[] xxʹ_γ = xxʹ_γs.get(γ);
			Double2D[] xyʹ_γ = xyʹ_γs.get(γ);
			for (int i = 0; i < N; i++) {
				int rᵢ = z[i];
				Double1D xᵢ_γ = X_γ.get(i);
				Double2D Σᵢ = Σ.get(rᵢ);
				Double2D Aᵢ = A.get(rᵢ);
				Double2D Aᵢ_γ = Aᵢ.subMatrix(γ_indices, null);
				Double2D Aᵢ_A0_γ = Aᵢ_γ.minus(A0_γ);
//				
				if (i > 0) {
					int[] zᵢ_ints = new int[i];
					System.arraycopy(z,0,zᵢ_ints,0,i); // FIXME - yeah, this is dumb - fix it later
					Integer1D zᵢ = new Integer1D(zᵢ_ints);
					Integer1D Rᵢ = zᵢ.items();
//		
					double ll_curr = 0;
					double nc = 0;
					for (int j_r=0; j_r<Rᵢ.size(); j_r++) {
						int r = Rᵢ.value()[j_r];
						int[] zᵣ = zᵢ.which(r).value();
						int Nᵣ = zᵣ.length; // FIXME: sub-cache
						
						Double1D Aᵣʹxᵢ_γ = Aᵢ_γ.transposeMult(xᵢ_γ);
						Double1D yᵢ_Aᵣʹxᵢ_γ = Y.get(i).minus(Aᵣʹxᵢ_γ);
							
						double llik_e_term = -0.5*yᵢ_Aᵣʹxᵢ_γ.mult(Σ.get(r).inverse().mult(yᵢ_Aᵣʹxᵢ_γ));
						double llik_nc_term = -l2π - 0.5*Math.log(Σ.get(r).det());
						double ll = Math.log(Nᵣ) + llik_nc_term + llik_e_term;
						
						if (r == rᵢ) {
							ll_curr = ll;
						}
						nc += Math.exp(ll);
					}
					// Compute the 'new category' term
					Double2D Ωᵤ_γ = Ωⁿ_γ.plus(xxʹ_γ[i]).inverse();
					Double2D ΩⁿA0_xᵢyᵢʹ_γ = ΩⁿA0_γ.plus(xyʹ_γ[i]);
					Double2D Ωᵤ_ΩⁿA0_xᵢyᵢʹ_γ = Ωᵤ_γ.mult(ΩⁿA0_xᵢyᵢʹ_γ);
					Double2D Ψᵤ_γ = Ψ_A0ʹΩⁿA0_γ.plus(yyʹ[i]).minus(ΩⁿA0_xᵢyᵢʹ_γ.transposeMult(Ωᵤ_ΩⁿA0_xᵢyᵢʹ_γ));
					double ldet_Ωᵤ_γ = 0.5*p*Math.log(Ωᵤ_γ.det());
					double ldet_Ψᵤ_γ = 0.5*κ_1*Math.log(Ψᵤ_γ.det());	
					double prior_lil = lΓκ_1 + ldet_Ψ + ldet_Ωⁿ_γ + ldet_Ωᵤ_γ - lπ - lΓκ - ldet_Ψᵤ_γ;
					
					double ll = lα + prior_lil;
					nc += Math.exp(ll);	
					if (zᵢ.which(rᵢ).size() == 0) {
						ll_curr = ll;
					}
					
					posterior += ll_curr - Math.log(nc);
				} else { // i == 0: probability 1
					posterior += 0;
				}
				posterior += ldet_Ψ - 0.5*p*κ*Math.log(2) - lΓκ - 0.5*(κ+p+1)*Math.log(Σᵢ.det()) - 0.5*Ψ.mult(Σᵢ.inverse()).trace().value();
				posterior += -h_γ*l2π - 0.5*h_γ*Math.log(Σᵢ.det()) - 0.5*(Σᵢ.inverse().mult(Aᵢ_A0_γ.transposeMult(Ωⁿ_γ.mult(Aᵢ_A0_γ)))).trace().value();
				
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
			}				
			
			// Set max
			if (posterior > max) {
				argmax = m;
				max = posterior;
			}
		}

		return c.get(argmax);			
	}

	public MLM_sample_params_varbsel_block_cache() {
		super();
	}
	
	public MLM_sample_params_varbsel_block_cache(Map<String, Flattenable> hypers, Map<String, Flattenable> init, Double2D data) {
		super(hypers);
		
		// Set up parameters
		this.params = new HashMap<String, RandomVar<? extends Flattenable>>();
		
		// Extend hyperparameters
		this.extended_hypers = MLM_sample_params_varbsel_block_cache.extend_hypers(hypers, data, (Double2D) init.get("X"));	
		
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

		// gamma
		this.params.put("gamma", new RandomVar<Integer1D>("gamma", null, (Integer1D) init.get("gamma")));

		// alpha
		this.params.put("alpha", new RandomVar<Double0D>("alpha", new PostAlphaDist(), (Double0D) init.get("alpha")));
		
		// Temperature (z)
		int deadline = ((Integer0D) hypers.get("deadline")).value();
		
		double T0z = ((Double0D) hypers.get("T0z")).value();
		double Tfz = ((Double0D) hypers.get("Tfz")).value();
		this.params.put("Tz", new RandomVar<Double0D>("Tz", new TemperatureSchedule(T0z, Tfz, deadline), (Double0D) hypers.get("T0z")));
		
		// Temperature (γ)
		double T0g = ((Double0D) hypers.get("T0g")).value();
		double Tfg = ((Double0D) hypers.get("Tfg")).value();
		this.params.put("Tγ", new RandomVar<Double0D>("Tγ", new TemperatureSchedule(T0g, Tfg, deadline), (Double0D) hypers.get("T0g")));	
	}

	@Override
	public ChainLink getInitialLink() {
		ChainLink cl = new ChainLink();
		cl.add(this.params.get("z"));
		cl.add(this.params.get("gamma"));
		cl.add(this.params.get("A"));
		cl.add(this.params.get("Sigma"));
		cl.add(this.params.get("Omega"));
		cl.add(this.params.get("A0"));
		cl.add(this.params.get("alpha"));
		cl.add(this.params.get("Tz"));
		cl.add(this.params.get("Tγ"));
		return cl;
	}
}
