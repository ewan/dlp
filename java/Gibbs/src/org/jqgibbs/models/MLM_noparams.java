package org.jqgibbs.models;

import java.lang.reflect.Array;
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
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Integer1D;
import org.jqgibbs.mathstat.probdist.BetaDist;
import org.jqgibbs.mathstat.probdist.CategoricalDist;
import org.jqgibbs.mathstat.probdist.GammaDist;
import org.jqgibbs.mathstat.probdist.ProbDistMC;

import cern.jet.stat.Gamma;

public class MLM_noparams extends Model {
	
//	private static double Γ(double x, int d) { // FIXME
//		if (d == 1) {
//			return Gamma.gamma(x);
//		} else { // FIXME check d0
//			return Math.pow(Math.PI, (d - 1) / 2) * MLM_noparams.Γ(x, d - 1)
//					* Gamma.gamma(x + (1 - d) / 2);
//		}
//	}
	
	private static double lΓ(double x, int p) {
		double r = 0.25*p*(p-1)*Math.log(Math.PI);
		for (int i=0; i<p; i++) {
			r += Gamma.logGamma(x+i/2);
		}
		return r;
	}
	
	private static Map<String,Object> extend_hypers(Map<String,Flattenable> hypers, Double2D Y, Double2D X) {
		Map<String,Object> extended_hypers = new HashMap<String,Object>();
		
		double κ = ((Double0D) hypers.get("kappa")).value();
		double κ_1 = κ+1;
		extended_hypers.put("κ", κ);
		
		int p = Y.numCols();
		extended_hypers.put("p", p);
		
		double lπ = 0.5*p*Math.log(Math.PI);
		extended_hypers.put("lπ", lπ);
			
		Double2D Ψ = ((Double2D) hypers.get("Psi"));
		Double2D Ωⁿ = ((Double2D) hypers.get("Omega")).inverse();
		Double2D A0 = (Double2D) hypers.get("A0");
		Double2D ΩⁿA0 = Ωⁿ.mult(A0);
		Double2D Ψ_A0ʹΩⁿA0 = Ψ.plus(A0.transposeMult(ΩⁿA0));
		double lΓκ = lΓ(κ, p);
		double lΓκ_1 = lΓ(κ_1, p);
		double ldet_Ψ = 0.5*κ*Math.log(Ψ.det());
		double ldet_Ωⁿ = 0.5*p*Math.log(Ωⁿ.det());	
			
		extended_hypers.put("Ωⁿ", Ωⁿ);
		extended_hypers.put("A0", A0);
		extended_hypers.put("ΩⁿA0", ΩⁿA0);
		extended_hypers.put("Ψ_A0ʹΩⁿA0", Ψ_A0ʹΩⁿA0);
	
		extended_hypers.put("Y", Y);
		extended_hypers.put("X", X);
	
		int N = Y.numRows();
		extended_hypers.put("N", N);
		
		Double2D[] xxʹ = new Double2D[N];
		Double2D[] yyʹ = new Double2D[N];
		Double2D[] xyʹ = new Double2D[N];
		double[] prior_lil = new double[N];
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
			double ldet_Ψᵤ = 0.5*κ_1*Math.log(Ψᵤ.det());				
			
			prior_lil[i] = lΓκ_1 + ldet_Ψ + ldet_Ωⁿ + ldet_Ωᵤ - lπ - lΓκ - ldet_Ψᵤ;
		}
		extended_hypers.put("xxʹ", xxʹ);
		extended_hypers.put("xyʹ", xyʹ);
		extended_hypers.put("yyʹ", yyʹ);
		extended_hypers.put("prior_lil", prior_lil);
		
		extended_hypers.put("alpha_a", (Double0D) hypers.get("alpha_a"));
		extended_hypers.put("alpha_b", (Double0D) hypers.get("alpha_b"));
		
		return extended_hypers;
	}
	
	class PostAlphaDist extends ProbDistMC<Double0D> {
		private BetaDist betaDist = new BetaDist();
		private GammaDist gammaDist = new GammaDist();
		private CategoricalDist catDist = new CategoricalDist();

		// Fixed parameters
		private Double0D alpha_a = (Double0D) MLM_noparams.this.extended_hypers.get("alpha_a");
		private Double0D alpha_b = (Double0D) MLM_noparams.this.extended_hypers.get("alpha_b");
		private Double0D N = new Double0D((Integer) MLM_noparams.this.extended_hypers.get("N")); 

		// Chain parameters
		private Double0D α;
		private Integer1D z;

		public PostAlphaDist() { }

		@Override
		public void setMCState(ChainLink l) {
			this.α = (Double0D) l.get("alpha").getNumericValue();
			this.z = (Integer1D) l.get("z").getNumericValue();
			this.setUpFromParms();
			this.initialized = true;
		}

		protected void setUpFromParms() {
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
		}

		@Override
		protected Double0D genVariate() {
			return gammaDist.variate();
		}
	}
	
	class CRPDist extends ProbDistMC<Integer1D> {
		// Fixed parameters
		private double κ = ((Double) MLM_noparams.this.extended_hypers.get("κ"));
		private Double2D Ωⁿ = (Double2D) MLM_noparams.this.extended_hypers.get("Ωⁿ");
		private Double2D Y = (Double2D) MLM_noparams.this.extended_hypers.get("Y");
		private Double2D X = (Double2D) MLM_noparams.this.extended_hypers.get("X");
		
		private int p = (Integer) MLM_noparams.this.extended_hypers.get("p");
		private int N = (Integer) MLM_noparams.this.extended_hypers.get("N");
		private double lπ = (Double) MLM_noparams.this.extended_hypers.get("lπ");
		private Double2D ΩⁿA0 = (Double2D) MLM_noparams.this.extended_hypers.get("ΩⁿA0");
		private Double2D Ψ_A0ʹΩⁿA0 = (Double2D) MLM_noparams.this.extended_hypers.get("Ψ_A0ʹΩⁿA0");
		private Double2D[] xxʹ = (Double2D[]) MLM_noparams.this.extended_hypers.get("xxʹ");
		private Double2D[] xyʹ = (Double2D[]) MLM_noparams.this.extended_hypers.get("xyʹ");
		private Double2D[] yyʹ = (Double2D[]) MLM_noparams.this.extended_hypers.get("yyʹ");
		private double[] prior_lil = (double[]) MLM_noparams.this.extended_hypers.get("prior_lil");
		
		private Map<Integer,Double2D> XᵣʹXᵣ = new HashMap<Integer,Double2D>();
		private Map<Integer,Double2D> XᵣʹYᵣ = new HashMap<Integer,Double2D>();
		private Map<Integer,Double2D> YᵣʹYᵣ = new HashMap<Integer,Double2D>();
		private Map<Integer,Double2D> Ψ_A0ʹΩⁿA0_YᵣʹYᵣ = new HashMap<Integer,Double2D>();
		private Map<Integer,Double2D> ΩⁿA0_XᵣʹYᵣ = new HashMap<Integer,Double2D>();
		private Map<Integer,Double2D> Ωⁿ_XᵣʹXᵣ = new HashMap<Integer,Double2D>();
		private Map<Integer,Double2D> Ωᵥ = new HashMap<Integer,Double2D>();
		private Map<Integer,Double2D> Ψᵥ = new HashMap<Integer,Double2D>();
		private Map<Integer,Double> ldet_Ψᵥ = new HashMap<Integer,Double>();
		private Map<Integer,Double> ldet_Ωᵥ = new HashMap<Integer,Double>();
		private Map<Integer,Integer> Nᵣ = new HashMap<Integer,Integer>();

		// Chain parameters
		private Integer1D z;
		private Double0D α;

		private Double1D postProb;
		private CategoricalDist catDist = new CategoricalDist();

		public CRPDist() { }

		@Override
		public void setMCState(ChainLink l) {
			this.z = (Integer1D) l.get("z").getNumericValue(); // FIXME - below
			this.α = (Double0D) l.get("alpha").getNumericValue();
			this.setUpFromParms();
			this.initialized = true;
		}

		protected Integer0D catVariate() {
			this.catDist.setParms(this.postProb);
			return this.catDist.variate();
		}

		protected void setUpFromParms() { }

		@Override
		protected Integer1D genVariate() {
			double lα = Math.log(α.value());
			
			Integer1D z = new Integer1D(this.z.value().clone());
			
			for (int i = 0; i < N; i++) {
				Integer1D R = z.items();
				int zi_old_value = z.value()[i];
				z.set(i, -1);
						
				int Nᵣᵢ;
				Double2D XᵣᵢʹXᵣᵢ = null;
				Double2D YᵣᵢʹYᵣᵢ = null;
				Double2D XᵣᵢʹYᵣᵢ = null;
				Double2D Ψ_A0ʹΩⁿA0_YᵣᵢʹYᵣᵢ = null;
				Double2D ΩⁿA0_XᵣᵢʹYᵣᵢ = null;
				Double2D Ωⁿ_XᵣᵢʹXᵣᵢ = null;
				Double2D Ωᵥᵢ = null;
				Double2D Ψᵥᵢ = null;
				double ldet_Ψᵥᵢ = 0;
				double ldet_Ωᵥᵢ = 0;
				if (XᵣʹXᵣ.containsKey(zi_old_value)) {
					// Just remove the current point
					Nᵣᵢ = Nᵣ.get(zi_old_value)-1;
					if (Nᵣᵢ > 0) {
						XᵣᵢʹXᵣᵢ = XᵣʹXᵣ.get(zi_old_value).minus(xxʹ[i]);
						YᵣᵢʹYᵣᵢ = YᵣʹYᵣ.get(zi_old_value).minus(yyʹ[i]);
						XᵣᵢʹYᵣᵢ = XᵣʹYᵣ.get(zi_old_value).minus(xyʹ[i]);
						Ψ_A0ʹΩⁿA0_YᵣᵢʹYᵣᵢ = Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.get(zi_old_value).minus(yyʹ[i]);
						ΩⁿA0_XᵣᵢʹYᵣᵢ = ΩⁿA0_XᵣʹYᵣ.get(zi_old_value).minus(xyʹ[i]);
						Ωⁿ_XᵣᵢʹXᵣᵢ = Ωⁿ_XᵣʹXᵣ.get(zi_old_value).minus(xxʹ[i]);
						Ωᵥᵢ = Ωⁿ_XᵣᵢʹXᵣᵢ.inverse();
						Ψᵥᵢ = Ψ_A0ʹΩⁿA0_YᵣᵢʹYᵣᵢ.minus(ΩⁿA0_XᵣᵢʹYᵣᵢ.transposeMult(Ωᵥᵢ.mult(ΩⁿA0_XᵣᵢʹYᵣᵢ)));
						ldet_Ψᵥᵢ = 0.5*(κ+Nᵣᵢ)*Math.log(Ψᵥᵢ.det());
						ldet_Ωᵥᵢ = 0.5*p*Math.log(Ωᵥᵢ.det());
					}
				} else {
					// Copy the new values without the current point into the local variables,
					// but also add appropriate alternates *with* the current point - see below
					int[] zᵣ = z.which(zi_old_value).value();
					Nᵣᵢ = zᵣ.length;
					if (Nᵣᵢ > 0) {
						Double2D Xᵣᵢ = X.getAll(zᵣ);
						XᵣᵢʹXᵣᵢ = Xᵣᵢ.transposeMult(Xᵣᵢ);
						Double2D Yᵣᵢ = Y.getAll(zᵣ);
						YᵣᵢʹYᵣᵢ = Yᵣᵢ.transposeMult(Yᵣᵢ); 
						XᵣᵢʹYᵣᵢ = Xᵣᵢ.transposeMult(Yᵣᵢ); 
						Ψ_A0ʹΩⁿA0_YᵣᵢʹYᵣᵢ = Ψ_A0ʹΩⁿA0.plus(YᵣᵢʹYᵣᵢ);
						ΩⁿA0_XᵣᵢʹYᵣᵢ = ΩⁿA0.plus(XᵣᵢʹYᵣᵢ);
						Ωⁿ_XᵣᵢʹXᵣᵢ = Ωⁿ.plus(XᵣᵢʹXᵣᵢ);
						Ωᵥᵢ = Ωⁿ_XᵣᵢʹXᵣᵢ.inverse();
						Ψᵥᵢ = Ψ_A0ʹΩⁿA0_YᵣᵢʹYᵣᵢ.minus(ΩⁿA0_XᵣᵢʹYᵣᵢ.transposeMult(Ωᵥᵢ.mult(ΩⁿA0_XᵣᵢʹYᵣᵢ)));
						ldet_Ψᵥᵢ = 0.5*(κ+Nᵣᵢ)*Math.log(Ψᵥᵢ.det());
						ldet_Ωᵥᵢ = 0.5*p*Math.log(Ωᵥᵢ.det());
						
						Nᵣ.put(zi_old_value,Nᵣᵢ+1);
						XᵣʹXᵣ.put(zi_old_value, XᵣᵢʹXᵣᵢ.plus(xxʹ[i]));
						YᵣʹYᵣ.put(zi_old_value, YᵣᵢʹYᵣᵢ.plus(yyʹ[i]));
						XᵣʹYᵣ.put(zi_old_value, XᵣᵢʹYᵣᵢ.plus(xyʹ[i]));
						Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.put(zi_old_value, Ψ_A0ʹΩⁿA0_YᵣᵢʹYᵣᵢ.plus(yyʹ[i]));
						ΩⁿA0_XᵣʹYᵣ.put(zi_old_value, ΩⁿA0_XᵣᵢʹYᵣᵢ.plus(xyʹ[i]));
						Ωⁿ_XᵣʹXᵣ.put(zi_old_value, Ωⁿ_XᵣᵢʹXᵣᵢ.plus(xxʹ[i]));
					} else {
						// FIXME: Will this ever be used?
						Nᵣ.put(zi_old_value,1);
						XᵣʹXᵣ.put(zi_old_value, xxʹ[i]);
						YᵣʹYᵣ.put(zi_old_value, yyʹ[i]); 
						XᵣʹYᵣ.put(zi_old_value, xyʹ[i]); 
						Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.put(zi_old_value, Ψ_A0ʹΩⁿA0.plus(yyʹ[i]));
						ΩⁿA0_XᵣʹYᵣ.put(zi_old_value, ΩⁿA0.plus(xyʹ[i]));
						Ωⁿ_XᵣʹXᵣ.put(zi_old_value, Ωⁿ.plus(xxʹ[i]));
					}
					Ωᵥ.put(zi_old_value, Ωⁿ_XᵣʹXᵣ.get(zi_old_value).inverse());
					Ψᵥ.put(zi_old_value, Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.get(zi_old_value).minus(ΩⁿA0_XᵣʹYᵣ.get(zi_old_value).transposeMult(Ωᵥ.get(zi_old_value).mult(ΩⁿA0_XᵣʹYᵣ.get(zi_old_value)))));
					ldet_Ψᵥ.put(zi_old_value, 0.5*(κ+Nᵣ.get(zi_old_value))*Math.log(Ψᵥ.get(zi_old_value).det()));
					ldet_Ωᵥ.put(zi_old_value, 0.5*p*Math.log(Ωᵥ.get(zi_old_value).det()));							
				}
				
				// Existing categories
				double[] logP = new double[R.size() + 1];
				double maxLogP = -Double.MAX_VALUE;
				int rj = 0;
				for (rj=0; rj<R.size(); rj++) {
					int r = R.value()[rj];
					int[] zᵣ = z.which(r).value();
					
					if (r == zi_old_value) {
						if (Nᵣᵢ == 0) {
							logP[rj] = -Double.MAX_VALUE;
						} else {
							// Use the -sub_i values as v (i.e., with this point removed)
							// and the values that are still in storage as uv (i.e., that still imagine
							// that this point is in the category)
							// (... this is sort of backwards...)
							double κ_Nᵣ = κ+Nᵣᵢ;
							double κ_Nᵣ_1 = κ_Nᵣ+1;
							double lΓκ_Nᵣ = lΓ((κ_Nᵣ)/2, p);
							double lΓκ_Nᵣ_1 = lΓ((κ_Nᵣ_1)/2, p);
						
							double ldet_Ψᵤᵥ = ldet_Ψᵥ.get(r);
							double ldet_Ωᵤᵥ = ldet_Ωᵥ.get(r);

							logP[rj] = Math.log(Nᵣᵢ) + lΓκ_Nᵣ_1 + ldet_Ωᵤᵥ + ldet_Ψᵥᵢ - lπ - lΓκ_Nᵣ - ldet_Ωᵥᵢ - ldet_Ψᵤᵥ;
						}
					} else {
						if (!XᵣʹXᵣ.containsKey(r)) {
							Nᵣ.put(r, zᵣ.length);
							if (Nᵣ.get(r) > 0) {
								Double2D Xᵣ = X.getAll(zᵣ);
								XᵣʹXᵣ.put(r, Xᵣ.transposeMult(Xᵣ));
								Double2D Yᵣ = Y.getAll(zᵣ);
								YᵣʹYᵣ.put(r, Yᵣ.transposeMult(Yᵣ)); 
								XᵣʹYᵣ.put(r, Xᵣ.transposeMult(Yᵣ)); 
								Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.put(r, Ψ_A0ʹΩⁿA0.plus(YᵣʹYᵣ.get(r)));
								ΩⁿA0_XᵣʹYᵣ.put(r, ΩⁿA0.plus(XᵣʹYᵣ.get(r)));
								Ωⁿ_XᵣʹXᵣ.put(r, Ωⁿ.plus(XᵣʹXᵣ.get(r)));
								Ωᵥ.put(r, Ωⁿ_XᵣʹXᵣ.get(r).inverse());
								Ψᵥ.put(r, Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.get(r).minus(ΩⁿA0_XᵣʹYᵣ.get(r).transposeMult(Ωᵥ.get(r).mult(ΩⁿA0_XᵣʹYᵣ.get(r)))));
								ldet_Ψᵥ.put(r, 0.5*(κ+Nᵣ.get(r))*Math.log(Ψᵥ.get(r).det()));
								ldet_Ωᵥ.put(r, 0.5*p*Math.log(Ωᵥ.get(r).det()));							
							}
	
						}	
						double κ_Nᵣ = κ+Nᵣ.get(r);
						double κ_Nᵣ_1 = κ_Nᵣ+1;
						double lΓκ_Nᵣ = lΓ((κ_Nᵣ)/2, p);
						double lΓκ_Nᵣ_1 = lΓ((κ_Nᵣ_1)/2, p);
						
						Double2D ΩⁿA0_XᵣʹYᵣ_xᵢyᵢʹ = ΩⁿA0_XᵣʹYᵣ.get(r).plus(xyʹ[i]);
						
						Double2D Ωᵤᵥ = Ωⁿ_XᵣʹXᵣ.get(r).plus(xxʹ[i]).inverse();
						Double2D Ψᵤᵥ = Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.get(r).plus(yyʹ[i]).minus(ΩⁿA0_XᵣʹYᵣ_xᵢyᵢʹ.transposeMult(Ωᵤᵥ.mult(ΩⁿA0_XᵣʹYᵣ_xᵢyᵢʹ)));
						
						double ldet_Ψᵤᵥ = 0.5*κ_Nᵣ_1*Math.log(Ψᵤᵥ.det());
						double ldet_Ωᵤᵥ = 0.5*p*Math.log(Ωᵤᵥ.det());
						
						logP[rj] = Math.log(Nᵣ.get(r)) + lΓκ_Nᵣ_1 + ldet_Ωᵤᵥ + ldet_Ψᵥ.get(r) - lπ - lΓκ_Nᵣ - ldet_Ωᵥ.get(r) - ldet_Ψᵤᵥ;	
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
				if (logP[rj] == Double.POSITIVE_INFINITY) {
					logP[rj] = Double.MAX_VALUE;
				}				
				if (logP[rj] > maxLogP) {
					maxLogP = logP[rj];
				}
				
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
				if (zi_new_value != zi_old_value) {
					// update new category
					if (zi_new_rjindex == newc_rjindex) {
						Nᵣ.put(zi_new_value, 1);
						XᵣʹXᵣ.put(zi_new_value, xxʹ[i]);
						YᵣʹYᵣ.put(zi_new_value, yyʹ[i]); 
						XᵣʹYᵣ.put(zi_new_value, xyʹ[i]); 
						Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.put(zi_new_value, Ψ_A0ʹΩⁿA0.plus(yyʹ[i]));
						ΩⁿA0_XᵣʹYᵣ.put(zi_new_value, ΩⁿA0.plus(xyʹ[i]));
						Ωⁿ_XᵣʹXᵣ.put(zi_new_value, Ωⁿ.plus(xxʹ[i]));
						Ωᵥ.put(zi_new_value, Ωⁿ_XᵣʹXᵣ.get(zi_new_value).inverse());
						Ψᵥ.put(zi_new_value, Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.get(zi_new_value).minus(ΩⁿA0_XᵣʹYᵣ.get(zi_new_value).transposeMult(Ωᵥ.get(zi_new_value).mult(ΩⁿA0_XᵣʹYᵣ.get(zi_new_value)))));
						ldet_Ψᵥ.put(zi_new_value, 0.5*(κ+Nᵣ.get(zi_new_value))*Math.log(Ψᵥ.get(zi_new_value).det()));
						ldet_Ωᵥ.put(zi_new_value, 0.5*p*Math.log(Ωᵥ.get(zi_new_value).det()));						
					} else {
						Nᵣ.put(zi_new_value,Nᵣ.get(zi_new_value)+1);
						XᵣʹXᵣ.put(zi_new_value, XᵣʹXᵣ.get(zi_new_value).plus(xxʹ[i]));
						YᵣʹYᵣ.put(zi_new_value, YᵣʹYᵣ.get(zi_new_value).plus(yyʹ[i]));
						XᵣʹYᵣ.put(zi_new_value, XᵣʹYᵣ.get(zi_new_value).plus(xyʹ[i]));
						Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.put(zi_new_value, Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.get(zi_new_value).plus(yyʹ[i]));
						ΩⁿA0_XᵣʹYᵣ.put(zi_new_value, ΩⁿA0_XᵣʹYᵣ.get(zi_new_value).plus(xyʹ[i]));
						Ωⁿ_XᵣʹXᵣ.put(zi_new_value, Ωⁿ_XᵣʹXᵣ.get(zi_new_value).plus(xxʹ[i]));
						Ωᵥ.put(zi_new_value, Ωⁿ_XᵣʹXᵣ.get(zi_new_value).inverse());
						Ψᵥ.put(zi_new_value, Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.get(zi_new_value).minus(ΩⁿA0_XᵣʹYᵣ.get(zi_new_value).transposeMult(Ωᵥ.get(zi_new_value).mult(ΩⁿA0_XᵣʹYᵣ.get(zi_new_value)))));
						ldet_Ψᵥ.put(zi_new_value, 0.5*(κ+Nᵣ.get(zi_new_value))*Math.log(Ψᵥ.get(zi_new_value).det()));
						ldet_Ωᵥ.put(zi_new_value, 0.5*p*Math.log(Ωᵥ.get(zi_new_value).det()));						
					}
					// update old category
					if (z.which(zi_old_value).size() == 1) {
						Nᵣ.remove(zi_old_value);
						XᵣʹXᵣ.remove(zi_old_value);
						YᵣʹYᵣ.remove(zi_old_value);
						XᵣʹYᵣ.remove(zi_old_value);
						Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.remove(zi_old_value);
						ΩⁿA0_XᵣʹYᵣ.remove(zi_old_value);
						Ωⁿ_XᵣʹXᵣ.remove(zi_old_value);
						Ωᵥ.remove(zi_old_value);
						Ψᵥ.remove(zi_old_value);
						ldet_Ψᵥ.remove(zi_old_value);
						ldet_Ωᵥ.remove(zi_old_value);
					} else {
						Nᵣ.put(zi_old_value,Nᵣᵢ);
						XᵣʹXᵣ.put(zi_old_value, XᵣᵢʹXᵣᵢ);
						YᵣʹYᵣ.put(zi_old_value, YᵣᵢʹYᵣᵢ);
						XᵣʹYᵣ.put(zi_old_value, XᵣᵢʹYᵣᵢ);
						Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.put(zi_old_value, Ψ_A0ʹΩⁿA0_YᵣᵢʹYᵣᵢ);
						ΩⁿA0_XᵣʹYᵣ.put(zi_old_value, ΩⁿA0_XᵣᵢʹYᵣᵢ);
						Ωⁿ_XᵣʹXᵣ.put(zi_old_value, Ωⁿ_XᵣᵢʹXᵣᵢ);
						Ωᵥ.put(zi_old_value, Ωᵥᵢ);
						Ψᵥ.put(zi_old_value, Ψᵥᵢ);
						ldet_Ψᵥ.put(zi_old_value, ldet_Ψᵥᵢ);
						ldet_Ωᵥ.put(zi_old_value, ldet_Ωᵥᵢ);
					}
				}
				z.set(i, zi_new_value);
			}
			return z;
		}		
	}
				
	@SuppressWarnings("unchecked")
	public static ChainLink pointEstimate(Chain c, Map<String, Object> extended_hypers) {
		double max = Double.NEGATIVE_INFINITY;
		int argmax = c.size() - 1;
		
		GammaDist gammaDist = new GammaDist((Double0D) extended_hypers.get("alpha_a"), (Double0D) extended_hypers.get("alpha_b"));	
		
		double[] prior_lil = ((Double1D) extended_hypers.get("prior_lil")).value();
		int N = prior_lil.length;
		
		double κ = ((Double0D) extended_hypers.get("κ")).value();
		double lπ = ((Double0D) extended_hypers.get("lπ")).value();
		
		Double2D Α0 = (Double2D) extended_hypers.get("A0");
		Double2D Ωⁿ = (Double2D) extended_hypers.get("Ωⁿ");
		Double2D ΩⁿA0 = (Double2D) extended_hypers.get("ΩⁿA0");
		Double2D Ψ_A0ʹΩⁿA0 = (Double2D) extended_hypers.get("Ψ_A0ʹΩⁿA0");
		int p = Α0.numCols();
		
		Double2D[] xxʹ = ((Double2D[]) extended_hypers.get("xxʹ"));
		Double2D[] yyʹ = ((Double2D[]) extended_hypers.get("yyʹ"));
		Double2D[] xyʹ = ((Double2D[]) extended_hypers.get("xyʹ"));
		
		Double2D X = ((Double2D) extended_hypers.get("X"));
		Double2D Y = ((Double2D) extended_hypers.get("Y"));
		
		int[] previous_z = null;
		
		Map<Integer,Double>[] log_cumul_lil = (Map<Integer,Double>[]) Array.newInstance((Class<HashMap<Integer,Double>>)(Class<?>)HashMap.class, N);
		Map<Integer,Double>[] cumul_lil = (Map<Integer,Double>[]) Array.newInstance((Class<HashMap<Integer,Double>>)(Class<?>)HashMap.class, N);
		
		for (int m = 0; m < c.size(); m++) {
			double posterior = 0;
			ChainLink cl_m = c.get(m);
			
			// α
			Double0D α = ((RandomVar<Double0D>) cl_m.get("alpha")).getNumericValue();
			posterior += gammaDist.logDensity(α);
			
			// z
			double lα = Math.log(α.value());
			int[] z = ((RandomVar<Integer1D>) cl_m.get("z")).getNumericValue().value();
			Map<Integer,Object> cat_changed = new HashMap<Integer,Object>();
			for (int i = 0; i < N; i++) {
				if (m == 0) {
					log_cumul_lil[i] = new HashMap<Integer,Double>();
					log_cumul_lil[i] = new HashMap<Integer,Double>();
				}
				
				int rᵢ = z[i];
				
				if (i > 0) {
					int[] zᵢ_ints = new int[i];
					System.arraycopy(z,0,zᵢ_ints,0,i); // FIXME - yeah, this is dumb - fix it later
					Integer1D zᵢ = new Integer1D(zᵢ_ints);
					Integer1D Rᵢ = zᵢ.items();
					int[] zᵣᵢ = zᵢ.which(rᵢ).value();
					int Nᵣᵢ = zᵣᵢ.length;
		
					if (cat_changed.containsKey(rᵢ)) {
						// Category has changed since last sample
						// The likelihood integral for this point needs
						// to be recomputed
						if (Nᵣᵢ > 0) {
							// Category still has points before i
							Double2D Xᵣᵢ = X.getAll(zᵣᵢ);
							Double2D Yᵣᵢ = Y.getAll(zᵣᵢ);
							Double2D XᵣᵢʹXᵣᵢ = Xᵣᵢ.transposeMult(Xᵣᵢ);
							Double2D YᵣᵢʹYᵣᵢ = Yᵣᵢ.transposeMult(Yᵣᵢ); 
							Double2D XᵣᵢʹYᵣᵢ = Xᵣᵢ.transposeMult(Yᵣᵢ);
							
							double κ_Nᵣᵢ = κ+Nᵣᵢ;
							double κ_Nᵣᵢ_1 = κ_Nᵣᵢ+1;
							
							Double2D Ωⁿ_XᵣᵢʹXᵣᵢ = Ωⁿ.plus(XᵣᵢʹXᵣᵢ);
							Double2D Ψ_A0ʹΩⁿA0_YᵣᵢʹYᵣᵢ = Ψ_A0ʹΩⁿA0.plus(YᵣᵢʹYᵣᵢ);
							Double2D ΩⁿA0_XᵣᵢʹYᵣᵢ = ΩⁿA0.plus(XᵣᵢʹYᵣᵢ);
							Double2D ΩⁿA0_XᵣᵢʹYᵣᵢ_xᵢyᵢʹ = ΩⁿA0_XᵣᵢʹYᵣᵢ.plus(xyʹ[i]);
							
							Double2D Ωᵥ = Ωⁿ_XᵣᵢʹXᵣᵢ.inverse();
							Double2D Ωᵤᵥ = Ωⁿ_XᵣᵢʹXᵣᵢ.plus(xxʹ[i]).inverse();
							Double2D Ψᵥ = Ψ_A0ʹΩⁿA0_YᵣᵢʹYᵣᵢ.minus(ΩⁿA0_XᵣᵢʹYᵣᵢ.transposeMult(Ωᵥ.mult(ΩⁿA0_XᵣᵢʹYᵣᵢ)));
							Double2D Ψᵤᵥ = Ψ_A0ʹΩⁿA0_YᵣᵢʹYᵣᵢ.plus(yyʹ[i]).minus(ΩⁿA0_XᵣᵢʹYᵣᵢ_xᵢyᵢʹ.transposeMult(Ωᵤᵥ.mult(ΩⁿA0_XᵣᵢʹYᵣᵢ_xᵢyᵢʹ)));
							
							double ldet_Ωᵥ = 0.5*p*Math.log(Ωᵥ.det());
							double ldet_Ωᵤᵥ = 0.5*p*Math.log(Ωᵤᵥ.det());
							double ldet_Ψᵥ = 0.5*κ_Nᵣᵢ*Math.log(Ψᵥ.det());
							double ldet_Ψᵤᵥ = 0.5*κ_Nᵣᵢ_1*Math.log(Ψᵤᵥ.det());	
	
							double lΓκ_Nᵣᵢ = lΓ((κ_Nᵣᵢ)/2, p);
							double lΓκ_Nᵣᵢ_1 = lΓ((κ_Nᵣᵢ_1)/2, p);
							
							double lil = Math.log(Nᵣᵢ) + lΓκ_Nᵣᵢ_1 + ldet_Ωᵤᵥ + ldet_Ψᵥ - lπ - lΓκ_Nᵣᵢ - ldet_Ωᵥ - ldet_Ψᵤᵥ;
							
							log_cumul_lil[i].put(rᵢ, lil);
							cumul_lil[i].put(rᵢ, Math.exp(lil));						
						} else {
							// Category now has no points before i: use the "new category" term
							double lil = lα + prior_lil[i];
							log_cumul_lil[i].put(rᵢ, lil);
							cumul_lil[i].put(rᵢ, Math.exp(lil));
						}
					}					
					// Compute normalizing constant
					double nc = cumul_lil[i].get(rᵢ);
					for (int j_r=0; j_r<Rᵢ.size(); j_r++) {
						int r = Rᵢ.value()[j_r];
						if (r != rᵢ) {
							if (cat_changed.containsKey(r)) {
								// Another category has changed: recompute posterior integrated likelihood
								int[] zᵣ = zᵢ.which(r).value();
								int Nᵣ = zᵣ.length;
					
								Double2D Xᵣ = X.getAll(zᵣ);
								Double2D Yᵣ = Y.getAll(zᵣ);
								Double2D XᵣʹXᵣ = Xᵣ.transposeMult(Xᵣ);
								Double2D YᵣʹYᵣ = Yᵣ.transposeMult(Yᵣ); 
								Double2D XᵣʹYᵣ = Xᵣ.transposeMult(Yᵣ); 
								
								double κ_Nᵣ = κ+Nᵣ;
								double κ_Nᵣ_1 = κ_Nᵣ+1;
								
								Double2D Ωⁿ_XᵣʹXᵣ = Ωⁿ.plus(XᵣʹXᵣ);
								Double2D Ψ_A0ʹΩⁿA0_YᵣʹYᵣ = Ψ_A0ʹΩⁿA0.plus(YᵣʹYᵣ);
								Double2D ΩⁿA0_XᵣʹYᵣ = ΩⁿA0.plus(XᵣʹYᵣ);
								Double2D ΩⁿA0_XᵣʹYᵣ_xᵢyᵢʹ = ΩⁿA0_XᵣʹYᵣ.plus(xyʹ[i]);
								
								Double2D Ωᵥ = Ωⁿ_XᵣʹXᵣ.inverse();
								Double2D Ωᵤᵥ = Ωⁿ_XᵣʹXᵣ.plus(xxʹ[i]).inverse();
								Double2D Ψᵥ = Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.minus(ΩⁿA0_XᵣʹYᵣ.transposeMult(Ωᵥ.mult(ΩⁿA0_XᵣʹYᵣ)));
								Double2D Ψᵤᵥ = Ψ_A0ʹΩⁿA0_YᵣʹYᵣ.plus(yyʹ[i]).minus(ΩⁿA0_XᵣʹYᵣ_xᵢyᵢʹ.transposeMult(Ωᵤᵥ.mult(ΩⁿA0_XᵣʹYᵣ_xᵢyᵢʹ)));
								
								double ldet_Ωᵥ = 0.5*p*Math.log(Ωᵥ.det());
								double ldet_Ωᵤᵥ = 0.5*p*Math.log(Ωᵤᵥ.det());
								double ldet_Ψᵥ = 0.5*κ_Nᵣ*Math.log(Ψᵥ.det());
								double ldet_Ψᵤᵥ = 0.5*κ_Nᵣ_1*Math.log(Ψᵤᵥ.det());	
		
								double lΓκ_Nᵣ = lΓ((κ_Nᵣ)/2, p);
								double lΓκ_Nᵣ_1 = lΓ((κ_Nᵣ_1)/2, p);
								
								double lil = Math.log(Nᵣ) + lΓκ_Nᵣ_1 + ldet_Ωᵤᵥ + ldet_Ψᵥ - lπ - lΓκ_Nᵣ - ldet_Ωᵥ - ldet_Ψᵤᵥ;
								log_cumul_lil[i].put(r, lil);
								cumul_lil[i].put(r, Math.exp(lil));							
							}
							nc += cumul_lil[i].get(r);
						}
					}
					if (Nᵣᵢ > 0) {
						nc += Math.exp(lα + prior_lil[i]);
					}
					posterior += log_cumul_lil[i].get(rᵢ) - Math.log(nc);
				} else { // i == 0: probability 1
					posterior += 0;
				}

				if (m == 0) {
					cat_changed.put(rᵢ,null);
				} else {
					if (!cat_changed.containsKey(rᵢ)) {
						// Up to i-1, category had not changed since last sample...
						if (rᵢ != previous_z[i]) {
							// ... but category (and old counterpart) just changed: notify future points
							cat_changed.put(rᵢ,null);
							cat_changed.put(previous_z[i],null);
						}
					}					
				}
			}
			previous_z = z;
			
			// Set max
			if (posterior > max) {
				argmax = m;
				max = posterior;
			}
			
		}

		return c.get(argmax);
	}

	public MLM_noparams() {
		super();
	}

	public MLM_noparams(Map<String, Flattenable> hypers, Map<String, Flattenable> init, Double2D data) {
		super(hypers);
		
		// Set up parameters
		this.params = new HashMap<String, RandomVar<? extends Flattenable>>();
		
		// Extend hyperparameters
		this.extended_hypers = MLM_noparams.extend_hypers(hypers, data, (Double2D) init.get("X"));

		// Z
		CRPDist postZ = new CRPDist();
		Integer1D defaultZ = (Integer1D) init.get("z");
		RandomVar<Integer1D> rvZ = new RandomVar<Integer1D>("z", postZ, defaultZ);
		this.params.put("z", rvZ);

		// alpha
		ProbDistMC<Double0D> postAl = new PostAlphaDist();
		Double0D defaultAl = (Double0D) init.get("alpha");
		RandomVar<Double0D> rvAl = new RandomVar<Double0D>("alpha", postAl, defaultAl);
		this.params.put("alpha", rvAl);
	}

	@Override
	public ChainLink getInitialLink() {
		ChainLink cl = new ChainLink();
		cl.add(this.params.get("z"));
		cl.add(this.params.get("alpha"));
		return cl;
	}
}
