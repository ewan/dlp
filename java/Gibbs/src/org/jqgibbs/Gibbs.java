package org.jqgibbs;

import java.util.HashMap;
import java.util.Map;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Double3D;
import org.jqgibbs.mathstat.Integer1D;
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.models.FLGFDModel;

public class Gibbs {

	public static double[] parseDoubleArray(String[] ss) {
		double[] ssd = new double[ss.length];
		for (int i = 0; i < ss.length; i++) {
			try { 
				ssd[i] = Double.parseDouble(ss[i]);
			} catch (NumberFormatException e) {
				ssd[i] = 0;
			}
		}
		return ssd;
	}

	private static String dataFileName = "/Users/jsf/Desktop/eclipse/gibbs/dlp/java/Gibbs/spanish_mfiau_f0.txt";
	private static int maxIter = 50;//1000;
	private static int burnIn = 10;//500;
	private static int lag = 5;
	private static String outFileName = "dump.out";

	public static void main(String[] args) throws Exception {
		Double2D d;
		d = new SamplerData(Gibbs.dataFileName).getNumericValue();
		int h = 3;
		// predictors
		double[][] fb = new double[d.size()][h];	
		for (int i=0; i<d.size(); i++) {
			fb[i][0] = 1;
			fb[i][1] = (double) d.getCol(4).get(i).value();
			fb[i][2] = (double) d.getCol(0).get(i).value();
		}
		Double2D fixedB = new Double2D(fb);
		// responses
		int ddim = 3;
		d = d.getColAll(1,2,3);
		// hyperparameters
		Double2D cov = d.cov();
		Map<String, Numeric> hypers = new HashMap<String,Numeric>();
		hypers.put("lambda", new Double0D(22));
		hypers.put("beta", new Double0D(1));
		hypers.put("W", new Double2D(h,ddim));
		hypers.put("S", cov);
		double[][] psi = new double[ddim][ddim];
		for (int i=0; i<ddim; i++) {
			for (int j=0; j<ddim; j++) {
				if (i==j) {
					psi[i][j] = 10000;	
//				} else {
//					psi[i][j] = 0.01;
				}
			}
		}
		hypers.put("Psi", cov);
		hypers.put("kappa", new Double0D(4));
		double[][] phi = new double[h][h];
		phi[0][0] = 1;
		for (int i=1; i<h; i++) {
			phi[i][i] = 0.05;
		}
		hypers.put("Phi", new Double2D(phi));
		hypers.put("ala", new Double0D(2));
		hypers.put("alb", new Double0D(.5));
		hypers.put("xa", new Double0D(1));
		hypers.put("xb", new Double0D(1));			
		hypers.put("be", new Double0D(1));
		hypers.put("ga", new Double0D(1));			
		
		// initialization
		Map<String, Numeric> init = new HashMap<String,Numeric>();
		init.put("M", new Double2D(h,ddim));
		init.put("A", new Double3D(new Double2D(h,ddim)));
		double[] M = new double[ddim];
		for (int i=0; i<ddim; i++) {
			M[i] = 0;
		}
		double[][] sg = new double[ddim][ddim];
		for (int i=0; i<ddim; i++) {
			sg[i][i] = 1;
		}			
		init.put("Sg", new Double3D(new Double2D(sg)));
		double[][] omega = new double[h][h];
		for (int i=0; i<h; i++) {
			omega[i][i] = 1;
		}			
		init.put("Omega", new Double2D(omega));
		init.put("M", new Double2D(h,ddim));
		init.put("Z", new Integer1D(new int[d.size()]));
		init.put("B", fixedB);
		init.put("x", new Double0D(0.5));
		init.put("al", new Double0D(1));			
		int[] gamma = new int[2];
		for (int i=0; i<2; i++) {
			gamma[i] = 1;
		}
		// run gibbs sampler
		Sampler s = SamplerFactory.getSampler(d, new FLGFDModel(hypers, init, d.numCols()));
		for (int i=0; i<burnIn; i++) {
			if (i % 25 == 0) {
				System.err.println("Burnin iteration " + String.valueOf(i+1));
			}
			s.variateFast();
		}
		Chain c = new Chain();
		ChainLink cl = s.variateFast();
		for (int i=0; i<maxIter; i++) {
			cl = s.variateFast();
			c.addLink(cl);
			System.err.println(String.valueOf(i+1));
			for (int j=0; j<lag; j++) {
				s.variateFast();
			}
			//System.out.println(cl.toString());
		}
		// finish up
		ChainWriter cw = new ChainWriter(Gibbs.outFileName, s.getModel());
		cw.write(c);
		ChainLink pte = FLGFDModel.pointEstimate(c, d);
		c = new Chain();
		c.addLink(pte);
	}
}
