package org.jqgibbs.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jqgibbs.Chain;
import org.jqgibbs.ChainLink;
import org.jqgibbs.ChainWriter;
import org.jqgibbs.Flattenable;
import org.jqgibbs.GenericSampler;
import org.jqgibbs.Model;
import org.jqgibbs.RandomEngineSelector;
import org.jqgibbs.Sampler;
import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Double3D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Integer1D;
import org.jqgibbs.mathstat.probdist.CategoricalDist;
import org.jqgibbs.models.MLM_sample_params;
import org.jqgibbs.models.MLM_sample_params_varbsel_block;

import au.com.bytecode.opencsv.CSVReader;

public class GibbsEwan {

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

	private static String dataFileName = "/Users/emd/Work/active_projects/dlp_project/dlp/java/Gibbs/spanish_mfiau_f0.txt";
//	private static String dataFileName = "/Users/ewan/Work/School/Projects/DLEN/Data/hillsmall.txt";
	private static int maxIter = 500;
	private static int burnIn = 0;
	private static int lag = 5;
	private static String outFileName = "dump.out";
	private static String pteFileName = "pte.out";

	public static Double2D readDataFile(String fileName) throws IOException {
		CSVReader cr = new CSVReader(new BufferedReader(new InputStreamReader(
				new FileInputStream(new File(fileName)))));
		List<String[]> rows = cr.readAll();
		cr.close();
		int nrow = rows.size() - 1;
		int ncol = Arrays.asList(rows.get(0)).size();
		double[][] valueDouble = new double[nrow][ncol];
		for (int i = 1; i <= nrow; i++) {
			double[] rowDouble = GibbsEwan.parseDoubleArray(rows.get(i));
			valueDouble[i - 1] = rowDouble;
		}
		return new Double2D(valueDouble);
	}

	public static void main(String[] args) throws Exception {
		RandomEngineSelector.setFixedSeed(true);
		
		Double2D d = GibbsEwan.readDataFile(GibbsEwan.dataFileName);
//		int h = 3;
		int h = 4;
		// predictors
		CategoricalDist catDist = new CategoricalDist(new Double1D(new double[] {0.5,0.5}));
		double[][] fb = new double[d.size()][h];
		for (int i = 0; i < d.size(); i++) {
			fb[i][0] = 1;
//			fb[i][1] = (double) d.getCol(3).get(i).value();
//			fb[i][2] = catDist.variate().value();
			fb[i][1] = d.getCol(4).get(i).value();
			fb[i][2] = d.getCol(0).get(i).value();
			fb[i][3] = catDist.variate().value();
		}
		Double2D X = new Double2D(fb);
		// responses
		// int ddim = 3;
		Double2D Y = d.getColAll(1, 2, 3);
//		Double2D Y = d.getColAll(0, 1, 2);
		int p = 3;
		
		// hyperparameters
		Map<String, Flattenable> hypers = new HashMap<String, Flattenable>();
		Double2D cov = Y.cov();
		Double2D XʹX = X.transposeMult(X);
		
		hypers.put("deadline", new Integer0D(1));
		hypers.put("T0z", new Double0D(1));
		hypers.put("Tfz", new Double0D(1));
//		hypers.put("beta", new Double0D(1));
		hypers.put("W", new Double2D(h, p));
		hypers.put("S", cov);
//		double[][] psi = new double[ddim][ddim];
//		for (int i = 0; i < ddim; i++) {
//			for (int j = 0; j < ddim; j++) {
//				if (i == j) {
//					psi[i][j] = 10000;
//					// } else {
//					// psi[i][j] = 0.01;
//				}
//			}
//		}
//		hypers.put("Psi", cov);
//		hypers.put("kappa", new Double0D(4));
		double[][] phi = new double[h][h];
//		phi[0][0] = 1;
		phi[0][0] = 1;
		for (int i = 1; i < h; i++) {
			phi[i][i] = 0.05;
//			phi[i][i] = 9E-4;
		}
//		Double2D Φ = Double2D.ident(h).plus(XʹX).inverse();
		Double2D Φ = new Double2D(phi).mult(1E-0);
//		Double2D Φ = Double2D.ident(h).plus(XʹX).inverse().mult(1E-3);
//		Double2D Φ = Double2D.ident(h).mult(10000);
		hypers.put("Phi", Φ);
//		hypers.put("lambda", new Double0D(22));
		hypers.put("lambda", new Double0D(22));
//		hypers.put("be", new Double0D(1));
//		hypers.put("ga", new Double0D(1));
//
		
//		Double2D Ω = Double2D.ident(h).plus(XʹX).inverse().mult(0.5);
		Double2D Ω = Double2D.ident(h);
//		Double2D Ω = Double2D.ident(h).mult(5000);
//		Double2D A0 = Ω.mult(X.transposeMult(Y));
		Double2D A0 = new Double2D(h,p);
		hypers.put("Omega", Ω);
		hypers.put("A0", A0);
		hypers.put("Psi", cov);
//		hypers.put("Psi", cov);
		hypers.put("kappa", new Double0D(4));
//		hypers.put("kappa", new Double0D(7E1));
//		hypers.put("alpha_a", new Double0D(20000));
		hypers.put("alpha_a", new Double0D(2));
		hypers.put("alpha_b", new Double0D(.5));
//		hypers.put("alpha_b", new Double0D(0.001));
		hypers.put("tau",new Double0D(0.5));

		// initialization
		Map<String, Flattenable> init = new HashMap<String, Flattenable>();
		init.put("z", new Integer1D(new int[d.size()]));
		init.put("X", X);
		init.put("alpha", new Double0D(1));
		init.put("A0", A0);
		init.put("A", new Double3D(new Double2D(h, p)));
//		double[] M = new double[ddim];
//		for (int i = 0; i < ddim; i++) {
//			M[i] = 0;
//		}
		double[][] sg = new double[p][p];
		for (int i = 0; i < p; i++) {
			sg[i][i] = 1;
		}
		init.put("Sigma", new Double3D(new Double2D(sg)));
//		double[][] omega = new double[h][h];
//		for (int i = 0; i < h; i++) {
//			omega[i][i] = 1;
//		}
		init.put("Omega", Ω);
		int[] gamma = new int[h];
		gamma[0] = 1;
		for (int i = 1; i < h; i++) {
			gamma[i] = 1;
		}
		Integer1D γ = new Integer1D(gamma);
		init.put("gamma", γ);
		
		
		// run gibbs sampler
		Model m = new MLM_sample_params(hypers, init, Y); 
		//Model m = new MLM_sample_params_varbsel_block(hypers, init, Y); 
		//Model m = new MLM_noparams(hypers, init, Y); 
		Sampler s = new GenericSampler(m);
		
		System.err.println(System.currentTimeMillis());
		
		for (int i = 0; i < burnIn; i++) {
			if (i % 25 == 0) {
				System.err.println("Burnin iteration " + String.valueOf(i + 1));
			}
			s.variate();
		}
		Chain c = new Chain();
		ChainLink cl = s.variate();
		for (int i = 0; i < maxIter; i++) {
			cl = s.variate();
			c.addLink(cl);
			System.err.println(String.valueOf(i + 1));
			for (int j = 0; j < lag; j++) {
				s.variate();
			}
		}
		
		System.err.println(System.currentTimeMillis());
		
		// finish up
		ChainWriter cw = new ChainWriter(GibbsEwan.outFileName, m);
		cw.write(c);
		// ChainLink pte = MLM_sample_params.pointEstimate(c, d);
		ChainLink pte = MLM_sample_params_varbsel_block.pointEstimate(c, MLM_sample_params_varbsel_block.extend_hypers(hypers, Y, X));
		
		System.err.println(System.currentTimeMillis());
		
		c = new Chain();
		c.addLink(pte);
		cw = new ChainWriter(GibbsEwan.pteFileName, m);
		cw.write(c);
	}
}
