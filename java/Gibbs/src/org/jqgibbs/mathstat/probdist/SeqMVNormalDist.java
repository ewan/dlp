//package org.jqgibbs.mathstat.probdist;
//
//import java.util.List;
//
//import org.jqgibbs.mathstat.Double1D;
//import org.jqgibbs.mathstat.Double2D;
//import org.jqgibbs.mathstat.Numeric;
//
//public class SeqMVNormalDist extends SequentialProbDist<Double2D, Double1D, MVNormalDist> {
//	
//	public SeqMVNormalDist(Numeric... parms) throws ProbDistParmException {
//		super(parms);
//		List<Numeric[]> settingsAll = this.getSettingsAll();
//		assert settingsAll.size() > 0; // FIXME - right now this is *not* a safe
//										// assertion!!
//		// FIXME - Maybe we should be able to make 0-args calls to PDs without exceptions
//		this.setProbDist(new MVNormalDist(settingsAll.get(0)));	
//	}
//	
//	protected Double2D genVariate() throws ProbDistParmException {
//		List<Numeric[]> s = this.getSettingsAll();
//		int n = s.size();
//		double[][] sequence = null;
//		int i = 0;
//		for (Numeric[] settings : s) {
//			double[] v = this.getProbDist().variate(settings).value();			
//			if (sequence == null) {
//				sequence = new double[n][v.length];
//			}
//			for (int j=0; j<v.length; j++) {
//				sequence[i][j] = v[j];
//			}
//			i++;
//		}
//		return new Double2D(sequence);
//	}
//
//}
