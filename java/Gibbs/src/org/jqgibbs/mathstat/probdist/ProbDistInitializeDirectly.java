//package org.jqgibbs.mathstat.probdist;
//
//import java.util.List;
//
//import org.jqgibbs.mathstat.Numeric;
//
//public abstract class ProbDist<T extends Numeric> extends ProbDist<T> {
//
//	protected Numeric[] parms;
//	
//	@Override
//	/*
//	 * Classes must be defined such that parameters unspecified in parameter list are
//	 * at the end of getParmNames(), getParmClasses() [as null], and getParmCheck()
//	 * [as null].
//	 */
//	
//	public T variateFast(Numeric... parms) throws ProbDistParmException {
//		return this.genVariate();
//	}
//}
