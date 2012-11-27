package org.jqgibbs.mathstat.probdist;

import org.jqgibbs.mathstat.Numeric;


public abstract class ProbDistParmCheck {
	public abstract boolean test(Numeric parms);
	public abstract String message();
}