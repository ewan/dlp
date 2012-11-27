package org.jqgibbs.models;

import java.util.ArrayList;
import java.util.List;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDist;
import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;

import cern.jet.random.Gamma;
import cern.jet.random.engine.RandomEngine;

public class Deterministic<T extends Numeric> extends ProbDist<T> {
	private Numeric answer;
	
	public Deterministic(Numeric answer) throws ProbDistParmException {
		this.answer = answer;
	}
	
	private Numeric getAnswer() {
		return answer;
	}

	@Override
	protected double getDensity(T pt) {
		if (pt.equals(this.getAnswer())) {
			return 1;
		} else {
			return 0;
		}
	}	
}
