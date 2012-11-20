package org.jqgibbs.mathstat.probdist;

import java.util.ArrayList;
import java.util.List;

import org.jqgibbs.mathstat.AbstractSequence;
import org.jqgibbs.mathstat.Numeric;

public abstract class SequentialProbDist<U extends AbstractSequence<U, T>, T extends Numeric<T>, V extends ProbDist<T>>
		extends ProbDist<U> {

	private List<Numeric<?>[]> settingsAll;

	protected V probDist;
	protected void setProbDist(V v) {
		this.probDist = v;
	}
	protected V getProbDist() {
		return this.probDist;
	}

	protected List<Numeric<?>[]> getSettingsAll() {
		return this.settingsAll;
	}

	public SequentialProbDist(Numeric<?>... parms) throws ProbDistParmException {
		super(parms);		
	}

	@Override
	protected void installParmChecks() {
		return;
	}

	@Override
	protected void initializeParms(Numeric<?>... parms)
			throws ProbDistParmException {
		boolean allSequences = true;
		for (Numeric<?> p : parms) {
			if (!(p instanceof AbstractSequence<?, ?>)) {
				allSequences = false;
				break;
			}
		}
		if (allSequences) {
			AbstractSequence<?, ?> sequence;
			// FIXME - check
			int seqLength = ((AbstractSequence<?, ?>) parms[0]).size();
			this.settingsAll = new ArrayList<Numeric<?>[]>(seqLength);
			for (int i = 0; i < parms.length; i++) {
				sequence = (AbstractSequence<?, ?>) parms[i];
				if (sequence.size() != seqLength) {
					throw new ProbDistParmException(
							"Received parameter lists of different lengths");
				}
				for (int j = 0; j < seqLength; j++) {
					if (this.settingsAll.size() < j+1) {
						assert this.settingsAll.size() == j;
						this.settingsAll.add(new Numeric<?>[parms.length]);
					}
					this.settingsAll.get(j)[i] = sequence.get(j);
				}
			}
		} else {
			// Didn't get a sequence - just set up a sequence with length 1
			this.settingsAll = new ArrayList<Numeric<?>[]>(1);
			this.getSettingsAll().add(parms);
		}
		this.initialized = true;
	}

	@Override
	protected void initializeParmsInitial(Numeric<?>... parms)
			throws ProbDistParmException {
		this.initializeParms(parms);
	}

	@Override
	protected double getDensity(U pt) throws ProbDistParmException {
		throw new UnsupportedOperationException("Too lazy, come back later");
	}
}
