package org.jqgibbs.mathstat.probdist;

import java.util.List;

import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Integer1D;
import org.jqgibbs.mathstat.Numeric;

public class SeqCategoricalDist extends
		SequentialProbDist<Integer1D, Integer0D, CategoricalDistInitializeByP> {

	public SeqCategoricalDist(Numeric<?>... params)
			throws ProbDistParmException {
		super(params);
		List<Numeric<?>[]> settingsAll = this.getSettingsAll();
		assert settingsAll.size() > 0; // FIXME - right now this is *not* a safe
										// assertion!!
		// FIXME - Maybe we should be able to make 0-args calls to PDs without exceptions
		this.setProbDist(new CategoricalDistInitializeByP(settingsAll.get(0)));		
	}

	protected Integer1D genVariate() throws ProbDistParmException {
		assert this.getProbDist() != null;
		Integer1D sequence = null;
		for (Numeric<?>[] settings : this.getSettingsAll()) {
			if (sequence == null) {
				sequence = this.getProbDist().variate(settings).sequence();
			} else {
				sequence.add(this.getProbDist().variate(settings));
			}
		}
		return sequence;
	}
}
