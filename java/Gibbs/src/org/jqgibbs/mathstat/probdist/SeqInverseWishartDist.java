package org.jqgibbs.mathstat.probdist;

import java.util.List;

import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Double3D;
import org.jqgibbs.mathstat.Numeric;

public class SeqInverseWishartDist extends
		SequentialProbDist<Double3D, Double2D, InverseWishartDist> {

	public SeqInverseWishartDist(Numeric<?>... parms)
			throws ProbDistParmException {
		super(parms);
		List<Numeric<?>[]> settingsAll = this.getSettingsAll();
		assert settingsAll.size() > 0; // FIXME - right now this is *not* a safe
										// assertion!!
		// FIXME - Maybe we should be able to make 0-args calls to PDs without exceptions
		this.setProbDist(new InverseWishartDist(settingsAll.get(0)));			
	}

	protected Double3D genVariate() throws ProbDistParmException {
		Double3D sequence = null;
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
