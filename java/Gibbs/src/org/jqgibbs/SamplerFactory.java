package org.jqgibbs;

import java.util.List;

import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Numeric;

public class SamplerFactory {
	public static Sampler getSampler(Double2D d, Model m) throws
			GibbsException {
		return new GenericSampler(m, d);
	}
}
