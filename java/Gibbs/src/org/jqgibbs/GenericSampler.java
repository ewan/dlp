package org.jqgibbs;


import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.RandomVar;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;

public class GenericSampler extends Sampler {
	
	public GenericSampler(Model m, Double2D d) throws ProbDistParmException {
		super(m, d);
		// TODO Auto-generated constructor stub
	}
	
	protected ChainLink cloneCurrent() {
		try {
			return (ChainLink) this.current.clone();
		} catch (CloneNotSupportedException e) {
			// FIXME
			throw new RuntimeException(e);
		}		
	}
	
	@Override
	public ChainLink variateFast() {
		RandomVar<?> gv;
		ChainLink newLink;
		newLink = this.cloneCurrent();
		for (String varName : newLink.getOrder()) {
			gv = newLink.get(varName);
			try {
				gv.updatePosteriorFast(newLink, this.getData());
			} catch (Exception e) {
				// FIXME
				System.err.println(e.getMessage());
				StackTraceElement[] st = e.getStackTrace();
				for (int i=0; i<st.length; i++) {
					System.err.println(st[i].toString());
				}				
				throw new RuntimeException(e);
			}
		}
		this.setChainLink(newLink);
		return newLink;		
	}

	@Override
	public ChainLink variate() {
		RandomVar<?> gv;
		ChainLink newLink;
		newLink = this.cloneCurrent();
		for (String varName : newLink.getOrder()) {
			gv = newLink.get(varName);
			try {
				gv.updatePosterior(newLink, this.getData());
			} catch (Exception e) {
				// FIXME
				System.err.println(e.getMessage());
				StackTraceElement[] st = e.getStackTrace();
				for (int i=0; i<st.length; i++) {
					System.err.println(st[i].toString());
				}
				throw new RuntimeException(e);
			}
		}
		this.setChainLink(newLink);
		return newLink;
	}
}
