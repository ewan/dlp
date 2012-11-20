package org.jqgibbs;

import java.util.List;
import java.util.logging.Logger;

import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.mathstat.probdist.ProbDistParmCheck;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;

//public abstract class Sampler implements Runnable {
public abstract class Sampler {
	public static double MIN_VALUE = Double.MIN_VALUE;
	
	protected Model model;
	private Double2D data;
	protected ChainLink current;
	protected Logger logger;
	protected List<? extends ChainMonitor> cms;
	
	public Sampler(Model m, Double2D d) throws ProbDistParmException {
		this.model = m;
		this.setData(d);
		this.setChainLink(m.getInitialLink());
		//this.installParmChecks();
		//this.initializeParms(parms);
//		this.cms = cms;
//		for (ChainMonitor cm:this.cms) {
//			this.getChain().addObserver(cm);
//		}
//		if (this.getChain().isEmpty()) {
//			this.initializeChain();
//		}
	}
	
	protected void setChainLink(ChainLink initialLink) {
		this.current = initialLink;
	}

	public Model getModel() {
		return this.model;
	}
	
	protected void setData(Double2D d) {
		this.data = d;
	}

	protected Double2D getData() {
		return data;
	}

	public abstract ChainLink variate();

	public ChainLink variateFast() {
		return this.variate();
	}
	
//	protected void setChain(Chain c) {
//		this.c = c;
//	}
//
//	protected Chain getChain() {
//		return c;
//	}
//	
//	protected abstract void setUpFromParms();
//	protected abstract void initializeChain();
//	
//	public abstract void runIteration() throws ProbDistParmException ;
//	
//	protected boolean chainMonitorsDone() {
//		for (ChainMonitor cm:this.cms){
//			boolean result = cm.done();
//			if (cm.fatal() && result == true) {
//				return true;
//			}
//			if (!cm.fatal() && result == false) {
//				return false;
//			}
//		}
//		return true;
//	}
//	
//	public void run() {
//		Thread.currentThread().setPriority(Thread.MIN_PRIORITY);
//		while (!Thread.interrupted() && !this.chainMonitorsDone()) {
//			try {
//				this.runIteration();
//			} catch (ProbDistParmException e) {
//				Thread.currentThread().interrupt();
//			}
//			Thread.yield();
//		}
//	}
//	
//// MAYBE SHOULD BE INHERITABLE?
//	
//	protected abstract List<ProbDistParmCheck[]> getParmCheck();
//
//	protected abstract List<Class<? extends Numeric<?>>> getParmClasses();
//
//	protected abstract List<String> getParmNames();
//	
//// TO BE INHERITED(?)
//	
//	protected abstract void installParmChecks();
//	
//	protected boolean initialized = false;
//	protected Numeric<?>[] parms;
//
//	protected void checkParm(String parmName,
//			Class<? extends Numeric<?>> parmType, Numeric<?> parm,
//			ProbDistParmCheck... cs) throws ProbDistParmException {
//		if (parm == null) {
//			if (parmType == null && cs == null) {
//				return; // OK
//			} else {
//				throw new NullPointerException(
//						"Parameter object cannot be null");
//			}
//		} else if (!parmType.isInstance(parm)) {
//			throw new IllegalArgumentException(
//					"Parameter object of wrong class: expected "
//							+ parmType.toString());
//		}
//		if (cs != null) {
//			for (ProbDistParmCheck c : cs) {
//				if (!c.test(parm)) {
//					throw new ProbDistParmException("Invalid parameter " + parmName
//							+ ": " + c.message());
//				}
//			}
//		}
//	}
//	
//	protected void initializeParms(Numeric<?>... parms)
//	throws ProbDistParmException {;
//		// Call all parameter checks and assignments
//		this.parms = new Numeric<?>[this.getParmCheck().size()];
//		if (parms.length > this.getParmCheck().size()) {
//			throw new IllegalArgumentException("Too many parameters");
//		}
//		for (int i=0; i<parms.length; i++) {
//			this.checkParm(this.getParmNames().get(i), this.getParmClasses().get(i), parms[i], this.getParmCheck().get(i));
//			this.parms[i] = this.getParmClasses().get(i).cast(parms[i]);
//		}
//		// Do remaining work
//		this.setUpFromParms();
//		// Set flag
//		this.initialized = true;
//	}
//
//	
	
}
