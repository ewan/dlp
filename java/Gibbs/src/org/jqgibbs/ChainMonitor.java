package org.jqgibbs;

import java.util.Observer;

public abstract class ChainMonitor implements Observer {
	public abstract boolean done();
	public abstract boolean fatal();
}
