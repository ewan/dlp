package org.jqgibbs;

import java.lang.management.ManagementFactory;
import java.util.Observable;

public class MemoryMonitor extends ChainMonitor {
	private double p;
	private int checkEvery;
	private boolean fatal;
	
	private int lastCheck;
	private boolean memoryFull;
	private long heapSize;
	private long memoryUsage;
	
	public MemoryMonitor(double p, int checkEvery, boolean fatal) {
		this.p = p;
		this.checkEvery = checkEvery;
		this.fatal = fatal;
		this.lastCheck = -checkEvery;
		this.memoryFull = false;
		this.heapSize = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage().getCommitted(); 
	}

	public void update(Observable obs, Object obj) {
		Chain c = (Chain) obs;
		if (c.size() - this.lastCheck >= this.checkEvery) {
			this.memoryUsage = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage().getUsed();
			if (this.memoryUsage > p*this.heapSize) {
				this.memoryFull = true;
			}
			this.lastCheck = c.size();
		}
	}
	
	@Override
	public boolean done() {
		return this.memoryFull;
	}

	@Override
	public boolean fatal() {
		return this.fatal;
	}

}
