package org.jqgibbs;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Integer1D;

/**
 * TODO document with an eye to moving towards coda and towards memory efficiency
 * 
 * @author ewan
 *
 */
public class Chain implements Iterable<ChainLink> {

	private boolean flattened;
	private Map<String,Integer1D> flatOm;
	private Double2D flatDouble2D;

	private List<ChainLink> c;

	public Chain() {
		this.flattened = false;
		this.setC(new LinkedList<ChainLink>());
	}

	public void addLink(ChainLink l) {
		this.c.add(l);
		this.flattened = false;
	}

	public Iterator<ChainLink> iterator() {
		return this.c.iterator();
	}

	public int size() {
		return this.c.size();
	}
	
	public ChainLink get(int i) {
		return this.c.get(i);
	}
	
	private void setC(List<ChainLink> c) {
		this.c = c;
		this.flattened = false;
	}

	private void flatten() {
		if (this.flattened) {
			return;
		} else {
			// Find longest link
			int max = -1;
			int argmax = -1;
			for (int i=0; i<this.size(); i++) {
				ChainLink cl = this.get(i);
				int length = cl.length1D();
				if (length > max) {
					max = length; 
					argmax = i;
				}
			}
			// Get the order map for the longest link
			this.flatOm = this.get(argmax).getOrderMap();
			// Create flattened structure
			int rows = this.size();
			int cols = max;
			double[][] f2d = new double[rows][cols];
			// Fill in the matrix
			for (int i=0; i<this.size(); i++) {
				double[] row = this.get(i).paddedVec(this.flatOm).value();
				for (int j=0; j<row.length; j++) {
					f2d[i][j] = row[j];
				}
			}
			// Set flattened
			this.flatDouble2D = new Double2D(f2d);			
			this.flattened = true;
		}	
	}
	
	public Map<String,Integer1D> getOrderMap() {
		this.flatten();
		return this.flatOm;
	}
	
	public Double2D getFlatDouble2D() {
		this.flatten();
		return this.flatDouble2D;
	}
}
