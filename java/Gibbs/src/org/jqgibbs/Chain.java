package org.jqgibbs;

import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;

import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Integer1D;

public class Chain  implements List<ChainLink> {

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

	public ChainLink first() {
		return this.c.get(0);
	}

	public ChainLink last() {
		return this.c.get(this.c.size() - 1);
	}

	public ChainLink current() {
		return this.last();
	}

	@Override
	public String toString() {
		String s = "";
		String prefix = "";
		synchronized (this) {
			for (int i = 0; i < this.c.size(); i++) {
//				System.out.println(String.valueOf(i));
				s += prefix + Integer.toString(i + 1) + ": " + this.c.get(i);
				prefix = "\n";
			}
		}
		return s;
	}

	public boolean isEmpty() {
		return this.c.isEmpty();
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

	protected List<ChainLink> getC() {
		return this.c;
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

	public boolean add(ChainLink arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public void add(int arg0, ChainLink arg1) {
		// TODO Auto-generated method stub
		
	}

	public boolean addAll(Collection<? extends ChainLink> arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public boolean addAll(int arg0, Collection<? extends ChainLink> arg1) {
		// TODO Auto-generated method stub
		return false;
	}

	public void clear() {
		// TODO Auto-generated method stub
		
	}

	public boolean contains(Object arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public boolean containsAll(Collection<?> arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public int indexOf(Object arg0) {
		// TODO Auto-generated method stub
		return 0;
	}

	public int lastIndexOf(Object arg0) {
		// TODO Auto-generated method stub
		return 0;
	}

	public ListIterator<ChainLink> listIterator() {
		// TODO Auto-generated method stub
		return null;
	}

	public ListIterator<ChainLink> listIterator(int arg0) {
		// TODO Auto-generated method stub
		return null;
	}

	public boolean remove(Object arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public ChainLink remove(int arg0) {
		// TODO Auto-generated method stub
		return null;
	}

	public boolean removeAll(Collection<?> arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public boolean retainAll(Collection<?> arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public ChainLink set(int arg0, ChainLink arg1) {
		// TODO Auto-generated method stub
		return null;
	}

	public List<ChainLink> subList(int arg0, int arg1) {
		// TODO Auto-generated method stub
		return null;
	}

	public Object[] toArray() {
		// TODO Auto-generated method stub
		return null;
	}

	public <T> T[] toArray(T[] arg0) {
		// TODO Auto-generated method stub
		return null;
	}
}
