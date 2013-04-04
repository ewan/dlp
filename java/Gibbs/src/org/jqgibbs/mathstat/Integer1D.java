package org.jqgibbs.mathstat;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.jqgibbs.Flattenable;

public class Integer1D implements Flattenable,Cloneable {

	private int[] ints;
	private boolean dirtyWhich;
	private Map<Integer,Integer1D> which;

	public Integer1D(Integer0D... i0Ds) {
		this.setInts(new int[0]);
		this.addAll(Arrays.asList(i0Ds));
	}

	public Integer1D(Integer[] integers) {
		int[] ints = new int[integers.length];
		for (int i = 0; i < integers.length; i++) {
			if (integers[i] == null) {
				throw new NullPointerException();
			}
			ints[i] = integers[i].intValue();
		}
		this.setInts(ints);
	}

	public Integer1D(int[] ints) {
		this.setInts(ints);
	}

	public Double1D rowVec() {
		double[] ds = new double[this.size()];
		for (int i=0; i<this.size(); i++) {
			ds[i] = this.ints[i];
		}
		return new Double1D(ds);		
	}
	
	public int length1D() {
		return this.size();
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof Integer1D)) {
			return false;
		}		
		return Arrays.equals(this.value(),((Integer1D) o).value());
	}
	
	@Override
	public int hashCode() {
		int h = 0;
		for (int i=0; i<this.size(); i++) {
			h = h + this.value()[i];
		}
		return h;		
	}
	
	@Override
	public String toString() {
		String s = "";
		String prefix = "";
		for (int i=0; i < this.size(); i++) {
			s = s + prefix + String.valueOf(this.value()[i]);
			prefix = " ";
		}
		return s;
	}
	
	@Override
	public Object clone() {
		return new Integer1D(this.value().clone());
	}
	
	private synchronized void setInts(int[] ints) {
		this.ints = ints;
		this.dirtyWhich = true;
	}

	public boolean add(int i0D) {
		int[] ints = Arrays.copyOf(this.ints, this.size() + 1);
		ints[this.size()] = i0D;
		this.setInts(ints);
		this.dirtyWhich = true;		
		return true;
	}

	public boolean addAll(Collection<? extends Integer0D> c) {
		int[] ints = Arrays.copyOf(this.ints, this.size() + c.size());
		int i = this.size();
		for (Integer0D curr : c) {
			if (curr == null) {
				throw new NullPointerException();
			}
			ints[i] = curr.value();
			i = i + 1;
		}
		this.setInts(ints);
		this.dirtyWhich = true;		
		return true;
	}

	public int size() {
		return this.ints.length;
	}

	public int[] value() {
		return this.ints;
	}

	public Integer1D getAll(int[] dis) {
		int[] d = this.value();
		int[] all = new int[dis.length];
		int ai = 0;
		for (int di : dis) {
			all[ai] = d[di]; 
			ai++;
		}
		return new Integer1D(all);
	}

	/**
	 * Returns a new Integer1D object that contains exactly those indices
	 * in this Integer1D with value equal to some particular number.
	                          
	@param  n  the integer to find against

	 */
	// FIXME
	public Integer1D which(Integer n) {
		if (this.dirtyWhich) {
			int[] active = this.items().value();
			this.which = new HashMap<Integer,Integer1D>();
			for (int k : active) {
				int[] w = new int[this.size()];
				int wi = 0;
				for (int i=0; i<this.size(); i++) { // Could swap in some other search
					if (this.value()[i] == k) {
						w[wi] = i;
						wi++;
					}
				}
				this.which.put(k, new Integer1D(Arrays.copyOf(w, wi)));
			}
			this.dirtyWhich = false;
		}
		Integer1D whichN = this.which.get(n);
		if (whichN == null) {
			return new Integer1D();
		}
		return whichN;
	}

	public boolean remove(int i) {
		int[] newValue = new int[this.size()];
		boolean found = false;
		int k=0;
		for (int j=0; j<this.size(); j++) {
			if (!found && this.value()[j] == i) {
				found = true;
			} else {
				newValue[k] = this.value()[j];
				k++;
			}
		}
		if (found) {
			this.ints = Arrays.copyOf(newValue, this.size()-1);
			this.dirtyWhich = true;
			return true;
		}
		return false;
	}
	
	public int set(int i, int t) throws IndexOutOfBoundsException {
		if (i > this.size() || i < 0) {
			throw new IndexOutOfBoundsException("Tried to set index out of range"); // FIXME
		}
		Integer oldValue = null;
		if (i == this.size()) {
			this.add(t);
		} else {
			oldValue = this.ints[i];
			this.value()[i] = t;
		}
		if (this.which != null && !this.dirtyWhich) {
			if (oldValue != null && this.which.containsKey(oldValue)) {
				boolean removed = this.which.get(oldValue).remove(i);
				if (!removed) {
					throw new IllegalStateException("Value wasn't in which table when we said it should be");
				}
				if (this.which.get(oldValue).size() == 0) {
					this.which.remove(oldValue);
				}
			}
			if (!this.which.containsKey(t)) {
				this.which.put(t, new Integer1D());
			}
			this.which.get(t).add(i);
		} else {
			this.dirtyWhich = true;
		}
		return t;
	}
	
	/**
	 * @return
	 * An vector containing only the unique items in the vector.
	 */
	public Integer1D items() {
		if (this.size() == 0) {
			return new Integer1D(new int[0]);
		}
		int[] sorted = Arrays.copyOf(this.value(), this.size());
		Arrays.sort(sorted);
		int iunique = 0;
		int[] unique = new int[this.size()];
		int prev = sorted[0] - 1;
		int N = this.size();
		for (int i=0; i<N; i++) {
			if (sorted[i] != prev) {
				unique[iunique] = sorted[i];
				prev = sorted[i];
				iunique++;		
			}
		}
		return new Integer1D(Arrays.copyOf(unique, iunique));
	}

	public int minNotIn(int a, int b) {
		int[] values = Arrays.copyOf(this.value(), this.size());
		Arrays.sort(values);
		for (int v=a; v<=b; v++) {
			if (Arrays.binarySearch(values, v) < 0) {
				return v;
			}
		}
		return b+1;
	}

}
