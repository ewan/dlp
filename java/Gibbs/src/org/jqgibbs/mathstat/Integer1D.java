package org.jqgibbs.mathstat;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;

public class Integer1D implements Numeric, Iterable<Integer0D> {

	private int[] ints;
	private boolean dirtyWhich;
	private Map<Integer0D,Integer1D> which;

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

	private synchronized void setInts(int[] ints) {
		this.ints = ints;
		this.dirtyWhich = true;
	}

	private synchronized int[] getInts() {
		return this.ints;
	}

	public boolean add(Integer0D i0D) {
		int[] ints = Arrays.copyOf(this.getInts(), this.size() + 1);
		if (i0D == null) {
			throw new NullPointerException();
		}
		ints[this.size()] = i0D.value();
		this.setInts(ints);
		this.dirtyWhich = true;		
		return true;
	}

	public boolean addAll(Collection<? extends Integer0D> c) {
		int[] ints = Arrays.copyOf(this.getInts(), this.size() + c.size());
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

	public boolean contains(Object o) {
		int d;
		if (o == null) {
			throw new NullPointerException();
		} else if (o instanceof Integer0D) {
			d = ((Integer0D) o).value();
		} else if (o instanceof Integer) {
			d = ((Integer) o).intValue();
		} else {
			return false;
		}
		for (int i = 0; i < this.size(); i++) {
			if (this.get(i).value() == d) {
				return true;
			}
		}
		return false;
	}

	public int size() {
		return this.getInts().length;
	}

	public int[] value() {
		return this.getInts();
	}

	//@Override
	public Integer0D get(int i) {
		return new Integer0D(this.getInts()[i]);
	}
	
	//@Override
	public Integer1D getAll(int... dis) {
		int[] d = this.value();
		int[] all = new int[dis.length];
		int ai = 0;
		for (int di : dis) {
			all[ai] = d[di]; 
			ai++;
		}
		return new Integer1D(all);
	}
	
	//@Override
	public Object clone() throws CloneNotSupportedException {
		return new Integer1D(this.getInts().clone());
	}

	/**
	 * Returns a new Integer1D object that contains exactly those indices
	 * in this Integer1D with value equal to some particular number.
	                          
	@param  n  the integer to find against

	 */

	//@Override
	public Integer1D which(Integer0D n) {
		if (this.dirtyWhich) {
			Integer1D active = this.items();
			this.which = new HashMap<Integer0D,Integer1D>();
			for (Integer0D k : active) {
				int[] w = new int[this.size()];
				int wi = 0;
				for (int i=0; i<this.size(); i++) {
					if (this.value()[i] == k.value()) {
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

	public Integer0D max() {
		int max = Integer.MIN_VALUE;
		for (Integer0D i : this) {
			if (i.value() > max) {
				max = i.value();
			}
		}
		return new Integer0D(max);
	}

	//@Override
	public Integer1D cloneFromVector(Double1D v) {
		int[] ds = new int[this.size()];
		for (int i=0; i<this.size(); i++) {
			if (i >= v.size()) {
				return new Integer1D(ds);
			}
			ds[i] = (int) v.get(i).value();
		}
		return new Integer1D(ds);
	}
	
	//@Override
	public boolean equals(Object o) {
		if (!(o instanceof Integer1D)) {
			return false;
		}		
		return this.value() == ((Integer1D) o).value();
	}
	
	//@Override
	public int hashCode() {
		int h = 0;
		for (int i=0; i<this.size(); i++) {
			h = h + this.value()[i];
		}
		return h;		
	}

	public Integer1D intersect(Integer1D other) {
		int[] i1, i2;
		if (this.size() >= other.size()) {
			i1 = Arrays.copyOf(other.value(), other.size());
			i2 = Arrays.copyOf(this.value(), this.size());
		} else {
			i1 = Arrays.copyOf(this.value(), this.size());
			i2 = Arrays.copyOf(other.value(), other.size());		
		}
		Arrays.sort(i1);
		Arrays.sort(i2);
		int[] intersection = new int[Math.min(i1.length, i2.length)];
		int n=0;		
		if (i1.length > 0) {
			int last = i1[0]-1;
			for (int i=0; i<i1.length; i++) {
				if (i1[i] != last) {
					if (Arrays.binarySearch(i2, i1[i]) >= 0) {
						intersection[n] = i1[i];
						n++;
					}
					last = i1[i];
				}
			}
		}
		return new Integer1D(Arrays.copyOf(intersection, n));
	}

	//@Override
	public boolean remove(Object o) {
		if (!(o instanceof Integer0D)) {
			throw new UnsupportedOperationException();
		}
		int i = ((Integer0D) o).value();
		return this.remove(i);
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
	
	//@Override
	public Integer0D set(int i, Integer0D t) throws IndexOutOfBoundsException {
		if (i > this.size() || i < 0) {
			throw new IndexOutOfBoundsException("Tried to set index out of range"); // FIXME
		}
		Integer0D oldValue = null;
		if (i == this.size()) {
			this.add(t);
		} else {
			oldValue = this.get(i);
			this.value()[i] = t.value();
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
			this.which.get(t).add(new Integer0D(i));
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

	public Integer0D minNotIn(int a, int b) {
		int[] values = Arrays.copyOf(this.value(), this.size());
		Arrays.sort(values);
		for (int v=a; v<=b; v++) {
			if (Arrays.binarySearch(values, v) < 0) {
				return new Integer0D(v);
			}
		}
		return new Integer0D(b+1);
	}
	

	//@Override
	public String toString() {
		String s = "";
		String prefix = "";
		for (int i=0; i < this.size(); i++) {
			s = s + prefix + String.valueOf(this.value()[i]);
			prefix = " ";
		}
		return s;
	}

	public Double1D toDouble1D() {
		double[] ds = new double[this.size()];
		for (int i=0; i<this.size(); i++) {
			ds[i] = this.ints[i];
		}
		return new Double1D(ds);
	}
	
	public Double1D mult(Double2D o) {
		return this.toDouble1D().mult(o);
	}
	
	//@Override
	public Double1D rowVec() {
		return this.toDouble1D();
	}
	
	//@Override
	public int length1D() {
		return this.size();
	}

	public int getValue(int i) {
		return this.value()[i];
	}

	public Integer1D samePairs() {
		int n = this.size();
		int nPairs = n*n - n*(n+1)/2;
		int[] v = this.getInts();
		int[] same = new int[nPairs];
		int si = 0;
		for (int i=0; i<n; i++) {
			for (int j=(i+1); j<n; j++) {
				if (v[i] == v[j]) {
					same[si] = 1;
				}
				si++;
			}
		}
		return new Integer1D(same);
	}
	
	public Integer0D mode() {
		int max = 0;
		Integer1D items = this.items();
		Integer0D mode = items.get(0);
		for (Integer0D item : items) {
			int size = this.getAll(this.which(item).value()).size();
			if (size > max) {
				mode = item;
				max = size;
			}
		}
		return mode;
	}

	public Integer1D getAllBut(int j) {
		int[] d = this.value();
		int[] all = new int[d.length-1];
		for (int i=0; i<this.size(); i++) {
			if (i < j) {
				all[i] = d[i];
			} else if (i > j) {
				all[i-1] = d[i];
			}
		}
		return new Integer1D(all);	
	}
	
	public Integer2D outer(Integer1D o) {
		double[] ot = new double[this.size()];
		for (int i=0; i<this.ints.length; i++) {
			ot[i] = this.ints[i];
		}
		double[] od = new double[o.size()];
		int[] ov = o.value();
		for (int i=0; i<ov.length; i++) {
			od[i] = ov[i];
		}
		DoubleMatrix2D r = AlgebraStatic.multOuterStatic(new DenseDoubleMatrix1D(ot),
														 new DenseDoubleMatrix1D(od));
		int[][] ri = new int[this.size()][o.size()];
		double[][] rd = r.toArray();
		for (int i=0; i<this.ints.length; i++) {
			for (int j=0; j<ov.length; j++) {
				ri[i][j] = (int) rd[i][j];
			}
		}
		return new Integer2D(ri);
	}

	public Integer0D sum() {
		int s = 0;
		for (int i=0; i<this.ints.length; i++) {
			s += this.ints[i];
		}
		return new Integer0D(s);
	}	
	
	public Iterator<Integer0D> iterator() {
		return new Iterator<Integer0D>() {
			private int curr = 0;

			public boolean hasNext() {
				return (this.curr < Integer1D.this.size());
			}

			public Integer0D next() {
				Integer0D d = Integer1D.this.get(this.curr);
				this.curr++;
				return d;
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}
}
