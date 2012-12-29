package org.jqgibbs;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Integer1D;

public class ChainLink implements Iterable<RandomVar<?>> {

	private Map<String, RandomVar<?>> m = new HashMap<String, RandomVar<?>>();
	private List<String> order = new LinkedList<String>();
	
	public ChainLink() { }

	public RandomVar<?> get(String name) {
		return this.m.get(name);
	}
	
	public boolean add(RandomVar<?> gv) {
		if (this.m.containsKey(gv.getName())) {
			return false;
		}
		this.m.put(gv.getName(), gv);
		return this.getOrder().add(gv.getName());
	}

	public String toString() {
		String s = "";
		String prefix = "";
		List<String> names = this.getOrder();
		for (String name : names) {
			s += prefix + this.get(name).toString();
			prefix = " ";
		}
		return s;
	}

	public Iterator<RandomVar<?>> iterator() {
		return new Iterator<RandomVar<?>>() {
			private Iterator<String> i = ChainLink.this.getOrder().iterator();
			public boolean hasNext() {
				return this.i.hasNext();
			}
			public RandomVar<?> next() {
				return ChainLink.this.m.get(this.i.next());
			}
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	public int size() {
		assert this.m.size() == this.order.size();
		return this.order.size();
	}

	public List<String> getOrder() {
		return order;
	}

	public int length1D() {
		int n = 0;
		for (RandomVar<?> rv : this) {
			n += rv.getNumericValue().length1D();
		}		
		return n;
	}
	
	public Map<String,Integer1D> getOrderMap() {
		Map<String,Integer1D> om = new HashMap<String,Integer1D>();
		int k=0;
		for (RandomVar<?> rv : this) {
			int length = rv.getNumericValue().length1D();
			int[] is = new int[length];
			for (int i=0; i<length; i++) {
				is[i] = k;
				k++;
			}
			Integer1D indices = new Integer1D(is);
			om.put(rv.getName(),indices);
		}
		return om;
	}

	public Double1D paddedVec(Map<String, Integer1D> flatOm) {
		int length = 0;
		for (String n : flatOm.keySet()) {
			Integer1D ind = flatOm.get(n);
			length += ind.size();
		}
		double[] pv = new double[length];
		int actualLength = 0;
		for (String n : flatOm.keySet()) {
			int firstIndex = flatOm.get(n).value()[0];			
			RandomVar<?> v = this.get(n);
			double[] fv = v.getNumericValue().rowVec().value();
			if (actualLength + fv.length > length) {
				throw new IllegalArgumentException("Provided inappropriate order map");				
			} else {
				for (int j=0; j<fv.length; j++) {
					pv[firstIndex+j] = fv[j];
				}	
				actualLength += fv.length;				
			}
		}
		return new Double1D(pv);
	}
}
