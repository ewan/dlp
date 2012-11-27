package org.jqgibbs;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Integer1D;
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.mathstat.RandomVar;

/*
 * TODO
 * Is it appropriate to have a new Exception subclass here?
 */
public class ChainLink implements Cloneable, Iterable<RandomVar<?>>, Numeric {

	// FIXME - shouldn't this be a ListSequence?
	private Map<String, RandomVar<?>> m;
	private List<String> order;
	
	public ChainLink(RandomVar<?>... c) {
		this.setM(new HashMap<String, RandomVar<?>>());
		this.setOrder(new LinkedList<String>());
		for (RandomVar<?> gv : c) {
			this.add(gv);
		}
	}

	public RandomVar<?> get(String name) {
		return this.getM().get(name);
	}
	
	public boolean add(RandomVar<?> gv) {
		if (this.getM().containsKey(gv.getName())) {
			return false;
		}
		this.getM().put(gv.getName(), gv);
		return this.getOrder().add(gv.getName());
	}
	
	//@Override
	public int hashCode() {
		int h = 1;
		for (RandomVar<?> gv : this.getM().values()) {
			h = h*gv.hashCode();
		}
		return h;
	}
	
	//@Override
	public boolean equals(Object o) {
		if (!(o instanceof ChainLink)) {
			return false;
		} else {
			ChainLink cl = (ChainLink) o;
			return this.getM().equals(cl.getM()) && this.getOrder().equals(cl.getOrder());
		}
	}
	
	//@Override
	public String toString() {
//		String cap = "]";
//		String s = "[";
		String s = "";
		String prefix = "";
		List<String> names = this.getOrder();
		for (String name : names) {
//			s += "\n" + Formatting.indentAll(this.get(name).toString(), 2);
//			cap = "\n]";
			s += prefix + this.get(name).toString();
			prefix = " ";
		}
//		s += cap;
		return s;
	}

	public boolean isEmpty() {
		return this.getM().isEmpty();
	}

	public RandomVar<?> get(int i) {
		return this.getM().get(this.getOrder().get(i));
	}

	public boolean addAll(Collection<? extends RandomVar<?>> c) {
		throw new UnsupportedOperationException();
	}

	public boolean contains(Object o) {
		return this.getOrder().contains(o);
	}

	public Iterator<RandomVar<?>> iterator() {
		return new Iterator<RandomVar<?>>() {
			private Iterator<String> i = ChainLink.this.getOrder().iterator();
			public boolean hasNext() {
				return this.i.hasNext();
			}
			public RandomVar<?> next() {
				return ChainLink.this.getM().get(this.i.next());
			}
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	public int size() {
		assert this.getM().size() == this.getOrder().size();
		return this.getOrder().size();
	}

	private void setM(Map<String, RandomVar<?>> m) {
		this.m = m;
	}

	private Map<String, RandomVar<?>> getM() {
		return m;
	}

	private void setOrder(List<String> order) {
		this.order = order;
	}

	public List<String> getOrder() {
		return order;
	}
	
	public Object clone() throws CloneNotSupportedException {
		ChainLink clone = new ChainLink();
		for (RandomVar<?> rv : this) {
			clone.add((RandomVar<?>) rv.clone());
		}
		return (Object) clone;
	}

	//@Override
	public ChainLink getAll(int... is) {
		// TODO Auto-generated method stub
		return null;
	}

	//@Override
	public ChainLink cloneFromVector(Double1D v) {
		// TODO Auto-generated method stub
		return null;
	}

	//@Override
	public int length1D() {
		int n = 0;
		for (RandomVar<?> rv : this) {
			n += rv.length1D();
		}		
		return n;
	}
	
	//@Override
	public Double1D rowVec() {
		List<Double1D> ld = new LinkedList<Double1D>();
		int ncol = 0;
		for (RandomVar<?> rv : this) {
			Double1D v = rv.rowVec();
			ld.add(v);
			ncol += v.size();
		}
		double dd[] = new double[ncol];
		int k = 0;
		for (Double1D v : ld) {
			double dv[] = v.value();
			for (int i=0; i<dv.length; i++) {
				dd[k] = dv[i];
				k++;
			}
		}
		return new Double1D(dd);
	}
	
	public Map<String,Integer1D> getOrderMap() {
		Map<String,Integer1D> om = new HashMap<String,Integer1D>();
		int k=0;
		for (RandomVar<?> rv : this) {
			int length = rv.length1D();
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

	//@Override
	public RandomVar<?> set(int i, RandomVar<?> t)
			throws IndexOutOfBoundsException {
		// TODO Auto-generated method stub
		return null;
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
			double[] fv = v.rowVec().value();
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
