//package org.jqgibbs.mathstat;
//
//import java.util.Arrays;
//import java.util.Collection;
//import java.util.Iterator;
//import java.util.List;
//import java.util.LinkedList;
//
//public class ListSequence<T extends Numeric<T>> extends AbstractSequence<ListSequence<T>, T> {
//	private static final long serialVersionUID = -3878927904883763658L;
//	
//	private List<T> l;
//	
//	public ListSequence(T... t) {
//		this.l = new LinkedList<T>(Arrays.asList(t));
//	}
//	
//	private ListSequence(LinkedList<T> l) {
//		this.l = l;
//	}
//	
//	@Override
//	public int size() {
//		return this.l.size();
//	}
//	
//	@Override
//	public synchronized boolean add(T t) {
//		return this.l.add(t);
//	}
//	
//	@Override
//	public synchronized T get(int i) {
//		return this.l.get(i);
//	}
//	
//	@Override
//	public boolean contains(Object o) {
//		return this.l.contains(o);
//	}
//	
//	@Override
//	public boolean addAll(Collection<? extends T> c) {
//		return this.l.addAll(c);
//	}
//	
//	@Override
//	public Iterator<T> iterator() {
//		return this.l.iterator();
//	}
//	
//	@SuppressWarnings("unchecked")
//	@Override
//	public Object clone() throws CloneNotSupportedException {
//		ListSequence<T> ls = new ListSequence<T>();
//		for (T t : this) {
//			ls.add((T) t.clone());
//		}
//		return ls;
//	}
//	
//	@Override
//	public ListSequence<T> getAll(int... is) {
//		LinkedList<T> l = new LinkedList<T>();
//		for (int i : is) {
//			l.add(this.get(i));
//		}
//		return new ListSequence<T>(l);
//	}
//
//	@Override
//	public ListSequence<T> cloneFromVector(Double1D v) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//	
//	@SuppressWarnings("unchecked")
//	@Override
//	public boolean equals(Object o) {
//		if (!(this.getClass().isAssignableFrom(o.getClass()))) {
//			return false;
//		}
//		return this.l.equals(((ListSequence<T>) o).l); // FIXME - how the hell do I do this?
//	}
//	
//	@Override
//	public int hashCode() {
//		return this.l.hashCode();
//	}
//
//	@Override
//	public T set(int i, T t) {
//		if (i == this.size()) { // Add
//			this.l.add(t);
//			return t;
//		} else {
//			return this.l.set(i, t);
//		}
//	}
//	
//	@Override
//	public Double1D rowVec() {
//		List<Double1D> ld = new LinkedList<Double1D>();
//		int ncol = 0;
//		for (Numeric<? extends Numeric<?>> n : this) {
//			Double1D v = n.rowVec();
//			ld.add(v);
//			ncol += v.size();
//		}
//		double dd[] = new double[ncol];
//		int k = 0;
//		for (Double1D v : ld) {
//			double dv[] = v.value();
//			for (int i=0; i<dv.length; i++) {
//				dd[k] = dv[i];
//				k++;
//			}
//		}
//		return new Double1D(dd);
//	}
//	
//	@Override
//	public int length1D() {
//		int n=0;
//		for (T t : this) {
//			n += t.length1D();
//		}
//		return n;
//	}
//}
