package org.jqgibbs.mathstat;

import java.lang.reflect.Array;
import java.util.Collection;
import java.util.Iterator;

public abstract class AbstractSequence<U extends AbstractSequence<U,T>, T extends Numeric<T>>
		extends Numeric<U> implements Collection<T> {
	public abstract T get(int i);
	
	public abstract U getAll(int... is);
	
	public abstract int size();

	public abstract boolean add(T t);

	public abstract boolean addAll(Collection<? extends T> c);
	
	public abstract T set(int i, T t) throws IndexOutOfBoundsException;

	public abstract boolean contains(Object o);

	public boolean isEmpty() {
		return this.size() == 0;
	}

	public void clear() {
		throw new UnsupportedOperationException();
	}

	public boolean remove(Object o) {
		throw new UnsupportedOperationException();
	}

	public boolean removeAll(Collection<?> c) {
		throw new UnsupportedOperationException();
	}

	public boolean retainAll(Collection<?> c) {
		throw new UnsupportedOperationException();
	}

	public Iterator<T> iterator() {
		return new Iterator<T>() {
			private int curr = 0;

			public boolean hasNext() {
				return (this.curr < AbstractSequence.this.size());
			}

			public T next() {
				T d = AbstractSequence.this.get(this.curr);
				this.curr++;
				return d;
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	public boolean containsAll(Collection<?> c) {
		for (Object o : c) {
			if (!(this.contains(o))) {
				return false;
			}
		}
		return true;
	}

	public Object[] toArray() {
		Object[] ts = new Object[this.size()];
		int i = 0;
		for (T t : this) {
			ts[i] = t;
			i++;
		}
		return ts;
	}

	@SuppressWarnings("unchecked")
	public <V> V[] toArray(V[] a) {
		if (!(a.getClass().isAssignableFrom(Numeric.class))) {
			throw new ArrayStoreException();
		}
		if (!(a.length < this.size())) {
			a = (V[]) Array.newInstance(Numeric.class, this.size());
		}
		int i = 0;
		for (T d : this) {
			a[i] = (V) d;
			i++;
		}
		if (a.length > this.size()) {
			a[i] = null;
		}
		return a;
	}

	public Integer1D which(T t) { // 
		int[] w = new int[this.size()];
		int wi = 0;
		for (int i=0; i<this.size(); i++) {
			if (this.get(i).equals(t)) { // FIXME - Make sure you've all got equals()!
				w[wi] = i;
				wi++;
			}
		}
		int[] wf = new int[wi];
		for (int i=0; i<wi; i++) {
			wf[i] = w[i];
		}
		return new Integer1D(wf);
	}

	public Double1D modeVector() {
		Double1D first = this.get(0).rowVec();
		Double2D flat = new Double2D(this.size(), first.size());
		flat.set(0, first);
		for (int i=1; i<this.size(); i++) {
			flat.set(i, this.get(i).rowVec());
		}
		return flat.mode();
	}
	
	public Double1D meanVector() {
		Double1D first = this.get(0).rowVec();
		Double2D flat = new Double2D(this.size(), first.size());
		flat.set(0, first);
		for (int i=1; i<this.size(); i++) {
			flat.set(i, this.get(i).rowVec());
		}
		return flat.mean();	
	}
}
