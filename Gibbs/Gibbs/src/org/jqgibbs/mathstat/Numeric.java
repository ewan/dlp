package org.jqgibbs.mathstat;

/*
 * TODO
 * Couldn't we get rid of a lot of copy-paste code? Or is it all stuck there
 * because of primitives?
 */
public abstract class Numeric<T extends Numeric<T>> implements Cloneable {
	@SuppressWarnings("unchecked")
	public AbstractSequence<? extends AbstractSequence<?, T>, T> sequence() {
		return new ListSequence<T>((T) this);
	}
	
	@Override
	public abstract Object clone() throws CloneNotSupportedException;

	@Override
	public abstract boolean equals(Object o);
	
	@Override
	public abstract int hashCode();
	
	public abstract T cloneFromVector(Double1D v);

	// FIXME
	// You need to get your dimensions straight and just
	// make this vec().
	public abstract Double1D rowVec();
	
	public abstract int length1D();
}
