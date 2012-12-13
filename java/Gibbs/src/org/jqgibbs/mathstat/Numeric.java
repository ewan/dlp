package org.jqgibbs.mathstat;

public abstract interface Numeric extends Cloneable {
	public Double1D rowVec();
	public int length1D();
	public Object clone() throws CloneNotSupportedException;
	public <T extends Numeric> T cloneFromVector(Double1D d);
}
