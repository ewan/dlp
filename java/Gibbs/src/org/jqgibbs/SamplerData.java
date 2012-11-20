package org.jqgibbs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

import au.com.bytecode.opencsv.*;

import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.AbstractSequence;
import org.jqgibbs.mathstat.Numeric;

/*
 * TODO
 * Should we really wrap the IOException?
 * Should this extend Numeric?
 * 
 * This class seems like a *mistake*. Just use Double2D and set up a reader class.
 */
public class SamplerData extends AbstractSequence<SamplerData, Double1D> {

	private List<String> colNames;
	private Double2D numericValue;

	public SamplerData(Double2D numericValue) {
		this.setNumericValue(numericValue);
	}

	public SamplerData(Double2D numericValue, List<String> colNames) {
		this.setNumericValue(numericValue);
		this.setColNames(colNames);
	}

	public SamplerData(String dataFileName) throws GibbsException {
		this(new File(dataFileName));
	}

	/*
	 * TODO Expects header FIXME Seems to be inserting a row zeroes at the head
	 */
	public SamplerData(File dataFile) throws DataIOException {
		try {
			CSVReader cr = new CSVReader(new BufferedReader(
					new InputStreamReader(new FileInputStream(dataFile))));
			List<String[]> rows = cr.readAll();
			cr.close();
			this.setColNames(Arrays.asList(rows.get(0)));
			int nrow = rows.size() - 1;
			int ncol = this.getColNames().size();
			double[][] valueDouble = new double[nrow][ncol];
			for (int i = 1; i <= nrow; i++) {
				double[] rowDouble = Gibbs.parseDoubleArray(rows.get(i));
				valueDouble[i-1] = rowDouble;
			}
			this.setNumericValue(new Double2D(valueDouble));
		} catch (IOException e) {
			throw new DataIOException("Unable to read from data file: "
					+ e.getMessage());
		}
	}

	public List<String> getColNames() {
		return this.colNames;
	}

	private void setColNames(List<String> colNames) {
		this.colNames = colNames;
	}

	public int getNumDataPoints() {
		return this.getNumericValue().numRows();
	}

	private void setNumericValue(Double2D numericValue) {
		this.numericValue = numericValue;
	}

	public Double2D getNumericValue() {
		return numericValue;
	}

	@Override
	public Object clone() throws CloneNotSupportedException {
		return new SamplerData(this.getNumericValue(), this.getColNames());
	}

	@Override
	public Double1D get(int i) {
		return this.getNumericValue().get(i);
	}

	@Override
	public SamplerData getAll(int... dis) {
		return new SamplerData(this.getNumericValue().getAll(dis), this
				.getColNames());
	}

	@Override
	public boolean add(Double1D e) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean addAll(Collection<? extends Double1D> c) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean contains(Object o) {
		return this.getNumericValue().contains(o);
	}

	@Override
	public Iterator<Double1D> iterator() {
		return this.getNumericValue().iterator();
	}

	@Override
	public int size() {
		return this.getNumericValue().size();
	}

	public Double2D cov() {
		return this.getNumericValue().cov();
	}
	
	public Double1D mean() {
		return this.getNumericValue().mean();
	}

	public Double1D sum() {
		return this.getNumericValue().sum();
	}

	@Override
	public SamplerData cloneFromVector(Double1D v) {
		return new SamplerData(this.getNumericValue()
				.cloneFromVector(v));
	}

	public int dim() {
		// FIXME
		return this.get(0).size();
	}

	@Override
	public boolean equals(Object o) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public int hashCode() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public Double1D set(int i, Double1D t) throws IndexOutOfBoundsException {
		// TODO Auto-generated method stub
		return null;
	}

	public void removeCol(int col) {
		//FIXME checks
		double[][] c = new double[this.size()][this.getNumericValue().numCols()-1];
		for (int i=0; i<this.size(); i++) {
			int jC = 0;
			for (int j=0; j<this.getNumericValue().numCols(); j++) {
				if (j != col) {
					c[i][jC] = this.getNumericValue().value()[i][j];
					jC++;
				}
			}
		}
		this.setNumericValue(new Double2D(c));
	}

	@Override
	public int length1D() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public Double1D rowVec() {
		// TODO Auto-generated method stub
		return null;
	}

}
