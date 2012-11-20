package org.jqgibbs.tests;

import static org.junit.Assert.*;
import au.com.bytecode.opencsv.*;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;

import org.jqgibbs.GibbsData;
import org.jqgibbs.GibbsException;
import org.jqgibbs.Sampler;
import org.junit.Test;

public class GibbsDataTest {

	private static String testGibbsDataResourceName = "test_gibbs_data_file_000.csv";
	private static List<int[]> testGibbsDataCells = Arrays.asList(
		new int[] {3, 3},
		new int[] {1, 1},
		new int[] {15, 2},
		new int[] {103, 0}
	);
	private static List<int[]> testGibbsHeaderCells = Arrays.asList(
			new int[] {0, 1},
			new int[] {0, 3}
	);
	
	@Test
	public void testGibbsData() throws GibbsException {
		String testFileName = this.getClass().getResource(GibbsDataTest.testGibbsDataResourceName).getFile();
		new GibbsData(testFileName);
	}
	
	private interface CSVTest {
		public void test(List<String[]> rows, GibbsData gd);
	}
	
	private void testCSV(CSVTest ct) throws GibbsException {
		try {
			String testFileName = this.getClass().getResource(GibbsDataTest.testGibbsDataResourceName).getFile();
			GibbsData gd = new GibbsData(testFileName);	
			CSVReader cr = new CSVReader(new BufferedReader(new InputStreamReader(new FileInputStream(testFileName))));
			List<String[]> rows = cr.readAll();
			ct.test(rows, gd);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

//	@Test
//	public void testNumCols() {
//		CSVTest ct = new CSVTest() {
//			public void test(List<String[]> rows, GibbsData gd) {
//				String[] row = rows.iterator().next();
//				assertEquals("Incorrect number of columns", gd.numCols(), row.length);
//			}
//		};
//		this.testCSV(ct);
//	}
	
	@Test
	public void testHeader() throws GibbsException {
		CSVTest ct = new CSVTest() {
			public void test(List<String[]> rows, GibbsData gd) {
				String[] header = rows.iterator().next();
				for (int[] cell : GibbsDataTest.testGibbsHeaderCells) {
					int col = cell[1];
					assertEquals("Header contained incorrect information", gd.colNames().get(col), header[col]);
//					if (!(gd.get(header[col]).name() == header[col])) {
//						fail("Name of variable indexed by name " + header[col] + " did not match header");
//					}
//					if (!(gd.vars().get(col).name() == header[col])) {
//						fail("Name of variable with numerical index " + String.valueOf(col) + " did not match header");	
//					}
				}
			}
		};	
		this.testCSV(ct);
	}

	@Test
	public void testDataIntegrity() throws GibbsException {
		CSVTest ct = new CSVTest() {
			public void test(List<String[]> rows, GibbsData gd) {
				for (int[] cell : GibbsDataTest.testGibbsDataCells) {
					int row = cell[0];
					int col = cell[1];
					double csvCell = Double.parseDouble(rows.get(row)[col]);
//					if (!(gd.vars().get(col).get(row) == csvCell)) {
//						fail("Variable " + String.valueOf(col) + " had the wrong data in row " + String.valueOf(row));	
//					}
//					if (!(gd.points().get(col).get(row) == csvCell)) {
//						fail("Data point " + String.valueOf(row) + " had the wrong data in column " + String.valueOf(col));	
//					}
					assertEquals("Name of variable with numerical index " + String.valueOf(col) + " did not match header",
								 gd.value().get(row, col),
								 csvCell,
								 Sampler.MIN_VALUE);
				}
			}
		};	
		this.testCSV(ct);
	}
}
