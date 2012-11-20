package org.jqgibbs.tests;
import static org.junit.Assert.*;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

import org.jqgibbs.mathstat.rv.RandomVar2D;
import org.junit.Test;


public class GibbsVariate2DTest {

//	@Test
//	public void testSize() {
//		GibbsVariable g = new GibbsVariable((double) 1, (double) 3, (double) 56);
//		DenseDoubleMatrix2D d = new DenseDoubleMatrix2D(new double[][] { new double[] {1,  3, 56} });
//		assertEquals(g.size(), d.size());
//	}

	private static double[][] testMatrix1 = new double[][] {{1,2},{3,4}};
	private static double[][] testMatrix2 = new double[][] {{2,2},{3,4}};	
//	private static double[] testVector1 = new double[] {5,6};
//	private static double[][] testAddRowResult1 = new double[][] {{1,2},{3,4},{5,6}};
	private static String testName1 = "V";
	private static String testName2 = "W";
	private static String testString1 = "V: {1, 2},\n   {3, 4}";
	
	@Test
	public void testGibbsVariable2D() {
		new RandomVar2D(GibbsVariate2DTest.testMatrix1);
	}
	
	@Test
	public void testHashCode() {
		RandomVar2D g1 = new RandomVar2D(
				GibbsVariate2DTest.testName1,
				GibbsVariate2DTest.testMatrix1);
		RandomVar2D g2 = new RandomVar2D(
				GibbsVariate2DTest.testName1,
				GibbsVariate2DTest.testMatrix1);
		assertEquals(g1.hashCode(), g2.hashCode());
	}

	@Test
	public void testEquals() {
		RandomVar2D g1 = new RandomVar2D(
				GibbsVariate2DTest.testName1,
				GibbsVariate2DTest.testMatrix1);
		RandomVar2D g2 = new RandomVar2D(
				GibbsVariate2DTest.testName1,
				GibbsVariate2DTest.testMatrix1);
		RandomVar2D g3 = new RandomVar2D(
				GibbsVariate2DTest.testName2,
				GibbsVariate2DTest.testMatrix1);
		RandomVar2D g4 = new RandomVar2D(
				GibbsVariate2DTest.testName1,
				GibbsVariate2DTest.testMatrix2);
		assertEquals(g1, g2);
		assertFalse(g1.equals(g3));
		assertFalse(g1.equals(g4));
	}
	
	@Test
	public void testValueColt() {
		RandomVar2D g = new RandomVar2D(GibbsVariate2DTest.testMatrix1);
		DenseDoubleMatrix2D d = new DenseDoubleMatrix2D(GibbsVariate2DTest.testMatrix1);
		assertEquals(g.toColt(), d);
	}
	
	@Test
	public void testName() {
		RandomVar2D g = new RandomVar2D(
				GibbsVariate2DTest.testName1,
				GibbsVariate2DTest.testMatrix1);
		assertEquals(g.getName(), GibbsVariate2DTest.testName1);
	}
	
//	@Test
//	public void testAddRow() {
//		GibbsVariable2D g = new GibbsVariable2D(GibbsVariable2DTest.testMatrix1);
//		g.addRow(GibbsVariable2DTest.testVector1);
//		assertArrayEquals(g.valueRaw(), GibbsVariable2DTest.testAddRowResult1);
//	}

	@Test
	public void testToString() {
		RandomVar2D g = new RandomVar2D(
				GibbsVariate2DTest.testName1,
				GibbsVariate2DTest.testMatrix1);
		assertEquals(g.toString(), GibbsVariate2DTest.testString1);
	}
}
