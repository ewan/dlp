package org.jqgibbs.tests;

import static org.junit.Assert.*;

import org.jqgibbs.Sampler;
import org.jqgibbs.mathstat.rv.RandomVar1D;
import org.junit.Test;

import cern.colt.matrix.impl.DenseDoubleMatrix1D;

public class GibbsVariate1DTest {

	private static double[] testMatrix1 = new double[] {1,2,3,4};
	private static double[] testMatrix2 = new double[] {1,1,3,4};
	private static String testName1 = "V";
	private static String testName2 = "W";
	private static String testString1 = "V: {1, 2, 3, 4}";
	
	@Test
	public void testValueColt() {
		RandomVar1D g = new RandomVar1D(GibbsVariate1DTest.testMatrix1);
		DenseDoubleMatrix1D d = new DenseDoubleMatrix1D(GibbsVariate1DTest.testMatrix1);
		assertEquals(g.toColt(), d);
	}

	@Test
	public void testName() {
		RandomVar1D g = new RandomVar1D(
				GibbsVariate1DTest.testName1,
				GibbsVariate1DTest.testMatrix1);
		assertEquals(g.getName(), GibbsVariate1DTest.testName1);
	}
	
	@Test
	public void testHashCode() {
		RandomVar1D g1 = new RandomVar1D(
				GibbsVariate1DTest.testName1,
				GibbsVariate1DTest.testMatrix1);
		RandomVar1D g2 = new RandomVar1D(
				GibbsVariate1DTest.testName1,
				GibbsVariate1DTest.testMatrix1);
		assertEquals(g1.hashCode(), g2.hashCode());
	}

	@Test
	public void testEquals() {
		RandomVar1D g1 = new RandomVar1D(
				GibbsVariate1DTest.testName1,
				GibbsVariate1DTest.testMatrix1);
		RandomVar1D g2 = new RandomVar1D(
				GibbsVariate1DTest.testName1,
				GibbsVariate1DTest.testMatrix1);
		RandomVar1D g3 = new RandomVar1D(
				GibbsVariate1DTest.testName2,
				GibbsVariate1DTest.testMatrix1);
		RandomVar1D g4 = new RandomVar1D(
				GibbsVariate1DTest.testName1,
				GibbsVariate1DTest.testMatrix2);
		assertEquals(g1, g2);
		assertFalse(g1.equals(g3));
		assertFalse(g1.equals(g4));
	}

	@Test
	public void testGet() {
		RandomVar1D g = new RandomVar1D(GibbsVariate1DTest.testMatrix1);
		for (int i=0; i<GibbsVariate1DTest.testMatrix1.length; i++) {
			assertEquals(g.get(i), GibbsVariate1DTest.testMatrix1[i], Sampler.MIN_VALUE);
		}
	}

	@Test
	public void testValueRaw() {
		RandomVar1D g = new RandomVar1D(GibbsVariate1DTest.testMatrix1);
		assertEquals(g.toDouble(), GibbsVariate1DTest.testMatrix1);
	}
	
	@Test
	public void testToString() {
		RandomVar1D g = new RandomVar1D(GibbsVariate1DTest.testName1, GibbsVariate1DTest.testMatrix1);
		assertEquals(g.toString(), GibbsVariate1DTest.testString1);
	}

}
