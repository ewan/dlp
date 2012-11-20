package org.jqgibbs.tests;


import org.jqgibbs.Chain;
import org.jqgibbs.GibbsException;
import org.jqgibbs.Sampler;
import org.jqgibbs.mathstat.rv.RandomVar;
import org.jqgibbs.mathstat.rv.RandomVar1D;
import org.junit.Test;

import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import static org.junit.Assert.*;

public class ChainTest {

	private static String testGibbsVariable1Name = "X";
	private static double[] testGibbsVariable1_1Mat = new double[] {2,4,3,4};
	private static RandomVar testGibbsVariable1_1 = new RandomVar1D(
		testGibbsVariable1Name, testGibbsVariable1_1Mat
	);
	private static double[] testGibbsVariable1_2Mat = new double[] {2,3,3,4};
	private static RandomVar testGibbsVariable1_2 = new RandomVar1D(
			testGibbsVariable1Name, testGibbsVariable1_2Mat
		);

	private static String testChainLink1Str = "1: [\n  " + testGibbsVariable1_1.toString() + "\n]";
	private static String testChainLink2Str = "2: [\n  " + testGibbsVariable1_2.toString() + "\n]";

	private static String testChainStr = testChainLink1Str + "\n" + testChainLink2Str;
	
	@Test
	public void testChain() {
		new Chain();
	}
	
	@Test
	public void testAddLinkLastFirst() throws GibbsException {
		Chain c = new Chain();
		c.addLink(new RandomVar1D((double) 1));
		DenseDoubleMatrix1D d = new DenseDoubleMatrix1D(new double[] {1});
		assertEquals(c.last().get(RandomVar.defaultName()).toColt(), d);
		c.addLink(new RandomVar1D((double) 15, (double) 3));
		DenseDoubleMatrix1D e = new DenseDoubleMatrix1D(new double[] {15, 3});
		assertEquals(c.last().get(RandomVar.defaultName()).toColt(), e);
		assertEquals(c.first().get(RandomVar.defaultName()).toColt(), d);
		c.addLink(new RandomVar1D((double) 2, (double) 5));
		assertEquals(c.first().get(RandomVar.defaultName()).toColt(), d);
	}
	
	@Test
	public void testAddLinkVarargs() throws GibbsException {
		Chain c = new Chain();
		c.addLink(new RandomVar1D("X", (double) 1));
		c.addLink(new RandomVar1D("X", (double) 15, (double) 3), new RandomVar1D("Y", (double) 45));
		assertEquals(((RandomVar1D)c.last().get("Y")).get(0), (double) 45, Sampler.MIN_VALUE);
	}
	
	@Test
	public void testToString() throws GibbsException {
		Chain c = new Chain();
		c.addLink(ChainTest.testGibbsVariable1_1);
		c.addLink(ChainTest.testGibbsVariable1_2);
		assertEquals(c.toString(), ChainTest.testChainStr);
	}
	
	@Test
	public void testIsEmpty() throws GibbsException {
		Chain c = new Chain();
		assertTrue(c.isEmpty());
		c.addLink(ChainTest.testGibbsVariable1_1);
		assertFalse(c.isEmpty());
	}
}
