package org.jqgibbs.tests;

import static org.junit.Assert.*;

import org.jqgibbs.ChainLink;
import org.jqgibbs.GibbsException;
import org.jqgibbs.mathstat.rv.RandomVar;
import org.jqgibbs.mathstat.rv.RandomVar1D;
import org.junit.Test;

public class ChainLinkTest {

	private static String testGibbsVariable1Name = "X";
	private static double[] testGibbsVariable1Mat = new double[] {1,2,3,4};
	private static RandomVar testGibbsVariable1 = new RandomVar1D(
		testGibbsVariable1Name, testGibbsVariable1Mat
	);
	
	private static double[] testGibbsVariable2Mat = new double[] {2,4,3,4};
	private static RandomVar testGibbsVariable2 = new RandomVar1D(
		testGibbsVariable1Name, testGibbsVariable2Mat
	);
	//private static String testGibbsVariable2Str = "X: [2,4,3,4]";
	
	private static String testGibbsVariable3Name = "Y";
	private static RandomVar testGibbsVariable3 = new RandomVar1D(
		testGibbsVariable3Name, testGibbsVariable2Mat
	);
	
	private static String testChainLink13Str = "[\n  " + testGibbsVariable1.toString() + "\n  " + testGibbsVariable3.toString() + "\n]";
	
	@Test
	public void testChainLink() throws GibbsException {
		new ChainLink(ChainLinkTest.testGibbsVariable1);
	}
	
	@Test
	public void testChainLinkDuplicateName() {
		try {
			new ChainLink(
				ChainLinkTest.testGibbsVariable1,
				ChainLinkTest.testGibbsVariable2
			);
		} catch (GibbsException e) {
			assertTrue(true);
		}
	}
	
	@Test
	public void testChainLinkDistinctNames() {
		try {
			new ChainLink(
				ChainLinkTest.testGibbsVariable1,
				ChainLinkTest.testGibbsVariable3
			);
		} catch (GibbsException e) {
			assertTrue(true);
		}
	}

	@Test
	public void testGet() throws GibbsException {
		ChainLink cl = new ChainLink(
				ChainLinkTest.testGibbsVariable1,
				ChainLinkTest.testGibbsVariable3
			);
		assertEquals(cl.get(ChainLinkTest.testGibbsVariable1Name), ChainLinkTest.testGibbsVariable1);
		assertEquals(cl.get(ChainLinkTest.testGibbsVariable3Name), ChainLinkTest.testGibbsVariable3);
	}

	@Test
	public void testHashCode() throws GibbsException {
		ChainLink cl1 = new ChainLink(
				ChainLinkTest.testGibbsVariable1,
				ChainLinkTest.testGibbsVariable3
			);
		ChainLink cl2 = new ChainLink(
				ChainLinkTest.testGibbsVariable3,
				ChainLinkTest.testGibbsVariable1
			);
		assertEquals(cl1.hashCode(), cl2.hashCode());
	}

	@Test
	public void testEquals() throws GibbsException {
		ChainLink cl1 = new ChainLink(
				ChainLinkTest.testGibbsVariable1,
				ChainLinkTest.testGibbsVariable3
			);
		ChainLink cl2 = new ChainLink(
				ChainLinkTest.testGibbsVariable3,
				ChainLinkTest.testGibbsVariable1
			);
		ChainLink cl3 = new ChainLink(
				ChainLinkTest.testGibbsVariable1
			);
		ChainLink cl4 = new ChainLink(
				ChainLinkTest.testGibbsVariable2,
				ChainLinkTest.testGibbsVariable3
			);
		assertEquals(cl1, cl2);
		assertFalse(cl1.equals(cl3));
		assertFalse(cl1.equals(cl4));
	}

	@Test
	public void testToString() throws GibbsException {
		ChainLink cl = new ChainLink(
				ChainLinkTest.testGibbsVariable1,
				ChainLinkTest.testGibbsVariable3
			);
		assertEquals(cl.toString(), ChainLinkTest.testChainLink13Str);
	}

}
