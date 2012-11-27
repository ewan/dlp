package org.jqgibbs;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Observable;

import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.dom.DOMResult;

import org.jqgibbs.mathstat.Double0D;
import org.jqgibbs.mathstat.Double1D;
import org.jqgibbs.mathstat.Double2D;
import org.jqgibbs.mathstat.Integer0D;
import org.jqgibbs.mathstat.Integer1D;
import org.jqgibbs.mathstat.Numeric;
import org.jqgibbs.mathstat.RandomVar;
import org.jqgibbs.mathstat.probdist.ProbDistParmException;
import org.jqgibbs.util.DOMValidateDTD;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class Chain extends Observable implements List<ChainLink> {

	private boolean flattened;
	private Map<String,Integer1D> flatOm;
	private Double2D flatDouble2D;
	private static String localDTD = Chain.class.getResource("chain.dtd").getFile();

	private static enum NodeName {
		matrix, model, parametercolumns, name, hypers, namedrow, namedmatrix, dims, namedvalue
	}

	private static Numeric nodeToValue(Node n) {
		n.normalize();
		Node v = n.getFirstChild();
		Node type = n.getAttributes().getNamedItem("type");
		if (type != null && type.getNodeValue().equals("int")) {
			return new Integer0D(Integer.parseInt(v.getNodeValue()));
		} else {
			return new Double0D(Double.parseDouble(v.getNodeValue()));
		}
	}

	private static Double1D nodeToRow(Node n) throws GibbsException {
		n.normalize();
		NodeList nl = n.getChildNodes();
		Node v;
		Double1D d1D = new Double1D();
		for (int i = 0; i < nl.getLength(); i++) {
			v = nl.item(i);
			if (v.getNodeType() != Node.ELEMENT_NODE) {
				continue;
			}
			Numeric d = Chain.nodeToValue(v);
			if (d instanceof Double0D) {
				d1D.add((Double0D) d);
			} else {
				throw new GibbsException(
						"Error in chain file: double row cannot contain integer value");
			}
		}
		return d1D;
	}

	private static Double2D nodeToMatrix(Node n) throws GibbsException {
		NodeList nl = n.getChildNodes();
		Node r;
		Double2D d2D = new Double2D();
		for (int i = 0; i < nl.getLength(); i++) {
			r = nl.item(i);
			if (r.getNodeType() != Node.ELEMENT_NODE) {
				continue;
			}
			d2D.add(Chain.nodeToRow(r));
		}
		return d2D;
	}

	private static String nodeToModelName(Node n) {
		n.normalize();
		Node m = n.getFirstChild();
		return m.getNodeValue();
	}

	private static Model nodeToModel(Node n,
			Map<String, Class<? extends Model>> modelTable)
			throws GibbsException {
		NodeList nl = n.getChildNodes();
		Node m;
		String modelName = null; // FIXME - also below!
		Map<String, Numeric> hypers = null;
		Model model = null;
		int dims = 0; // FIXME -especially here..
		for (int i = 0; i < nl.getLength(); i++) {
			m = nl.item(i);
			if (m.getNodeType() != Node.ELEMENT_NODE) {
				continue;
			}
			try {
				switch (Chain.NodeName.valueOf(m.getNodeName())) {
				case name:
					modelName = Chain.nodeToModelName(m);
					break;
				case hypers:
					hypers = Chain.nodeToMap(m);
					break;
				case dims:
					NodeList mnl = m.getChildNodes();
					Node mm = null;
					for (int j=0; j<mnl.getLength(); j++) {
						mm = mnl.item(j);
						if (mm.getNodeType() != Node.ELEMENT_NODE) {
							continue;
						} else {
							break;
						}
					}
					Numeric dn = Chain.nodeToValue(mm);
					if (dn instanceof Double0D) {
						dims = (int) ((Double0D) dn).value();
					} else if (dn instanceof Integer0D) {
						dims = (int) ((Integer0D) dn).value();
					} else {
						assert false;
					}
					break;
				}
			} catch (IllegalArgumentException e) {
				assert false; // Should never get here, since we validated
			}
		}
		try {
			model = modelTable.get(modelName).getConstructor(Map.class,
					int.class).newInstance(hypers, dims);
		} catch (InvocationTargetException e) {
			throw (ProbDistParmException) e.getCause();
		} catch (IllegalAccessException e) {
			throw new GibbsException("Error initializing model " + modelName, e);
		} catch (InstantiationException e) {
			throw new GibbsException("Error initializing model " + modelName, e);
		} catch (NoSuchMethodException e) {
			throw new GibbsException("Error initializing model " + modelName, e);
		}
		return model;
	}

	private static Map<String, Numeric> nodeToMap(Node n)
			throws GibbsException {
		NodeList nl = n.getChildNodes();
		Node c;
		Map<String, Numeric> m = new HashMap<String, Numeric>();
		String s;
		for (int i = 0; i < nl.getLength(); i++) {
			c = nl.item(i);
			if (c.getNodeType() != Node.ELEMENT_NODE) {
				continue;
			}
			s = c.getAttributes().getNamedItem("name").getNodeValue();
			switch (Chain.NodeName.valueOf(c.getNodeName())) {
			case namedrow:
				m.put(s, Chain.nodeToRow(c));
				break;
			case namedmatrix:
				m.put(s, Chain.nodeToMatrix(c));
				break;
			case namedvalue:
				m.put(s, Chain.nodeToValue(c));
				break;
			}
		}
		return m;
	}

	private static Map<String, Integer1D> nodeToMapInteger1D(Node n)
			throws GibbsException {
		NodeList nl = n.getChildNodes();
		Node r;
		Map<String, Integer1D> m = new HashMap<String, Integer1D>();
		String s;
		Integer1D i1D;
		for (int i = 0; i < nl.getLength(); i++) {
			r = nl.item(i);
			if (r.getNodeType() != Node.ELEMENT_NODE) {
				continue;
			}
			s = r.getAttributes().getNamedItem("name").getNodeValue();
			i1D = Chain.nodeToRow(r).toInteger1D();
			m.put(s, i1D);
		}
		return m;
	}

	private List<ChainLink> c;
	private Model model;

	public Chain() {
		this.flattened = false;
		this.setC(new LinkedList<ChainLink>());
	}

	public Chain(String filename, Map<String, Class<? extends Model>> modelTable)
			throws GibbsException {
		// Create empty chain
		super();
		// // Replace the DTD (or add one if there isn't one) with the
		// // DTD for our file format
		// XMLEventFactory ef = XMLEventFactory.newInstance();
		// XMLEvent dtd = ef.createDTD(Chain.localDoctype);
		// XMLInputFactory inFactory = XMLInputFactory.newInstance();
		// inFactory.setProperty("javax.xml.stream.isValidating", "true");
		try {
			// // To replace the DTD, construct an EventReaderDelegate that
			// filters
			// // the output of the (XMLEventReader's) nextEvent() method
			// XMLEventReader reader = inFactory
			// .createXMLEventReader(new StreamSource(filename));
			// reader = new DTDReplacer(reader, dtd);
			// // Jump through some hoops to install the filter
			// StAXSource ss = new StAXSource(reader);
			// TransformerFactory tf = TransformerFactory.newInstance();
			// Transformer t = tf.newTransformer();
			// // Parse the XML
			// DOMResult d = new DOMResult();
			// t.transform(ss, d);
			File xmlFile = new File(filename);
			Document d = DOMValidateDTD.validatedDOM(xmlFile,
					Chain.localDTD);
			Node root = d.getDocumentElement();
			root.normalize();
			// Extract elements
			NodeList nl = root.getChildNodes();
			Node n;
			Double2D matrix = null;
			Map<String, Integer1D> namesToColumns = null;
			for (int i = 0; i < nl.getLength(); i++) {
				n = nl.item(i);
				if (n.getNodeType() != Node.ELEMENT_NODE) {
					continue;
				}
				try {
					switch (Chain.NodeName.valueOf(n.getNodeName())) {
					case matrix:
						matrix = Chain.nodeToMatrix(n);
						break;
					case model:
						this.setModel(Chain.nodeToModel(n, modelTable));
						break;
					case parametercolumns:
						namesToColumns = Chain.nodeToMapInteger1D(n);
						break;
					}
				} catch (IllegalArgumentException e) {
					assert false; // Should never get here, since we validated
				}
			}
			for (Double1D r : matrix) {
				this.getC()
						.add(this.getModel().getChainLink(r, namesToColumns));
			}
		} catch (FileNotFoundException e) {
			throw new GibbsException("Error opening chain file", e);
		} catch (IOException e) {
			throw new GibbsException("Error reading chain file", e);
		} catch (SAXException e) {
			throw new GibbsException("Error reading chain file", e);
		} catch (TransformerConfigurationException e) {
			throw new GibbsException("Error initializing XML reader", e);
		} catch (TransformerException e) {
			throw new GibbsException("Error parsing chain file", e);
		}
	}

	public synchronized void addLink(ChainLink l) {
		this.c.add(l);
		this.flattened = false;
		this.setChanged();
		this.notifyObservers();
	}

	public ChainLink addLink(RandomVar<?>... vs)
			throws DuplicateVariableNameException {
		ChainLink l = new ChainLink(vs);
		this.addLink(l);
		this.flattened = false;
		return l;
	}

	public Iterator<ChainLink> iterator() {
		return this.c.iterator();
	}

	public ChainLink first() {
		return this.c.get(0);
	}

	public ChainLink last() {
		return this.c.get(this.c.size() - 1);
	}

	public ChainLink current() {
		return this.last();
	}

	@Override
	public String toString() {
		String s = "";
		String prefix = "";
		synchronized (this) {
			for (int i = 0; i < this.c.size(); i++) {
//				System.out.println(String.valueOf(i));
				s += prefix + Integer.toString(i + 1) + ": " + this.c.get(i);
				prefix = "\n";
			}
		}
		return s;
	}

	public boolean isEmpty() {
		return this.c.isEmpty();
	}

	public int size() {
		return this.c.size();
	}
	
	public ChainLink get(int i) {
		return this.c.get(i);
	}
	
	private void setC(List<ChainLink> c) {
		this.c = c;
		this.flattened = false;
	}

	protected List<ChainLink> getC() {
		return this.c;
	}

	private void setModel(Model model) {
		this.model = model;
		this.flattened = false;
	}

	public Model getModel() {
		return this.model;
	}
	
	private void flatten() {
		if (this.flattened) {
			return;
		} else {
			// Find longest link
			int max = -1;
			int argmax = -1;
			for (int i=0; i<this.size(); i++) {
				ChainLink cl = this.get(i);
				int length = cl.length1D();
				if (length > max) {
					max = length; 
					argmax = i;
				}
			}
			// Get the order map for the longest link
			this.flatOm = this.get(argmax).getOrderMap();
			// Create flattened structure
			int rows = this.size();
			int cols = max;
			double[][] f2d = new double[rows][cols];
			// Fill in the matrix
			for (int i=0; i<this.size(); i++) {
				double[] row = this.get(i).paddedVec(this.flatOm).value();
				for (int j=0; j<row.length; j++) {
					f2d[i][j] = row[j];
				}
			}
			// Set flattened
			this.flatDouble2D = new Double2D(f2d);			
			this.flattened = true;
		}	
	}
	
	public Map<String,Integer1D> getOrderMap() {
		this.flatten();
		return this.flatOm;
	}
	
	public Double2D getFlatDouble2D() {
		this.flatten();
		return this.flatDouble2D;
	}

	public boolean add(ChainLink arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public void add(int arg0, ChainLink arg1) {
		// TODO Auto-generated method stub
		
	}

	public boolean addAll(Collection<? extends ChainLink> arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public boolean addAll(int arg0, Collection<? extends ChainLink> arg1) {
		// TODO Auto-generated method stub
		return false;
	}

	public void clear() {
		// TODO Auto-generated method stub
		
	}

	public boolean contains(Object arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public boolean containsAll(Collection<?> arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public int indexOf(Object arg0) {
		// TODO Auto-generated method stub
		return 0;
	}

	public int lastIndexOf(Object arg0) {
		// TODO Auto-generated method stub
		return 0;
	}

	public ListIterator<ChainLink> listIterator() {
		// TODO Auto-generated method stub
		return null;
	}

	public ListIterator<ChainLink> listIterator(int arg0) {
		// TODO Auto-generated method stub
		return null;
	}

	public boolean remove(Object arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public ChainLink remove(int arg0) {
		// TODO Auto-generated method stub
		return null;
	}

	public boolean removeAll(Collection<?> arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public boolean retainAll(Collection<?> arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	public ChainLink set(int arg0, ChainLink arg1) {
		// TODO Auto-generated method stub
		return null;
	}

	public List<ChainLink> subList(int arg0, int arg1) {
		// TODO Auto-generated method stub
		return null;
	}

	public Object[] toArray() {
		// TODO Auto-generated method stub
		return null;
	}

	public <T> T[] toArray(T[] arg0) {
		// TODO Auto-generated method stub
		return null;
	}
}
