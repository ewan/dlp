package org.jqgibbs.util;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.xml.sax.ErrorHandler;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

public abstract class DOMValidateDTD {
	public static Document validatedDOM(File xmlFile, String dtdName)
			throws SAXException, IOException, TransformerException {
		ErrorHandler eh = new ErrorHandler() {
			// Ignore the fatal errors
			public void fatalError(SAXParseException e) throws SAXException {
				throw e;
			}

			// Validation errors
			public void error(SAXParseException e) throws SAXParseException {
				throw e;
			}

			// Show warnings
			public void warning(SAXParseException e) throws SAXParseException {
				throw e;
			}
		};
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = null;
		try {
			builder = factory.newDocumentBuilder();
		} catch (ParserConfigurationException e) {
			assert false;
		}
		builder.setErrorHandler(eh);
		// Add DTD
		Document xmlDocumentNoDtd = builder.parse(xmlFile);
		DOMSource sourceNoDtd = new DOMSource(xmlDocumentNoDtd);
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		StreamResult resultWithDtd = new StreamResult(os);
		TransformerFactory tf = TransformerFactory.newInstance();
		Transformer transformer = tf.newTransformer();
		transformer.setOutputProperty(OutputKeys.DOCTYPE_SYSTEM, dtdName);
		transformer.transform(sourceNoDtd, resultWithDtd);
		// Parse new XML
		factory.setValidating(true);
		try {
			builder = factory.newDocumentBuilder();
		} catch (ParserConfigurationException e) {
			assert false;
		}
		builder.setErrorHandler(eh);
		InputStream xmlStreamWithDtd = new ByteArrayInputStream(os.toByteArray());
		Document xmlDocumentWithDtd = builder.parse(xmlStreamWithDtd);
		return xmlDocumentWithDtd;
	}
}