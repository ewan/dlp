package org.jqgibbs;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

public class ChainWriter {
	String filename;
	Model model;
	private FileOutputStream fos;

	public ChainWriter(String filename, Model model) throws IOException {
		this.filename = filename;
		this.model = model;
		this.fos = new FileOutputStream(new File(this.filename));
	}

	public void write(Chain c) throws IOException {
		for (ChainLink cl : c) {
			this.fos.write((cl.toString() + '\n').getBytes());
		}
	}
}
