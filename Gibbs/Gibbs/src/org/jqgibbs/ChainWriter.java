package org.jqgibbs;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

public class ChainWriter {
	String filename;
	Model model;
	private FileOutputStream fos;
	
	public ChainWriter(String filename, Model model) {
		this.filename = filename;
		this.model = model;
		try {
			this.setFos(new FileOutputStream(new File(this.filename)));
		} catch (FileNotFoundException e) {
			// FIXME
			e.printStackTrace();
		} 
	}
	
	public void write(Chain c) {
		try {
			for (ChainLink cl : c) {
				this.getFos().write((cl.toString() + '\n').getBytes());
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace(); // FIXME
		}	
	}

	private void setFos(FileOutputStream fos) { //FIXME
		this.fos = fos;
	}

	private FileOutputStream getFos() {
		return fos;
	}
}
