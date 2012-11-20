package org.jqgibbs;

import java.io.IOException;
import java.util.logging.Logger;

import jline.*;

public class UISampler implements Runnable {
	private Chain c;
	private Thread ts;
	private Logger logger;
	
	public UISampler(Thread ts, SamplerData d, Chain c) {
		this.c = c;
		this.ts = ts;
		this.logger = Logger.getLogger("org.jqgibbs.UISampler");		
	}
	
	public void run() {
		Thread.currentThread().setPriority(Thread.MAX_PRIORITY);
		ConsoleReader cr;
		try {
			cr = new ConsoleReader();
			String line = null;
			while (!Thread.interrupted()) {
				line = cr.readLine("> ");
				if (line == null) {
					this.ts.interrupt();
					Thread.currentThread().interrupt();
				} else if (line.equals("print chain")) {
					System.out.println("Right away, sir...");
					System.out.println(this.c);
				} else if (line.equals("print chain length")) {
					System.out.println("Right away, sir...");
					System.out.println(this.c.size());
				}
			}			
		} catch (IOException e) {
			this.logger.warning("Can't read user input: " + e.getMessage());
		}

	}
}
