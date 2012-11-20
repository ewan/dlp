package org.jqgibbs;

public class Formatting {
	public static String rep(String s, int n) {
		String rep = "";
		for (int i=0; i<n; i++) {
			rep += s;
		}
		return rep;
	}
	
	public static String catLinesWithSubsetIndented(String[] lines, int n, int first, int last) {
		String s = "";
		String prefix = "";
		for (int i=0; i<lines.length; i++) {
			if (i >= first && i <= last) {
				s += prefix + rep(" ", n) + lines[i];
			} else {
				s += prefix + lines[i];
			}
			prefix = "\n";
		}
		return s;
	}
	
	public static String catLinesWithSubsetDedented(String[] lines, int n, int first, int last) {
		String s = "";
		String prefix = "";
		for (int i=0; i<lines.length; i++) {
			if (i >= first && i <= last && n < lines[i].length()) {
				s += prefix + lines[i].substring(n);
			} else {
				s += prefix + lines[i];
			}
			prefix = "\n";
		}
		return s;
	}
	
	public static String indentAll(String s, int n) {
		String[] lines = s.split("\n");
		return catLinesWithSubsetIndented(lines, n, 0, lines.length);
	}
	
	public static String indentSubsequent(String s, int n) {
		String[] lines = s.split("\n");
		return catLinesWithSubsetIndented(lines, n, 1, lines.length);
	}

	public static String dedent(String s, int n) {
		String[] lines = s.split("\n");
		return catLinesWithSubsetDedented(lines, n, 0, lines.length);
	}
}
