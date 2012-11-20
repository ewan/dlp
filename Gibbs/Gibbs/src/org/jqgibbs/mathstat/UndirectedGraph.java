package org.jqgibbs.mathstat;

import java.util.LinkedList;

public class UndirectedGraph {
	private int[][] adjacencies;
	
	public UndirectedGraph(Integer1D adjacencies) {
		// 0 = N^2 - N - 2P
		// N = (1 + sqrt(1+8P))/2
		int P = adjacencies.size();
		if (Math.floor((1+Math.sqrt(1+8*P))/2) != (1+Math.sqrt(1+8*P))/2) {
			throw new IllegalArgumentException("Adjacency vector has invalid length");
		}
		int N = (int) (1+Math.sqrt(1+8*P))/2;
		// Set up adjacency matrix
		int[][] a = new int[N][N];
		int pi = 0;
		for (int i=0; i<N; i++) {
			for (int j=(i+1); j<N; j++) {
				a[i][j] = a[j][i] = adjacencies.getValue(pi);
				pi++;
			}
		}
		this.setAdjacencies(a);
	}
	
	protected void setAdjacencies(int[][] a) {
		this.adjacencies = a;
	}
	
	public int[][] getAdjacencies() {
		return this.adjacencies;
	}
	
	public int numNodes() {
		return this.getAdjacencies().length;
	}
	
	public Integer1D components() {
		int[] components = new int[this.numNodes()];
		boolean[] added = new boolean[this.numNodes()];
		boolean[] explored = new boolean[this.numNodes()];
		LinkedList<Integer0D> stack = new LinkedList<Integer0D>();
		int currComponent = -1;
		for (int i=0; i<this.numNodes(); i++) {
			if (!explored[i]) {
				currComponent++;
				stack.push(new Integer0D(i));
				added[i] = true;
				while (!stack.isEmpty()) {
					Integer0D node = stack.pop();
					if (!explored[node.value()]) {
						for (int j=0; j<this.numNodes(); j++) {
							if (this.getAdjacencies()[node.value()][j] != 0) {
								if (!added[j]) {
									stack.add(new Integer0D(j));
									added[j] = true;
								}
							}
						}
						components[node.value()] = currComponent;
						explored[node.value()] = true;
					}
				}
			}
		}
		return new Integer1D(components);
	}
}
