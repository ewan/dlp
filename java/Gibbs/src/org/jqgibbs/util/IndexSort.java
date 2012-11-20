package org.jqgibbs.util;

// Not my code! Conveniently taken from a forum post on stackoverflow

public class IndexSort {
	
	public static void sort(int[] main, int[] index) {
	    sort(main, index, 0, index.length - 1);
	}

	// quicksort a[left] to a[right]
	public static void sort(int[] a, int[] index, int left, int right) {
	    if (right <= left) return;
	    int i = partition(a, index, left, right);
	    sort(a, index, left, i-1);
	    sort(a, index, i+1, right);
	}

	// partition a[left] to a[right], assumes left < right
	private static int partition(int[] a, int[] index, 
	int left, int right) {
	    int i = left - 1;
	    int j = right;
	    while (true) {
	        while (less(a[++i], a[right]))      // find item on left to swap
	            ;                               // a[right] acts as sentinel
	        while (less(a[right], a[--j]))      // find item on right to swap
	            if (j == left) break;           // don't go out-of-bounds
	        if (i >= j) break;                  // check if pointers cross
	        exch(a, index, i, j);               // swap two elements into place
	    }
	    exch(a, index, i, right);               // swap with partition element
	    return i;
	}

	// is x < y ?
	private static boolean less(int x, int y) {
	    return (x < y);
	}

	// exchange a[i] and a[j]
	private static void exch(int[] a, int[] index, int i, int j) {
	    int swap = a[i];
	    a[i] = a[j];
	    a[j] = swap;
	    int b = index[i];
	    index[i] = index[j];
	    index[j] = b;
	}
	
}
