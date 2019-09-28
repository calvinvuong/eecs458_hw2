// Calvin Vuong ccv7
// A class to perform local and global alignments

import java.util.List;
import java.util.ArrayList;

public class Alignment {
    // Implements the Needleman-Wunsch algorithm to perform global alignment
    // Takes as input to Strings s1, s2 representing ATCG sequences and the scores for match, mismatch, and indel
    public static int globalAlignment(String s1, String s2, int match, int mismatch, int indel) {
	// length of input strings
	int m = s1.length();
	int n = s2.length();
	
	int[][] scores = new int[m+1][n+1]; 	//  score solution matrix
	int[][] backtrack = new int[m+1][n+1]; 	// backtrack matrix
	// uses 3 bit binary (represented as an int) to encode directions
	// up, left, diagonal in that order (e.g. 101 corresponds to diagonal and up)
	
	// initialize score matrix
	for ( int i = 0; i <= m; i++ ) // initialize F(i, 0)
	    scores[i][0] = i * indel;
	for ( int j = 0; j <= n; j++ ) // initialize F(0, j)
	    scores[0][j] = j * indel;

	// initialize backtrack matrix
	backtrack[0][0] = 0;
	for ( int i = 1; i <= m; i++ ) // initialize F(i, 0)
	    backtrack[i][0] = 4;
	for ( int j = 1; j <= n; j++ ) // initialize F(0, j)
	    backtrack[0][j] = 2;
	
	
	// fill in subproblem matrices
	for ( int i = 1; i <= m; i++ ) {
	    for ( int j = 1; j <= n; j++ ) {
		int sub = scores[i-1][j-1] + sub(s1.charAt(i-1), s2.charAt(j-1), match, mismatch); // opt. subproblem score for substitution
		int gS2 = scores[i-1][j] + indel; // opt. subproblem score for gap in s2
		int gS1 = scores[i][j-1] + indel; // opt. subproblem score for gap in s1

		// fill in solution
		int best = max(sub, gS2, gS1);
		scores[i][j] = best;

		// fill in backtrack data
		if ( best == sub ) // substitution gave you best optimal subsolution
		    backtrack[i][j] += 1; // point diagonal
		if ( best == gS2 ) // inserting gap in s2 gave best optimal subsolution
		    backtrack[i][j] += 4; // point up
		if ( best == gS1 ) // inserting gap in s1 gave best optimal subsolution
		    backtrack[i][j] += 2; // point left
	    }

	}

	
	printMatrix(scores);
	printMatrix(backtrack);
	System.out.println("Optimal score: " + scores[m][n]);
	System.out.println("Number of optimal alignments: " + countAlignments(backtrack, m, n));
	printAlignments(s1, s2, backtrack, m, n);
	
	return scores[m][n];
    }


    // Takes four inputs: chars x,y and ints match, mismatch
    // Returns match if chars match, mismatch if chars mismatch
    public static int sub(char x, char y, int match, int mismatch) {
	if ( x == y )
	    return match;
	else
	    return mismatch;
    }

    public static void printAlignments(String s1, String s2, int[][] backtrack, int i, int j) {
	// all alignments, stored in a list of String arrays
	List<String[]> alignments = getAlignments(s1, s2, backtrack, i, j);
	for ( int k = 0; k < alignments.size(); k++ ) {
	    System.out.println(alignments.get(k)[0]);
	    System.out.println(alignments.get(k)[1]);
	    System.out.println();
	}
	System.out.println();
    }

    // Takes the sequence Strings and matrix representing backtracking data
    // Also takes a starting index to backtrack from
    // Prints out all the paths to the optimal solution
    // Recursively defined
    public static List<String[]> getAlignments(String s1, String s2, int[][] backtrack, int i, int j) {
	if ( i == 0 & j == 0 ) {
	    List<String[]> newList = new ArrayList<String[]>();
	    newList.add(new String[] {"", ""});
	    return newList;
	}

	List<String[]> newAlignments = new ArrayList<String[]>();
	// if (i, j) points diagonal
	if ( backtrack[i][j] %  2 == 1 ) {
	    // get all the optimal alignment sequences so far to index (i-1, j-1)
	    List<String[]> setOfAlignments = getAlignments(s1, s2, backtrack, i-1, j-1);
	    for ( int k = 0; k < setOfAlignments.size(); k++ ) { // for all alignments in alignPrefix
		String[] alignment = setOfAlignments.get(k);
		alignment[0] = alignment[0] + s1.charAt(i-1);
		alignment[1] = alignment[1] + s2.charAt(j-1);
		newAlignments.add(alignment);
	    }
	}

	// if (i, j) points left
	if ( backtrack[i][j] >= 2 && backtrack[i][j] / 2 != 2 ) {
	    // get all the optimal alignment sequences so far to index (i, j-1)
	    List<String[]> setOfAlignments = getAlignments(s1, s2, backtrack, i, j-1);
	    for ( int k = 0; k < setOfAlignments.size(); k++ ) {
		String[] alignment = setOfAlignments.get(k);
		alignment[0] = alignment[0] + "-";
		alignment[1] = alignment[1] + s2.charAt(j-1);
		newAlignments.add(alignment);
	    }
	}

	// if ( i, j) points up
	if ( backtrack[i][j] >= 4 ) {
	    // get all the optimal alignment sequences so far to index (i-1, j)
	    List<String[]> setOfAlignments = getAlignments(s1, s2, backtrack, i-1, j);
	    for ( int k = 0; k < setOfAlignments.size(); k++ ) {
		String[] alignment = setOfAlignments.get(k);
		alignment[0] = alignment[0] + s1.charAt(i-1);
		alignment[1] = alignment[1] + "-";
		newAlignments.add(alignment);
	    }
	}

	return newAlignments;
    }

    // counts the number of optimal alignments
    // recursively defined
    public static int countAlignments(int[][] backtrack, int i, int j) {
	// base case
	if ( backtrack[i][j] == 0 ) // end of sequence
	    return 1;
	
	int numAlignments = 0;
	// if (i, j) points diagonal
	if ( backtrack[i][j] % 2 == 1 )
	    // count alignment sequences in the diagonal direction
	    numAlignments += countAlignments(backtrack, i-1, j-1);
	// if (i, j) points left
	if ( backtrack[i][j] >= 2 && backtrack[i][j] / 2 != 2 )
	    // count alignment sequences in the left directio
	    numAlignments += countAlignments(backtrack, i, j-1);
	// if (i, j) points up
	if ( backtrack[i][j] >= 4 )
	    // count alignment sequences in the up direction
	    numAlignments += countAlignments(backtrack, i-1, j);

	return numAlignments;
    }
    
    // Prints the input matrix
    // For debugging purposes
    public static void printMatrix(int[][] matrix) {
	for ( int i = 0; i < matrix.length; i++ ) {
	    for ( int j = 0; j < matrix[i].length; j++ ) 
		System.out.print(matrix[i][j] + "\t");
	    System.out.println();
	}
	System.out.println();
    }
    
    // Returns the maximum of the three input ints: x, y, z
    private static int max(int x, int y , int z) {
	if ( x >= y && x >= z )
	    return x;
	if ( y >= x && y >= z )
	    return y;
	else
	    return z;
    }

   
    public static void main(String[] args) {
	globalAlignment("ATTTGG", "ATCG", 1, -1, -1);
	globalAlignment("GATTACA", "GCATGCT", 1, -1, -1);
	/*
	System.out.println( max(5, 5, 5) );
	System.out.println( max(3, 5, 5) );
	System.out.println( max(5, 3, 5) );
	System.out.println( max(5, 3, 3) );
	System.out.println( max(4, 6, 8) );
	System.out.println( max(7, 7, 3) );
	System.out.println( max(7, 10, 5) );
	*/
    }

}

