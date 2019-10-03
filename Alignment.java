// Calvin Vuong ccv7
// A class to perform local and global alignments

import java.util.List;
import java.util.ArrayList;
import java.io.File;
import java.util.Scanner;

public class Alignment {

    // Implements the Smith-Waterman algorithm to perform local alignment
    // Takes as input two Strings s1, s2 representing ATCG sequences and the scores for match, mismatch, and indel
    public static int localAlignment(String s1, String s2, int match, int mismatch, int indel) {
	// length of input strings
	int m = s1.length();
	int n = s2.length();

	int[][] scores = new int[m+1][n+1];    // score solution matrix
	int[][] backtrack = new int[m+1][n+1]; // backtrack matrix
	// uses 3 bit binary (represented as an int) to encode directions
	// up, left, diagonal in that order (e.g. 101 corresponds to diagonal and up)
	
	// initialize score matrix (same as Needleman-Wunsch)
	// initialize score matrix
	for ( int i = 0; i <= m; i++ ) // initialize F(i, 0)
	    scores[i][0] = i * indel;
	for ( int j = 0; j <= n; j++ ) // initialize F(0, j)
	    scores[0][j] = j * indel;

	// initialize backtrack matrix (different from Needleman-Wunsch)
	backtrack[0][0] = 0;
	for ( int i = 1; i <= m; i++ ) // initialize F(i, 0)
	    backtrack[i][0] = 0;
	for ( int j = 1; j <= n; j++ ) // initialize F(0, j)
	    backtrack[0][j] = 0;

	// fill in subproblem matrices
	for ( int i = 1; i <= m; i++ ) {
	    for ( int j = 1; j <= n; j++ ) {
		int sub = scores[i-1][j-1] + sub(s1.charAt(i-1), s2.charAt(j-1), match, mismatch); // opt. subproblem score for substitution
		int gS2 = scores[i-1][j] + indel; // opt. subproblem score for gap in s2
		int gS1 = scores[i][j-1] + indel; // opt. subproblem score for gap in s1
		
		// fill in solution
		int best = max(sub, gS2, gS1);
		// "restart" alignment if best of subproblems gives negative score 
		if ( best < 0 ) 
		    scores[i][j] = 0;
		else 
		    scores[i][j] = best;

		// fill in backtrack data
		if ( scores[i][j] == sub ) // subs. gave you best opt. subsoln.
		    backtrack[i][j] += 1; // point diagonal
		if ( scores[i][j] == gS2 ) // inserting gap in s2 gave best opt. subsoln.
		    backtrack[i][j] += 4; // point up
		if ( scores[i][j] == gS1 ) // inserting gap in s1 gave best opt. subsoln.
		    backtrack[i][j] += 2; // point left
	    }
	}

	// find optimal score in matrix
	int maxScore = scores[0][0];
	// the indices (i, j) of the optimal solutions (there can be more than one)
	// List<int[]> optIndices = new ArrayList<int[]>();
	for ( int i = 0; i <= m; i++ ) {
	    for ( int j = 0; j <= n; j++ ) {
		if ( scores[i][j] > maxScore ) 
		    maxScore = scores[i][j];
	    }
	}

	// find the indices where the max scores are
	List<int[]> maxScoreIndices = new ArrayList<int[]>();
	for ( int i = 0; i <= m; i++ ) {
	    for ( int j = 0; j <= n; j++ ) {
		if ( scores[i][j] == maxScore )
		    maxScoreIndices.add(new int[] {i, j});
	    }
	}
	
	// count the number of optimal alignments
	int numOptimalAlignments = 0;
	for ( int k = 0; k < maxScoreIndices.size(); k++ ) {
	    int i = maxScoreIndices.get(k)[0];
	    int j = maxScoreIndices.get(k)[1];
	    numOptimalAlignments += countAlignments(scores, backtrack, i, j, true);
	}
	
	printMatrix(scores);
	printMatrix(backtrack);
	System.out.println("Optimal score: " + maxScore);
	System.out.println("Number of optimal alignments: " + numOptimalAlignments);

	// print all the optimal alignments
	for ( int k = 0; k < maxScoreIndices.size(); k++ ) {
	    int i = maxScoreIndices.get(k)[0];
	    int j = maxScoreIndices.get(k)[1];
	    printAlignments(s1, s2, scores, backtrack, i, j, true);
	}

	return maxScore;
    }

    
    // Implements the Needleman-Wunsch algorithm to perform global alignment
    // Takes as input two Strings s1, s2 representing ATCG sequences and the scores for match, mismatch, and indel
    // Returns the optimal alignment score
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
	System.out.println("Number of optimal alignments: " + countAlignments(scores, backtrack, m, n, false));
	printAlignments(s1, s2, scores, backtrack, m, n, false);
	
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

    public static void printAlignments(String s1, String s2, int[][] scores, int[][] backtrack, int i, int j, boolean local) {
	// all alignments, stored in a list of String arrays
	List<String[]> alignments = getAlignments(s1, s2, scores, backtrack, i, j, local);
	for ( int k = 0; k < alignments.size(); k++ ) {
	    System.out.println(alignments.get(k)[0]);
	    System.out.println(alignments.get(k)[1]);
	    System.out.println();
	}
	System.out.println();
    }

    // Takes the sequence Strings and matrix representing backtracking data
    // Also takes a starting index to backtrack from
    // Also takes a flag that indicates whether or not it is getting a global or local alignment
    // Prints out all the paths to the optimal solution
    // Recursively defined
    public static List<String[]> getAlignments(String s1, String s2, int[][] scores, int[][] backtrack, int i, int j, boolean local) {
	if ( i == 0 & j == 0 ) { // base case
	    List<String[]> newList = new ArrayList<String[]>();
	    newList.add(new String[] {"", ""});
	    return newList;
	}

	List<String[]> newAlignments = new ArrayList<String[]>();
	// if (i, j) points diagonal
	if ( backtrack[i][j] %  2 == 1 ) {
	    // get all the optimal alignment sequences so far to index (i-1, j-1)
	    List<String[]> setOfAlignments = getAlignments(s1, s2, scores, backtrack, i-1, j-1, local);
	    for ( int k = 0; k < setOfAlignments.size(); k++ ) { // for all alignments in setOfAlignments
		String[] alignment = setOfAlignments.get(k);
		alignment[0] = alignment[0] + s1.charAt(i-1);
		alignment[1] = alignment[1] + s2.charAt(j-1);
		newAlignments.add(alignment);
	    }
	}

	// if (i, j) points left
	if ( backtrack[i][j] >= 2 && backtrack[i][j] / 2 != 2 ) {
	    // get all the optimal alignment sequences so far to index (i, j-1)
	    List<String[]> setOfAlignments = getAlignments(s1, s2, scores, backtrack, i, j-1, local);
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
	    List<String[]> setOfAlignments = getAlignments(s1, s2, scores, backtrack, i-1, j, local);
	    for ( int k = 0; k < setOfAlignments.size(); k++ ) {
		String[] alignment = setOfAlignments.get(k);
		alignment[0] = alignment[0] + s1.charAt(i-1);
		alignment[1] = alignment[1] + "-";
		newAlignments.add(alignment);
	    }
	}

	// if (i, j) can also be a new alignment
	if ( local && scores[i][j] == 0 ) {
	    newAlignments.add(new String[] {"", ""});
	} 

	return newAlignments;
    }

    // Counts the number of optimal alignments
    // Also takes a flag that indicates whether or not it is getting a global or local alignment
    // recursively defined
    public static int countAlignments(int[][] scores, int[][] backtrack, int i, int j, boolean local) {
	// base case
	if ( backtrack[i][j] == 0 ) // end of sequence
	    return 1;
	
	int numAlignments = 0;
	// if (i, j) points diagonal
	if ( backtrack[i][j] % 2 == 1 )
	    // count alignment sequences in the diagonal direction
	    numAlignments += countAlignments(scores, backtrack, i-1, j-1, local);
	// if (i, j) points left
	if ( backtrack[i][j] >= 2 && backtrack[i][j] / 2 != 2 )
	    // count alignment sequences in the left directio
	    numAlignments += countAlignments(scores, backtrack, i, j-1, local);
	// if (i, j) points up
	if ( backtrack[i][j] >= 4 )
	    // count alignment sequences in the up direction
	    numAlignments += countAlignments(scores, backtrack, i-1, j, local);

	// if (i, j) can also be the start of a new sequence
	if ( local && scores[i][j] == 0 )
	    numAlignments += 1;

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


    // Read commands from an input file
    public static void main(String[] args) throws Exception {
	boolean local = false; // true if input specifies local alignment
	int match = 0;
	int mismatch = 0;
	int indel = 0;
	String s1 = "";
	String s2 = "";

	if ( args.length > 0 ) {
	    File commands = new File(args[0]);
	    Scanner scan = new Scanner(commands);
	    int lineNum = 1;
	    // read commands from file line by line
	    while ( scan.hasNextLine() ) {
		String line = scan.nextLine();
		if ( lineNum == 1 && line.toLowerCase().startsWith("g") )
		    local = false;
		else if ( lineNum == 1 && line.toLowerCase().startsWith("l") )
		    local = true;

		else if ( lineNum == 2 ) {
		    String[] penalties = line.split(", | |,");
		    match = Integer.parseInt(penalties[0]);
		    mismatch = Integer.parseInt(penalties[1]);
		    indel = Integer.parseInt(penalties[2]);
		}

		else if ( lineNum == 3 ) 
		    s1 = line;

		else if ( lineNum == 4 )
		    s2 = line;

		lineNum += 1;
	    }

	    // perform alignment based on input file params
	    if ( local ) 
		localAlignment(s1, s2, match, mismatch, indel);
	    else
		globalAlignment(s1, s2, match, mismatch, indel);
	}
	
	else {
	
	    localAlignment("ATTCG", "ATCG", 1, -1, -1);
	    
	/*
	//globalAlignment("ATTTGG", "ATCG", 1, -1, -1);
	globalAlignment("GATTACA", "GCATGCT", 1, -1, -1);
	//localAlignment("ATTTGG", "ATCG", 1, -1, -1);
	localAlignment("GATTACA", "GCATGCT", 1, -1, -1);
	localAlignment("ACTAGG", "AGTAGG", 1, -1, -1);
	localAlignment("GGCACGTTCACC", "TACAGC", 1, -1, -1);
	*/
	}
	
    }

}

