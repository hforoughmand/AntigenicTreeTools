/*************************************************************************
*
*    This source file is part of the software to infer antigenic trees.
*    Copyright (C) 2012  Lars Steinbrueck
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
**************************************************************************/

import java.util.*;
import java.io.*;
import Jama.*;

public class NNLSsolver {
	// native method
	native float[] solveBVLS (float A[], float b[], float bnd[], int n, int m);
	// load the library
	static {
		System.loadLibrary("bvlslib");
	}
	
	public NNLSsolver () {}
	
	// Calculate X matrix ----------------------------------------- //
	
	public int setBranchIndices_updown (node curNode, int index) {
		if (curNode.parent != null) {
			curNode.branchIndex_up = index;
			curNode.branchIndex_down = index * (-1);
			index++;
		}
		for (int i = 0; i < curNode.numChilds; i++) index = this.setBranchIndices_updown(curNode.childArray [i], index);
		return index;
	}
	
	private void setSubtreeLeafs (node curNode) {
		if (curNode.numChilds == 0) {
			curNode.subtreeLeafs = new node [1];
			curNode.subtreeLeafs [0] = curNode;
		}
		else {
			int	numSubtreeLeafs = 0,
				index = 0;
			for (int i = 0; i < curNode.numChilds; i++) {
				this.setSubtreeLeafs(curNode.childArray [i]);
				numSubtreeLeafs += curNode.childArray [i].subtreeLeafs.length;
			}
			curNode.subtreeLeafs = new node [numSubtreeLeafs];
			for (int i = 0; i < curNode.numChilds; i++) {
				for (int j = 0; j < curNode.childArray [i].subtreeLeafs.length; j++) {
					curNode.subtreeLeafs [index+j] = curNode.childArray [i].subtreeLeafs [j];
				}
				index += curNode.childArray [i].subtreeLeafs.length;
			}
		}
	}
	
	private byte[][] calculateXMatrix (node curNode, byte X [][], HashMap<String,Integer> labelMap, HashMap<String,Integer> leafIdentifiers, HashMap<Integer,String> leafIDs) {
		boolean	subtreeLeafsArray [] = null;
		
		if (curNode.numChilds > 0) {
			if (curNode.parent != null) {
				subtreeLeafsArray = new boolean [leafIdentifiers.size()];
				for (int i = 0; i < curNode.subtreeLeafs.length; i++) {
					subtreeLeafsArray [leafIdentifiers.get(curNode.subtreeLeafs[i].identifier)] = true;
				}
				for (int i = 0; i < leafIdentifiers.size(); i++) {
					if (!subtreeLeafsArray [i]) continue;
					for (int j = 0; j < leafIdentifiers.size(); j++) {
						if (subtreeLeafsArray [j] || i == j) continue;
						
						if (labelMap.containsKey (leafIDs.get(i) + "-" + leafIDs.get(j))) X [labelMap.get(leafIDs.get(i) + "-" + leafIDs.get(j))][curNode.branchIndex_up] = 1;
						if (labelMap.containsKey (leafIDs.get(j) + "-" + leafIDs.get(i))) X [labelMap.get(leafIDs.get(j) + "-" + leafIDs.get(i))][X[0].length/2 + (-1) * curNode.branchIndex_down] = 1;
					}
				}
			}
			for (int i = 0; i < curNode.numChilds; i++) X = this.calculateXMatrix(curNode.childArray [i], X, labelMap, leafIdentifiers, leafIDs);
		}
		else {
			for (int i = 0; i < leafIDs.size(); i++) {
				if (labelMap.containsKey (curNode.identifier + "-" + leafIDs.get(i))) X [labelMap.get(curNode.identifier + "-" + leafIDs.get(i))][curNode.branchIndex_up] = 1;
				if (labelMap.containsKey (leafIDs.get(i) + "-" + curNode.identifier)) X [labelMap.get(leafIDs.get(i) + "-" + curNode.identifier)][X[0].length/2 + (-1) * curNode.branchIndex_down] = 1;
			}
		}
		
		return X;
	}
	
	private List<Object> getMainData (node curNode, List<Object>data) {
		// ArrayList; 0 = labels (String[]), 1 = distance matrix (Matrix), 2 = isThreshold (boolean[]), 3 = isSerum (HasHSet<String>), 4 = isAntigen (HashSet<String>)
		int	numBranches = 0,
			numLeaves = 0,
			numUsed = ((String[]) data.get(0)).length;
		byte	X [][] = null;
		String	labels [] = (String []) data.get(0);
		HashMap<String,Integer>	labelMap = new HashMap<String,Integer> (),
								leafIdentifiers = new HashMap<String,Integer> ();
		HashMap<Integer,String>	leafIDs = new HashMap<Integer,String> ();
		List<Object>	resList = new ArrayList<Object> (),
						prunedData = null;
		Matrix	distances = (Matrix) data.get(1);

		numBranches = this.setBranchIndices_updown(curNode, 0);
		
		// set leafs in subtrees
		this.setSubtreeLeafs (curNode);
		numLeaves = curNode.subtreeLeafs.length;
		for (int i = 0; i < numLeaves; i++) {
			leafIdentifiers.put(curNode.subtreeLeafs [i].identifier, i);
			leafIDs.put(i, curNode.subtreeLeafs [i].identifier);
		}
		
		// initialize X matrix
		X = new byte [numUsed][2 * numBranches];
		
		// fill labelMap
		for (int i = 0; i < numUsed; i++) labelMap.put(labels[i], i);
		
		// get matrix
		X = this.calculateXMatrix(curNode, X, labelMap, leafIdentifiers, leafIDs);
		
		// separate data
		String	labels_known [] = null,
				labels_threshold [] = null;
		int	num_known = 0,
			num_threshold = 0,
			index_known = 0,
			index_threshold = 0;
		byte	X_known [][] = null,
				X_threshold [][] = null;
		boolean	isThreshold [] = (boolean []) data.get(2);
		Matrix	distances_known = null,
				distances_threshold = null;
		
		for (int i = 0; i < numUsed; i++) {
			if (isThreshold [i]) num_threshold++;
		}
		
		// set known
		if (num_threshold == 0) {
			num_known = numUsed;
			labels_known = labels;
			X_known = X;
			distances_known = distances;
		}
		else {
			num_known = numUsed - num_threshold; 
			labels_known = new String [num_known];
			labels_threshold = new String [num_threshold];
			X_known = new byte [num_known][];
			X_threshold = new byte [num_threshold][];
			distances_known = new Matrix (num_known, 1);
			distances_threshold = new Matrix (num_threshold, 1);
			for (int i = 0; i < numUsed; i++) {
				if (!isThreshold [i]) {
					X_known [index_known] = X [i];
					distances_known.set(index_known, 0, distances.get(i, 0));
					labels_known [index_known++] = labels [i];
				}
				else {
					X_threshold [index_threshold] = X [i];
					distances_threshold.set(index_threshold, 0, distances.get(i, 0));
					labels_threshold [index_threshold++] = labels [i];
				}
			}
		}
		
		// prune known data
		prunedData = this.pruneXMatrix (X_known);
		
		// collect results
		resList.add(labels_known);
		resList.add(distances_known);
		resList.addAll(prunedData);
		resList.add(labels_threshold);
		resList.add(distances_threshold);
		if (num_threshold == 0) resList.add(null);
		else {
			Matrix	X_threshold_matrix = new Matrix (num_threshold, X [0].length);
			for (int i = 0; i < X_threshold.length; i++)
				for (int j = 0; j < X_threshold [i].length; j++)
					X_threshold_matrix.set(i, j, X_threshold [i][j]);
			resList.add(X_threshold_matrix);
		}
		
		return resList;
		// labels (known): String []
		// distances (known): Matrix n x 1
		// X (known, pruned): Matrix n x m
		// index mapping: int []
		// replacements: ArrayList<ArrayList<Integer>>
		// labels (threshold): String []
		// distances (threshold): Matrix n x 1
		// X (threshold): Matrix n x m
	}
	
	private List<Object> pruneXMatrix (byte X [][]) {
		int	numUsed = 0,
			n = X.length,
			m = X [0].length,
			index = 0,
			indexMapping [] = new int [m];
		boolean	used [] = new boolean [m];
		Matrix	resMat = null;
		List<Object> resList = new ArrayList<Object> ();
		
		for (int j = 0; j < m; j++) {
			used [j] = false;
			for (int i = 0; i < n; i++) {
				if (X [i][j] > 0) {
					used [j] = true;
					numUsed++;
					break;
				}
			}
		}
		
		for (int i = 0; i < m; i++) indexMapping [i] = -1;
		
		resMat = new Matrix (n, numUsed);
		for (int j = 0; j < m; j++) {
			if (!used [j]) continue;
			for (int i = 0; i < n; i++) {
				resMat.set(i,index,X [i][j]);
				indexMapping [j] = index;
			}
			index++;
		}
		
		// remove columns explained by others
		// collect ones and sort
		m = resMat.getColumnDimension();
		n = resMat.getRowDimension();
		int	order [] = new int [m],
			sortedByOnes [] = new int [m],
			tmpVal = 0;
		boolean	removed [] = new boolean [m];
		
		for (int i = 0; i < m; i++) {
			order [i] = i;
			sortedByOnes [i] = 0;
			removed [i] = false;
			for (int j = 0; j < n; j++) {
				if (resMat.get(j, i) != 0.0) sortedByOnes [i]++;
			}
		}
		// sort
		for (int i = 0; i < (m-1); i++) {
			for (int j = m-1; j > i; j--) {
				if (sortedByOnes [j] > sortedByOnes [j-1]) {
					// value
					tmpVal = sortedByOnes [j];
					sortedByOnes [j] = sortedByOnes [j-1];
					sortedByOnes [j-1] = tmpVal;
					// place
					tmpVal = order [j];
					order [j] = order [j-1];
					order [j-1] = tmpVal;
				}
			}
		}
		
		ArrayList<Integer>	curRange = null,
							finalRange = null;
		ArrayList<ArrayList<Integer>>	replacements = new ArrayList<ArrayList<Integer>> ();
		int	curIndex = 0,
			curMissing = 0;
		boolean	add = false,
				explains [] = new boolean [m],
				explained [] = new boolean [m];
		for (int i = 0; i < m; i++) {
			curRange = new ArrayList<Integer> ();
			curIndex = i+1;
			curMissing = sortedByOnes [i];
			while (true) {
				// end of the list?
				if (curIndex >= m) {
					// empty range???
					if (curRange.size() == 0) {
						break;
					}
					// go back to last point
					else {
						curIndex = curRange.get(curRange.size()-1);
						curRange.remove(curRange.size()-1);
						curMissing -= sortedByOnes[curIndex++];
						continue;
					}
				}
				// explains already a different variable?
				if (explains [curIndex]) {
					curIndex++;
					continue;
				}
				if (sortedByOnes [curIndex] <= curMissing) {
					// check for overlap
					add = true;
					for (int k = 0; k < n; k++) {
						if (!add) break;
						// check column to replace
						if (resMat.get(k, order[i]) == 0.0 && resMat.get(k,order[curIndex]) == 1.0) {
							curIndex++;
							add = false;
							break;
						}
						
						// check already chosen columns
						for (int l = 0; l < curRange.size(); l++) {
							if (resMat.get(k, order[curRange.get(l)]) == 1.0 && resMat.get(k,order[curIndex]) == 1.0) {
								curIndex++;
								add = false;
								break;
							}
						}
					}
					if (add) {
						for (int j = 0; j < curRange.size(); j++) explains [curRange.get(j)] = true;
						curRange.add(curIndex);
						curMissing -= sortedByOnes [curIndex++];
					}
					if (curMissing == 0) {
						finalRange = new ArrayList<Integer> ();
						// need original position for re-introduction in post-processing
						for (int j = 0; j < indexMapping.length; j++) if (indexMapping [j] == order[i]) {
							finalRange.add(j);
							break;
						}
						for (int j = 0; j < curRange.size(); j++) {
							for (int k = 0; k < indexMapping.length; k++) if (indexMapping [k] == order[curRange.get(j)]) {
								finalRange.add(k);
								break;
							}
						}
						replacements.add(finalRange);
						explained [order[i]] = true;
						break;
					}
				}
				else {
					curIndex++;
				}
			}
		}
		
		// set new matrix
		Matrix	finalMat = new Matrix (n, numUsed-replacements.size());
		index = 0;
		for (int i = 0; i < m; i++) {
			if (explained [i]) continue;
			for (int j = 0; j < n; j++) finalMat.set(j,index,resMat.get(j, i));
			index++;
		}
		// change indexMapping
		int	shift = 0;
		for (int i = 0; i < indexMapping.length; i++) {
			if (indexMapping [i] == -1) continue;
			if (explained [indexMapping[i]]) {
				indexMapping [i] = -1;
				shift++;
				continue;
			}
			if (indexMapping[i] != -1) indexMapping [i] -= shift;
		}
		
		resList.add(finalMat);
		resList.add(indexMapping);
		resList.add(replacements);
		
		return resList;
	}
	
	private Matrix retrieveFullVector (Matrix v, int [] indexMapping, ArrayList<ArrayList<Integer>> replacements) {
		Matrix	vFull = new Matrix (indexMapping.length,1);
		ArrayList<Integer>	tmpList = null;
		boolean	isNeg = false,
				sameSign = true;
		double	commonValue = 0.0;
		
		// fill with inferred values
		for (int i = 0; i < indexMapping.length; i++) {
			if (indexMapping [i] == -1) vFull.set(i, 0, Double.NaN);
			else vFull.set (i, 0, v.get(indexMapping[i], 0));
		}
		
		// resolve replacements
		for (int i = replacements.size()-1; i >= 0; i--) { // need to reverse (go from columns with less ones to those with more ones
			tmpList = replacements.get(i);
			// only one element?
			if (tmpList.size() == 2) {
				vFull.set(tmpList.get(0), 0, vFull.get(tmpList.get(1), 0)/2.0);
				vFull.set(tmpList.get(1), 0, vFull.get(tmpList.get(1), 0)/2.0);
				continue;
			}
			// check sign
			isNeg = vFull.get(tmpList.get(1), 0) < 0;
			sameSign = true;
			// same sign??
			if (sameSign) {
				// get common value
				commonValue = (isNeg ? -1.0 : 1.0) * Double.MAX_VALUE;
				for (int j = 1; j < tmpList.size(); j++) {
					commonValue = (isNeg ? Math.max (commonValue, vFull.get(tmpList.get(j), 0)) : Math.min (commonValue, vFull.get(tmpList.get(j), 0)));
				}
				// set new values
				if (commonValue != 0.0) {
					vFull.set(tmpList.get(0), 0, commonValue);
					for (int j = 1; j < tmpList.size(); j++) {
						vFull.set(tmpList.get(j), 0, vFull.get(tmpList.get(j), 0) - commonValue);
					}
				}
				else {
					vFull.set(tmpList.get(0), 0, 0.0);
				}
			}
		}
		
		return vFull;
	}
	
	// using BVLS (bounded variable least squares)
	public Matrix solveNNLS_external (Matrix X, Matrix d) {
		int	n = X.getRowDimension(),
			m = X.getColumnDimension();
		float	A [] = new float [n*m],
				b [] = new float [n],
				bnd [] = new float [m*2],
				resVec [];
		Matrix	v = new Matrix(m,1);
		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				A [(j*n) + i] = (float) X.get (i,j);
			}
			b [i] = (float) d.get (i,0);
		}
		for (int i = 0; i < m; i++) {
			bnd [(i*2)] = 0.0f;
			bnd [(i*2)+1] = Float.MAX_VALUE;
		}
		
		resVec = solveBVLS(A, b, bnd, n, m);
		
		for (int i = 0; i < m; i++) v.set (i, 0, (double) resVec[i]);
		
		return v;
	}
	
	private List<Object> transformDistanceMatrix (HIMat data) {
		String	labels [] = null;
		int	numDist = 0,
			c = 0;
		double distVec [] = null;
		boolean	isThreshold [] = null;
		HashSet<String>	isSerum = new HashSet<String>(),
						isAntigen = new HashSet<String>();
		List<Object>	resList = new ArrayList<Object>();
		
		// count distances
		for (int i = 0; i < data.dim; i++) {
			for (int j = 0; j < data.dim; j++) if (!Double.isNaN(data.table[i][j])) numDist++;
		}
		
		// initialize arrays
		distVec = new double[numDist];
		isThreshold = new boolean[numDist];
		labels = new String[numDist];
		
		// collect data
		for (int i = 0; i < data.dim; i++) {
			for (int j = 0; j < data.dim; j++) {
				if (Double.isNaN(data.table[i][j])) continue; // skip non-existent values
				
				distVec[c] = data.table[i][j];
				isThreshold[c] = data.isThreshold[i][j];
				labels[c] = data.names[i] + "-" + data.names[j];
				
				c++;
			}
		}
		
		for (int  i = 0; i < data.dim; i++) {
			if (data.isAntigen[i]) isAntigen.add(data.names[i]);
			if (data.isSerum[i]) isSerum.add(data.names[i]);
		}
		
		resList.add(labels);
		resList.add(new Matrix(distVec,numDist));
		resList.add(isThreshold);
		resList.add(isSerum);
		resList.add(isAntigen);
		
		return resList;
	}
	
	private String [] collectMutations_updown (node curNode, String mutations []) {
		if (curNode.parent != null) {
			mutations [curNode.branchIndex_up] = mutations [mutations.length/2 + curNode.branchIndex_up] = curNode.mutations.toString();
		}
		for (int i = 0; i < curNode.numChilds; i++) mutations = this.collectMutations_updown (curNode.childArray [i], mutations);
		return mutations;
	}
	
	private void getMutationImpact (node curNode, Matrix v, int indexMapping [], String fileName) {
		String	mutations [] = null;
		int	dim = indexMapping.length;
		
		mutations = this.collectMutations_updown (curNode, new String [dim]);
		
		try {
			BufferedWriter	bw = new BufferedWriter (new FileWriter (new File (fileName + ".mutationImpact")));
			String	newLine = System.getProperty("line.separator");
			
			bw.write("branch_id\tweight\tmutation_set" + newLine);
			
			for (int i = 0; i < mutations.length; i++) {
				bw.write((i < indexMapping.length/2 ? i : -(i-indexMapping.length/2)) + "\t" + v.get(i,0) + "\t" + mutations [i] + newLine);
			}
			bw.close();
		}
		catch (IOException e) {
			System.out.println ("Caught the following exception: " + e);
		}
	}
	
	private void setBranchLengths_print (node curNode, Matrix v) {
		if (curNode.parent != null) {
			double	up_dist = v.get(curNode.branchIndex_up, 0),
					down_dist = v.get(v.getRowDimension()/2 + (-1) * curNode.branchIndex_down, 0);
			String	up_label = (Double.isNaN(up_dist) ? "-" : "" + Math.round(up_dist*1000.0)/1000.0),
					down_label = (Double.isNaN(down_dist) ? "-" : "" + Math.round(down_dist*1000.0)/1000.0);
			
			up_dist = (Double.isNaN(up_dist) ? 0.0 : up_dist);
			down_dist = (Double.isNaN(down_dist) ? 0.0 : down_dist);
			
			curNode.distance = Math.max(Math.abs(up_dist), Math.abs(down_dist));
			curNode.distLabel = up_label + " / " + down_label;
			if (up_dist < down_dist) curNode.negative = true;
		}
		
		for (int i = 0; i < curNode.numChilds; i++) this.setBranchLengths_print (curNode.childArray [i], v);
	}
	
	private double[] looError (Matrix X, Matrix d, int [] indexMapping, ArrayList<ArrayList<Integer>> replacements) {
		Matrix	X_tmp = null,
				d_tmp = null,
				v_tmp = null,
				d_test = null;
		int	n = X.getRowDimension(),
			m = X.getColumnDimension();
		double	error [] = new double [2+d.getRowDimension()],
				ratio = ((double)n) / 40.0;
		
		X_tmp = new Matrix (n-1,m);
		d_tmp = new Matrix (n-1,1);
		
		for (int i = 1; i < n; i++) {
			for (int j = 0; j < m; j++) {
				X_tmp.set (i-1,j, X.get (i,j));
			}
			d_tmp.set (i-1,0, d.get(i, 0));
		}
		
		System.out.println();
		System.out.println(" 0%                50%               100%");
		System.out.println(" |..................|...................|");
		System.out.print(" ");
		
		for (int i = 0; i < n; i++) {
			v_tmp = this.solveNNLS_external(X_tmp,d_tmp);
			//v_tmp = this.solveNNLS(X_tmp,d_tmp);
			
			d_test = X.times(v_tmp);
			error [2+i] = d_test.get(i, 0);
			error [0] += (Math.abs(d_test.get(i, 0) - d.get(i,0)));
			error [1] += (Math.pow(d_test.get(i, 0) - d.get(i,0),2));
			
			if (i < (n-1)) {
				for (int k = 0; k < m; k++) {
					X_tmp.set (i, k, X.get(i, k));
				}
				d_tmp.set (i,0,d.get(i, 0));
			}
			
			for (int j = 0; j < (Math.floor(((double)(i+1))/ratio) - Math.floor(((double)i)/ratio)); j++ ) System.out.print("*");
		}
		System.out.println();
		
		error [0] = (error [0]/((double)n));
		error [1] = Math.sqrt(error [1]/((double)n));
		
		return error;
	}
	
	private double[] cvError (Matrix X, Matrix d, String settings, int [] indexMapping, ArrayList<ArrayList<Integer>> replacements) {
		String	tmpSplit [] = settings.split(",");
		Matrix	X_tmp = null,
				d_tmp = null,
				v_tmp = null,
				d_test = null;
		int	n = X.getRowDimension(),
			m = X.getColumnDimension(),
			parts = Integer.parseInt(tmpSplit [0]),
			times = (tmpSplit.length > 1 ? Integer.parseInt(tmpSplit [1]) : 1),
			indexSample [] = null,
			cvSplit [] = null,
			curParts [] = null,
			curInd = 0,
			indOld = 0,
			numInPart = (int) Math.floor((double) n / (double) parts);
		double	error [] = new double [2+n],
				current = 0.0,
				ratio = ((double) (times*parts)) / 40.0;
		Random	r = new Random();
		
		X_tmp = new Matrix (n-numInPart,m);
		d_tmp = new Matrix (n-numInPart,1);
		curParts = new int [n-numInPart];
		
		indexSample = new int [n];
		cvSplit = new int [n];

		System.out.println();
		System.out.println(" 0%                50%               100%");
		System.out.println(" |..................|...................|");
		System.out.print(" ");
		
		for (int t = 0; t < times; t++) {
		
		// set up indices
		for (int i = 0; i < n; i++) {
			indexSample [i] = i % parts;
			cvSplit [i] = -1;
		}
		curInd = 0;
		for (int i = 0; i < n; i++) {
			// randomly select an index
			curInd = r.nextInt(n-i);
			cvSplit [i] = indexSample [curInd];
			// replace chosen element by last one
			if (curInd < (n-i-1)) indexSample [curInd] = indexSample [n-i-1];
		}
		
		curInd = 0;
		
		// set starting values
		curInd = 0;
		for (int i = 0; i < n; i++) {
			if (cvSplit [i] == 0) continue;
			for (int j = 0; j < m; j++) {
				X_tmp.set (curInd,j, X.get (i,j));
			}
			d_tmp.set (curInd,0, d.get(i, 0));
			curParts [curInd] = cvSplit [i];
			curInd++;
		}
		if (curInd < (n-numInPart)) {
			for (int i = curInd; i < (n-numInPart); i++) {
				d_tmp.set(i, 0, 0.0);
				for (int j = 0; j < m; j++) X_tmp.set (i,j, 0.0);
				curParts [i] = 1;
			}
		}
		
		for (int i = 0; i < parts; i++) {
			v_tmp = this.solveNNLS_external(X_tmp,d_tmp);
			//v_tmp = this.solveNNLS(X_tmp,d_tmp);
			
			d_test = X.times(v_tmp);
			for (int j = 0; j < n; j++) {
				if (cvSplit [j] == i) {
					error [2+j] += d_test.get(j, 0);
					error [0] += (Math.abs(d_test.get(j, 0) - d.get(j,0)));
					error [1] += (Math.pow(d_test.get(j, 0) - d.get(j,0),2));
				}
			}
			
			// replace elements
			if (i < (parts-1)) {
				indOld = 0;
				for (int k = 0; k < (n-numInPart); k++) {
					if (curParts [k] == (i+1)) {
						while (indOld < n) {
							if (cvSplit [indOld] == i) {
								curParts [k] = i;
								d_tmp.set(k, 0, d.get(indOld, 0));
								for (int j = 0; j < m; j++) X_tmp.set(k, j, X.get(indOld, j));
								indOld++;
								break;
							}
							indOld++;
						}
						// set unused elements to zero
						if (indOld == n && curParts [k] == (i+1)) {
							curParts [k] = (i+2);
							d_tmp.set(k, 0, 0.0);
							for (int j = 0; j < m; j++) X_tmp.set(k, j, 0.0);
						}
					}
				}
			}
			
			// System.out.print("\b");
			for (int j = 0; j < (Math.floor((current+1.0)/ratio) - Math.floor(current/ratio)); j++ ) System.out.print("*");
			current++;
		}
		
		}
		System.out.println();
		
		for (int i = 0; i < n; i++) error [2+i] /= ((double) times);

		error [0] = (error [0]/((double)(n*times)));
		error [1] = Math.sqrt(error [1]/((double)(n*times)));
		
		return error;
	}
	
	private void setLabelSuffix (node curNode, HashSet<String> isAntigen, HashSet<String> isSerum) {
		if (curNode.numChilds == 0) {
			if (isAntigen.contains(curNode.identifier) && isSerum.contains(curNode.identifier)) curNode.label += " [AG/SR]";
			else if (isAntigen.contains(curNode.identifier)) curNode.label += " [AG]";
			else if (isSerum.contains(curNode.identifier)) curNode.label += " [SR]";
			else curNode.label += " []";
		}
		else {
			for (int i = 0; i < curNode.numChilds; i++) this.setLabelSuffix (curNode.childArray [i], isAntigen, isSerum);
		}
	}
	
	/*
	public Matrix solveNNLS (Matrix X, Matrix d) {
		// solve ||Xv - d||
		int	xm = X.getColumnDimension (),
			xn = X.getRowDimension(),
			maxIter = 1000,
			iter = 0,
			curMax = -1,
			tmpInd = 0,
			maxPos = 0,
			svd_p = 0;
		double	epsilon = 1e-10,
				max = Double.MAX_VALUE,
				tmpVal = 0.0,
				alpha = Double.MAX_VALUE,
				tmpAlpha = Double.MAX_VALUE,
				lambda = 1e-10;
		Matrix	v = new Matrix (xm,1),
				w = null,
				Xpos = null,
				X_t = X.transpose(),
				zprime = null,
				z = null;
		SingularValueDecomposition	Xpos_svd = null;
		List<Integer>	P = new ArrayList<Integer>(),
						Z = new ArrayList<Integer>();
		boolean	allPos = false;
		
		// STEP I: initialize z, p is empty
		for (int i = 0; i < xm; i++) Z.add(i);
		
		while ((iter++) < maxIter) {
			switch (iter%4) {
				case 0: System.out.print("|"); break;
				case 1: System.out.print("/"); break;
				case 2: System.out.print("-"); break;
				case 3: System.out.print("\\"); break;
			}
			
			// STEP II: compute the gradient vector w = Xt(d-Xv)
			w = X_t.times(d.minus(X.times(v)));
			
			// STEP III
			if (checkGradient (Z, w, epsilon)) {
				System.out.print("\b");
				break;
			}
			
			// STEP IV: find maximum element of w and move the index to P
			curMax = Z.get(0);
			max = w.get (curMax, 0);
			maxPos = 0;
			for (int i = 1; i < Z.size(); i++) {
				tmpInd = Z.get(i);
				tmpVal = w.get(tmpInd, 0);
				if (tmpVal > max) {
					curMax = tmpInd;
					max = tmpVal;
					maxPos = i;
				}
			}
			P.add(curMax);
			Z.remove(maxPos);
			
			// secondary loop
			while (true) {
				// STEP V: calculate z
				Xpos = new Matrix (xn, P.size());
				for (int i = 0; i < P.size(); i++) {
					tmpInd = P.get(i);
					for (int j = 0; j < xn; j++) Xpos.set(j,i,X.get(j, tmpInd));
				}
				
				// Matrix is singular therefore solve explicitly via ridge regression
				Xpos_svd = Xpos.svd();
				svd_p = Xpos_svd.getU().getColumnDimension();
				// V(S^2+lambdaI)-1SUtd
				zprime = (Xpos_svd.getV().times((((Xpos_svd.getS()).times(Xpos_svd.getS())).plus(Matrix.identity(svd_p, svd_p).times(lambda))).inverse())).times(Xpos_svd.getS().times(Xpos_svd.getU().transpose().times(d)));
				
				z = new Matrix (xm,1);
				allPos = true;
				for (int i = 0; i < P.size(); i++) {
					tmpInd = P.get(i);
					z.set(tmpInd, 0, zprime.get(i, 0));
					if (z.get(tmpInd, 0) <= 0) allPos = false; 
				}
				
				// STEP VI: if all elements in z with indices in P are positive
				if (allPos) {
					v = z;
					break;
				}
				
				// STEP VII
				alpha = Double.MAX_VALUE;
				for (int i = 0; i < z.getRowDimension(); i++) {
					if (z.get(i, 0) < 0) {
						tmpAlpha = ((double) v.get(i, 0)) / ((double) (v.get(i, 0) - z.get(i, 0)));
						if (tmpAlpha < alpha) alpha = tmpAlpha;
					}
				}
				
				// STEP VIII: calculate v = v + alpha (z - v)
				v = v.plus(z.minus(v).times(alpha));
				
				// STEP IX:
				for (int i = 0; i < P.size(); i++) {
					tmpInd = P.get(i);
					if (Math.abs(v.get(tmpInd, 0)) < epsilon) {
						P.remove(i);
						Z.add(tmpInd);
						i--;
					}
				}
			}
			System.out.print("\b");
		}
		
		return v;
	}
	
	private boolean checkGradient (List<Integer> Z, Matrix w, double threshold) {
		// check if z is empty or all elements of w with indices in z have values <= threshold (0 or very close to zero)
		if (Z.isEmpty()) return true;
		
		for (int i = 0; i < Z.size(); i++) if (w.get(((Integer) Z.get(i)).intValue(), 0) > threshold) return false;
		return true;
	}
	*/
	
	// driver ----------------------------------------- //
	
	public void computeBranchWeights (node curNode, HIMat him, String fileName, String loo, String cv, boolean map, String [][] leaveSeqs, int numNodes) {
		List<Object>	d_data = this.transformDistanceMatrix (him), // ArrayList; 0 = labels (String[]), 1 = distance matrix (Matrix), 2 = isThreshold (boolean[]), 3 = isSerum (HasHSet<String>), 4 = isAntigen (HashSet<String>)
						processed_data = this.getMainData (curNode, d_data);
		Matrix	v = null,
				dprime = null;
		double	squaredSum = 0.0,
				absoluteSum = 0.0,
				looError [] = new double [2];
		treeObj	tree = new treeObj ();
		String	newLine = System.getProperty("line.separator");
		
		// set label suffix
		this.setLabelSuffix (curNode, (HashSet<String>) d_data.get(4), (HashSet<String>) d_data.get(3));
		
		// compute solution
		v = this.solveNNLS_external ((Matrix) processed_data.get(2), (Matrix) processed_data.get(1));
		//v = solveNNLS((Matrix) processed_data.get(2), (Matrix) processed_data.get(1));
		
		// predict values
		dprime = ((Matrix) processed_data.get(2)).times(v);
		
		if (loo != null) {
			looError = this.looError ((Matrix) processed_data.get(2), (Matrix) processed_data.get(1), (int []) processed_data.get(3), (ArrayList<ArrayList<Integer>>) processed_data.get(4));
		}
		if (cv != null) {
			looError = this.cvError ((Matrix) processed_data.get(2), (Matrix) processed_data.get(1), cv, (int []) processed_data.get(3), (ArrayList<ArrayList<Integer>>) processed_data.get(4));
		}
		
		// complete solution
		v = this.retrieveFullVector(v, (int []) processed_data.get(3), (ArrayList<ArrayList<Integer>>) processed_data.get(4));
		
		// output tree
		this.setBranchLengths_print(curNode, v);
		tree.writeNexusTree(tree.getNewick(curNode,map).toString()+";",fileName + ".withMuts.tre",leaveSeqs,curNode,numNodes,map);
		
		// map mutations with weights
		this.getMutationImpact(curNode, v, (int []) processed_data.get(3), fileName);
		
		// write results to files
		try {
			BufferedWriter	bw_d = new BufferedWriter (new FileWriter (new File (fileName + ".distance")));
			Matrix	d_known = (Matrix) processed_data.get(1);
			String	labels_known [] = (String []) processed_data.get(0);
			
			bw_d.write("label\ttrue_val\ttrain_val\ttrain_err" + (loo != null || cv != null ? "\ttest_val\ttest_err" : "") + newLine);
			
			for (int  i = 0; i < dprime.getRowDimension(); i++) {
				// compute error
				squaredSum += (Math.pow((dprime.get(i, 0)) - (d_known.get(i, 0)),2));
				absoluteSum += Math.abs((dprime.get(i, 0)) - (d_known.get(i, 0)));
				// write to file
				bw_d.write(labels_known[i] + "\t" + d_known.get(i, 0) + "\t" + dprime.get(i, 0) + "\t" + Math.pow(dprime.get(i, 0) - d_known.get(i, 0),2) + (loo != null || cv != null ? "\t" + looError[2+i] + "\t" + Math.pow(looError[2+i] - d_known.get(i, 0),2) : "") + newLine);
			}
			
			bw_d.write("# squared " + squaredSum + newLine);
			bw_d.write("# absolute " + absoluteSum + newLine);
			if (loo != null || cv != null) {
				bw_d.write("# loo squared " + looError[1] + newLine);
				bw_d.write("# loo absolute " + looError[0] + newLine);
			}
			bw_d.close();
		}
		catch (IOException e) {
			System.out.println("The following exception was caught during file output: " + e);
			System.exit(1);
		}
	}
}
