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

import java.lang.reflect.Array;
import java.util.*;


interface TreeRecursiveUpdate {

	void pre(node n);
	void post(node n);
}

public class Sankoff {
	
	// private void fitchPhaseOne (node curNode, int pos, boolean gaps) {
	// 	if (curNode.numChilds == 0) {
	// 		curNode.intermAnc = new HashSet<String>();
	// 		curNode.fitchType = "child";
			
	// 		if (curNode.type.equals("DNA")) {
	// 				if (curNode.tmpSequence.charAt(pos) != '-') {
	// 					switch (curNode.tmpSequence.charAt(pos)) {
	// 						case 'R': case 'r': curNode.intermAnc.add("A"); curNode.intermAnc.add("G"); break;
	// 						case 'Y': case 'y': curNode.intermAnc.add("C"); curNode.intermAnc.add("T"); break;
	// 						case 'M': case 'm': curNode.intermAnc.add("A"); curNode.intermAnc.add("C"); break;
	// 						case 'K': case 'k': curNode.intermAnc.add("G"); curNode.intermAnc.add("T"); break;
	// 						case 'S': case 's': curNode.intermAnc.add("C"); curNode.intermAnc.add("G"); break;
	// 						case 'W': case 'w': curNode.intermAnc.add("A"); curNode.intermAnc.add("T"); break;
	// 						case 'H': case 'h': curNode.intermAnc.add("A"); curNode.intermAnc.add("C"); curNode.intermAnc.add("T"); break;
	// 						case 'B': case 'b': curNode.intermAnc.add("C"); curNode.intermAnc.add("G"); curNode.intermAnc.add("T"); break;
	// 						case 'V': case 'v': curNode.intermAnc.add("A"); curNode.intermAnc.add("C"); curNode.intermAnc.add("G"); break;
	// 						case 'D': case 'd': curNode.intermAnc.add("A"); curNode.intermAnc.add("G"); curNode.intermAnc.add("T"); break;
	// 						case 'N': case 'n': case '?': curNode.intermAnc.add("A"); curNode.intermAnc.add("C"); curNode.intermAnc.add("G"); curNode.intermAnc.add("T"); if (gaps) curNode.intermAnc.add("-"); break; // if (gaps) curNode.intermAnc.add("-");
	// 						default: curNode.intermAnc.add(""+Character.toUpperCase(curNode.tmpSequence.charAt(pos))); break;
	// 					}
	// 				}
	// 				else {
	// 					if (gaps) {
	// 						curNode.intermAnc.add(""+Character.toUpperCase(curNode.tmpSequence.charAt(pos)));
	// 					}
	// 					else {
	// 						curNode.intermAnc.add("A"); 
	// 						curNode.intermAnc.add("C"); 
	// 						curNode.intermAnc.add("G"); 
	// 						curNode.intermAnc.add("T");
	// 					}
	// 				}
	// 			}
	// 			else {
	// 				if (curNode.tmpSequence.charAt(pos) != '-') {
	// 					switch (curNode.tmpSequence.charAt(pos)) {
	// 						case 'B': case 'b': curNode.intermAnc.add("D"); curNode.intermAnc.add("N"); break;
	// 						case 'Z': case 'z': curNode.intermAnc.add("E"); curNode.intermAnc.add("Q"); break;
	// 						case 'J': case 'j': curNode.intermAnc.add("L"); curNode.intermAnc.add("I"); break;
	// 						case 'X': case 'x': case '?': for (int i = 65; i < 91; i++) { curNode.intermAnc.add("" + (char)i); } break; //if (gaps) curNode.intermAnc.add("-"); 
	// 						default: curNode.intermAnc.add(""+Character.toUpperCase(curNode.tmpSequence.charAt(pos))); break;
	// 					}
	// 				}
	// 				else {
	// 					if (gaps) {
	// 						curNode.intermAnc.add(""+Character.toUpperCase(curNode.tmpSequence.charAt(pos)));
	// 					}
	// 					else {
	// 						for (int i = 65; i < 91; i++) { curNode.intermAnc.add("" + (char)i); }
	// 					}
	// 				}
	// 			}
	// 	}
	// 	else {
	// 		for (int i = 0; i < curNode.numChilds; i++) this.fitchPhaseOne (curNode.childArray [i], pos, gaps);
	// 		if (curNode.numChilds == 2) {
	// 			curNode.intermAnc = new HashSet<String>(curNode.childArray[0].intermAnc);
	// 			// if one is empty build the union
	// 			if (curNode.childArray[0].intermAnc.isEmpty() || curNode.childArray[1].intermAnc.isEmpty()) {
	// 				curNode.intermAnc.addAll (curNode.childArray[1].intermAnc);
	// 				curNode.fitchType = (curNode.childArray[0].intermAnc.isEmpty() ? curNode.childArray[0].fitchType : curNode.childArray[1].fitchType);
	// 			}
	// 			else {
	// 				curNode.intermAnc.retainAll(curNode.childArray[1].intermAnc);	// intersection
	// 				if (curNode.intermAnc.isEmpty()) {
	// 					curNode.intermAnc.addAll(curNode.childArray[0].intermAnc);	// union
	// 					curNode.intermAnc.addAll(curNode.childArray[1].intermAnc);	// union
	// 					curNode.fitchType = "union";
	// 				}
	// 				else {
	// 					curNode.fitchType = "intersection";
	// 				}
	// 			}
	// 		}
	// 		else {	// multifurcating tree => majority voting
	// 			curNode.fitchType = "multi";
	// 			curNode.intermAnc = new HashSet<String>();
				
	// 			int	intermCount [] = new int [27];
	// 			int	max = 1,
	// 				numSet = 0;
				
	// 			for (int i = 0; i < curNode.numChilds; i++) {
	// 				if (!curNode.childArray [i].intermAnc.isEmpty()) {
	// 					for (String s : curNode.childArray[i].intermAnc) {
	// 						if (s.equals("-")) {
	// 								if(gaps) intermCount[intermCount.length-1]++;
	// 						}
	// 						else if (s.equals("?")) {
								
	// 						}
	// 						else {
	// 							intermCount[((int)s.charAt(0)) -65]++;
	// 						}
	// 					}
	// 				}
	// 			}
				
	// 			for(int i = 0; i < intermCount.length; i++) {
	// 				if (intermCount[i] > 0) numSet++;
	// 				if (intermCount [i] > max) {
	// 					curNode.intermAnc = new HashSet<String>();
	// 					curNode.intermAnc.add("" + (char) (i + (i==(intermCount.length-1) ? 19 : 65)));
	// 					max = intermCount [i];
	// 				}
	// 				else if (intermCount [i] == max) curNode.intermAnc.add("" + (char) (i + (i==(intermCount.length-1) ? 19 : 65)));
	// 			}
	// 		}
	// 	}
	// }
	
	// private Set<String> [] checkDownwards (node curNode) {
	// 	Set<String>	returnArray [] = new Set [2];
	// 	returnArray [0] = new HashSet<String>();
	// 	returnArray [1] = null;
		
	// 	if (curNode.numChilds == 0) {
	// 		returnArray [0].add("child");
	// 	}
	// 	else {
	// 		boolean	isNonEmptyChild [] = new boolean [curNode.numChilds];
	// 		int	sum = 0;
			
	// 		for (int i = 0; i < curNode.numChilds; i++) {
	// 			if (!curNode.childArray[i].intermAnc.isEmpty()) {
	// 				sum++;
	// 				isNonEmptyChild [i] = true;
	// 			}
	// 		}
			
	// 		if (sum == 1) {
	// 			for (int i = 0; i < isNonEmptyChild.length; i++) {
	// 				if (isNonEmptyChild [i]) {
	// 					returnArray = this.checkDownwards(curNode.childArray[i]);
	// 					break;
	// 				}
	// 			}
	// 		}
	// 		else if (sum == 2) {
	// 			for (int i = 0; i < isNonEmptyChild.length; i++) {
	// 				if (isNonEmptyChild[i]) {
	// 					if (returnArray [1] == null) returnArray [1] = new HashSet<String>(curNode.childArray[i].intermAnc);
	// 					else {
	// 						Set<String> tmpSet = new HashSet<String>(returnArray[1]);
	// 						tmpSet.retainAll(curNode.childArray[i].intermAnc);
	// 						if (tmpSet.isEmpty()) returnArray [0].add("union");
	// 						else returnArray [0].add("intersection");
	// 						returnArray [1].addAll(curNode.childArray[i].intermAnc);
	// 						break;
	// 					}
	// 				}
	// 			}
	// 		}
	// 		else if (sum > 2) {
	// 			returnArray[0].add("multi");
	// 		}
	// 	}
		
	// 	return returnArray;
	// }
	
	// private void fitchPhaseTwo (node curNode, int pos, boolean gaps) {
	// 	if (curNode.numChilds > 0) {
	// 		if (curNode.parent == null) for (int i = 0; i < curNode.numChilds; i++) fitchPhaseTwo (curNode.childArray[i], pos, gaps);
	// 		else if (!curNode.intermAnc.isEmpty()) {
	// 			Set<String>	finalSet = new HashSet<String>();
	// 			for (String sc : curNode.intermAnc) {
	// 				for (String sp : curNode.parent.intermAnc) {
	// 					if (sc.equals(this.toUpper(sp))) {
	// 						finalSet.add(sc);
	// 						break;
	// 					}
	// 				}
	// 			}
				
	// 			if (curNode.numChilds > 2) {
	// 				if (!finalSet.isEmpty()) curNode.intermAnc = new HashSet<String>(finalSet);
	// 			}
	// 			else {
	// 				boolean	fixed = true;
	// 				for (String s : curNode.parent.intermAnc) {
	// 					if (!finalSet.contains(s) && !finalSet.contains(this.toUpper(s))) {
	// 						fixed = false;
	// 						break;
	// 					}
	// 				}
	// 				if (fixed) {
	// 					curNode.intermAnc = new HashSet<String>(finalSet);
	// 				}
	// 				else {
	// 					if (curNode.fitchType.equals("union")) {
	// 						for (String s : curNode.parent.intermAnc) {
	// 							if (!curNode.intermAnc.contains(this.toUpper(s))) {
	// 								curNode.intermAnc.add(this.toLower(s));
	// 							}
	// 						}
	// 					}
	// 					else if (curNode.fitchType.equals("child")) {
	// 						;
	// 					}
	// 					else {
	// 						Set<String>	down [] = checkDownwards (curNode); // introduced
	// 						if (!down[0].contains("child") && !down[0].contains("multi")) {
	// 							for (String s : curNode.parent.intermAnc) {
	// 								if (!curNode.intermAnc.contains(this.toUpper(s)) && down[1].contains(this.toUpper(s))) curNode.intermAnc.add(this.toLower(s));
	// 							}
	// 						}
	// 					}
	// 				}
	// 			}
				
	// 			for (int i = 0; i < curNode.numChilds; i++) fitchPhaseTwo (curNode.childArray[i], pos, gaps);
	// 		}
	// 	}
	// }
	
	// private void reconstructAncestralStates (node curNode, String strategy, int pos, boolean gaps) {
	// 	String	finalSet [],
	// 		parentChar;
	// 	Random	r = new Random ();
	// 	int	ind = 0;
		
	// 	Set<String> alphabetSet = new HashSet<String> ();
	// 	for (int i = 65; i < 91; i++) {
	// 		alphabetSet.add(("" + (char)i).toUpperCase());
	// 	}
	// 	if (gaps) alphabetSet.add("-");
		
	// 	if (curNode.parent == null) {
	// 		if (curNode.intermAnc.isEmpty()) curNode.tmpSequence += "-";
	// 		else {
	// 			finalSet = new String [curNode.intermAnc.size()];
	// 			for (String s : curNode.intermAnc) finalSet[ind++] = s;
				
	// 			if (finalSet.length == 1) {
	// 				curNode.tmpSequence += this.toUpper(finalSet [0]);
	// 			}
	// 			else {
	// 				// choose randomly
	// 				curNode.tmpSequence += this.toUpper(finalSet [Math.abs(r.nextInt()) % finalSet.length]);
	// 			}
	// 		}
	// 	}
	// 	else if (curNode.numChilds > 0){
	// 		parentChar = this.toUpper("" + curNode.parent.tmpSequence.charAt(curNode.tmpSequence.length()));
	// 		if (curNode.intermAnc.isEmpty()) curNode.tmpSequence += parentChar;
	// 		else {
	// 			finalSet = new String [curNode.intermAnc.size()];
	// 			for (String s : curNode.intermAnc) finalSet[ind++] = s;
	// 			if (finalSet.length == 1) {
	// 				curNode.tmpSequence += this.toUpper(finalSet [0]);
	// 			}
	// 			else {
	// 				if(curNode.intermAnc.contains(parentChar) || curNode.intermAnc.contains(this.toLower(parentChar))) {
	// 					if (curNode.intermAnc.contains(parentChar)) curNode.tmpSequence += parentChar;
	// 					else {
	// 						if (strategy.equals("DelTran")) {
	// 							curNode.tmpSequence += parentChar;
	// 						}
	// 						else {
	// 							curNode.intermAnc.retainAll(alphabetSet);
	// 							ind = 0;
	// 							finalSet = new String [curNode.intermAnc.size()];
	// 							for (String s : curNode.intermAnc) finalSet[ind++] = s;
								
	// 							// choose randomly
	// 							curNode.tmpSequence += this.toUpper(finalSet [Math.abs(r.nextInt()) % finalSet.length]);
	// 						}
	// 					}
	// 				}
	// 				else {
	// 					ind = 0;
	// 					curNode.intermAnc.retainAll(alphabetSet);
	// 					finalSet = new String [curNode.intermAnc.size()];
	// 					for (String s : curNode.intermAnc) {
	// 						finalSet[ind++] = s;
	// 					}
						
	// 					// choose randomly
	// 					curNode.tmpSequence += this.toUpper(finalSet [Math.abs(r.nextInt()) % finalSet.length]);
	// 				}
	// 			}
	// 		}
	// 	}
		
	// 	if (curNode.numChilds > 0) for (int i = 0; i < curNode.numChilds; i++) this.reconstructAncestralStates(curNode.childArray[i], strategy, pos, gaps);

	// 	/*
	// 	Ni -> Ni || (Ni) -> Ni is obligatory
	// 	(Ni) -> (Ni) || Ni -> (Ni) is deltran
	// 	(Ni) -> Nj || Ni -> Nj is acctran
	// 	(Ni) -> (Nj) || Ni -> (Nj) is not permitted
	// 	*/
	// }
	
	// private String toUpper (String s) {
	// 	if (s.charAt(0) == '^') return "-";
	// 	else return s.toUpperCase();
	// }
	
	// private String toLower (String s) {
	// 	if (s.charAt(0) == '-') return "^";
	// 	else return s.toLowerCase();
	// }

	// private void assignFixedAncestrals (node curNode, int pos, String ancChar) {
	// 	if (curNode.numChilds != 0) {
	// 		for (int i = 0; i < curNode.numChilds; i++) this.assignFixedAncestrals (curNode.childArray [i], pos, ancChar);
	// 		curNode.tmpSequence += ancChar;
	// 	}
	// }
	
	// private void copyLeaveSeqs (node curNode) {
	// 	if (curNode.numChilds > 0) {
	// 		for (int i = 0; i < curNode.numChilds; i++) this.copyLeaveSeqs(curNode.childArray[i]);
	// 		curNode.tmpSequence = "";
	// 	}
	// 	else curNode.tmpSequence = curNode.mSequence;
	// }
	
	// private int [] checkSiteRecursion (node curNode, int pos, int counts [], boolean gaps) {
	// 	int	tmpCounts [] = new int [counts.length];
	// 	System.arraycopy(counts,0,tmpCounts,0,counts.length);
	// 	if (curNode.numChilds > 0) {
	// 		for (int i = 0; i < curNode.numChilds; i++) tmpCounts = checkSiteRecursion (curNode.childArray[i], pos, tmpCounts, gaps);
	// 	}
	// 	else {
	// 		try{
	// 			if (curNode.tmpSequence.charAt(pos) == '-' || curNode.tmpSequence.charAt(pos) == '/') {
	// 				if (gaps) tmpCounts [tmpCounts.length-2]++;
	// 			}
	// 			else {
	// 				if (curNode.tmpSequence.charAt(pos) != '?') tmpCounts [((int)Character.toUpperCase(curNode.tmpSequence.charAt(pos))) - 65]++;
	// 			}
	// 			tmpCounts [tmpCounts.length-1]++;
	// 		}
	// 		catch (Exception e) {
	// 			System.out.println(e);
	// 			System.out.println(curNode.identifier + " " + curNode.tmpSequence + " " + pos);
	// 			System.exit(1);
	// 		}
	// 	}
		
	// 	return tmpCounts;
	// }
	
	// private String checkSite (node curNode, int pos, boolean gaps) {
	// 	int	counts [] = new int [28],
	// 		sum = 0,
	// 		which = 0;
	// 	for (int i = 0; i < curNode.numChilds; i++) counts = checkSiteRecursion (curNode.childArray[i], pos, counts, gaps);
	// 	for (int i = 0; i < (counts.length-1); i++) {
	// 		if (counts [i] > 0) {
	// 			sum++;
	// 			which = i;
	// 		}
	// 	}
	// 	if (sum > 1) {
	// 		return "multi";
	// 	}
	// 	else {
	// 		return "" + (char)(which+(which == 26 ? 19: 65));
	// 	}
	// }
	
	// private void setIntermediates (node curNode, int pos) {
	// 	if (curNode.numChilds > 0) {
	// 		curNode.mSequence = (pos == -1 ? new String(curNode.tmpSequence) : curNode.mSequence.subSequence(0, pos) + curNode.tmpSequence + (pos < (curNode.mSequence.length()-1) ? curNode.mSequence.substring(pos+1) : ""));
	// 		curNode.tmpSequence = "";
	// 		curNode.intermAnc = null;
	// 	}
		
	// 	for (int i = 0; i < curNode.numChilds; i++) this.setIntermediates(curNode.childArray [i], pos);
	// }
	
	/**
	 * Fills all node tmpSequence parameter based on tmpSequence of the leaf nodes
	 * @param root
	 * @param seqLength
	 * @param gaps allow gap. If false, we consider gap as unknown.
	 */
	public void ancestralStateReconstruction(node root, int seqLength, boolean gaps) {
		// System.err.println("Length="+seqLength);
		// assert(gaps == true);
		final char GAP_CHARACTER = '-';
		final int INFINITY = 10000 * 2 + 1;

		// fill the options
		boolean[] charCount = new boolean[256];
		// fillCharCount(root, charCount);
		recursiveCall(root, new TreeRecursiveUpdate() {
			@Override
			public void pre(node curNode) {
				if (curNode.numChilds == 0) {
					//check if no bad character is in the sequences
					boolean[] allowedOptions = new boolean[256];
					if (curNode.type.equals("DNA")) {
						for (char c : "ATCG-".toCharArray()) {
							allowedOptions[(int)c] = true;
						}
					} else {
						for (char c : "ARNDCQEGHILKMFPSTWYV-".toCharArray()) {
							allowedOptions[(int)c] = true;
						}
					}
					for (char c : curNode.mSequence.toCharArray()) {
						if (!allowedOptions[Character.toUpperCase(c)]) {
							throw new RuntimeException("Invalid character in " + curNode.ID + " c='" + c + "' with type " + curNode.type);
						}
					}
					//end of check
		
					for (char c : curNode.mSequence.toCharArray()) {
						charCount[(int)c] = true;
					}
				}
			}

			@Override
			public void post(node n) {
			}
		});
		List<Character> optionsList = new ArrayList<Character>();

		for (int i = 0; i < charCount.length; i++) {
			char c = (char) i, cUpper = Character.toUpperCase(c);
			if (charCount[i] && cUpper != c) {
				charCount[i] = false;
				charCount[(int)cUpper] = true;
			}
			if (charCount[i]) {
				optionsList.add((char)i);
			}
		}

		int[] optionsIndex = new int[256];
		Arrays.fill(optionsIndex, -1);
		char[] options = new char[optionsList.size()];
		for (int i = 0; i < optionsList.size(); i++) {
			options[i] = optionsList.get(i);
			optionsIndex[(int) options[i]] = i;
		}

		// fill the costMatrix
		int[][] costMatrix = new int[options.length][options.length];
		for (int i = 0; i<options.length; i++) {
			for (int j = 0; j<options.length; j++) {
				if (i == j) {
					costMatrix[i][j] = 0;
				} else if (options[i] == GAP_CHARACTER || options[j] == GAP_CHARACTER) {
					if (options[i] == GAP_CHARACTER) {
						costMatrix[i][j] = 1;
					} else {
						costMatrix[i][j] = 1;
					}
				} else {
					costMatrix[i][j] = 1;
				}
			}
		}
		if (!gaps && optionsIndex[(int) GAP_CHARACTER] != -1) {
			int i = optionsIndex[(int) GAP_CHARACTER];
			for (int j = 0; j<options.length; j++) {
				costMatrix[i][j] = INFINITY;
			}
		}

		// System.err.print("\n  ");
		// for (int j = 0; j<options.length; j++) {
		// 	System.err.print(options[j] + " ");
		// }
		// System.err.println();
		// // for (int j = 0; j<options.length; j++) {
		// // 	System.err.print(optionsIndex[options[j]] + " ");
		// // }
		// // System.err.println();
		// for (int i = 0; i<options.length; i++) {
		// 	System.err.print(options[i] + " ");
		// 	for (int j = 0; j<options.length; j++) {
		// 		System.err.print(costMatrix[i][j] + " ");
		// 	}
		// 	System.err.println();
		// }

		// fill nodeIndex
		// TODO: use node.ID
		Map<node, Integer> nodeIndex = new HashMap<>();
		List<node> indexNode = new ArrayList<>();
		recursiveCall(root, new TreeRecursiveUpdate() {

			@Override
			public void pre(node n) {
			}

			@Override
			public void post(node n) {
				nodeIndex.put(n, indexNode.size());
				indexNode.add(n);
			}
			
		});

		// reset mSequence for intermediate nodes
		recursiveCall(root, new TreeRecursiveUpdate() {
			@Override
			public void pre(node n) {
				if (n.numChilds != 0) {
					// System.err.println("Reset " + n.ID + ": " + n.mSequence);
					n.mSequence = "";
				}
			}

			@Override
			public void post(node n) {
			}
		});




		for (int seqIndex_ = 0; seqIndex_ < seqLength; seqIndex_++) {
			final int seqIndex = seqIndex_;
			int[][] A = new int[indexNode.size()][options.length];

			recursiveCall(root, new TreeRecursiveUpdate() {				
				@Override
				public void pre(node n) {
				}

				@Override
				public void post(node n) {
					int ni = nodeIndex.get(n);
					if (n.numChilds == 0) {
						char c = Character.toUpperCase(n.mSequence.charAt(seqIndex));
						int cIndex = optionsIndex[(int) c];
						for (int i = 0; i<options.length; i++) {
							A[ni][i] = INFINITY;
						}
						A[ni][cIndex] = 0;
					} else {
						for (int cIndex = 0; cIndex < options.length; cIndex++) {
							int s = 0;
							for (int i = 0; i < n.numChilds; i++) {
								node child = n.childArray[i];
								int childIndex = nodeIndex.get(child);
								int sChild = INFINITY;
								for (int c2Index = 0; c2Index < options.length; c2Index++) {
									sChild = Math.min(sChild, costMatrix[cIndex][c2Index] + A[childIndex][c2Index]);
								}
								s += sChild;
							}
							A[ni][cIndex] = s;
						}
					}
				}
			});

			int[] selectedChar = new int[indexNode.size()];
			recursiveCall(root, new TreeRecursiveUpdate() {

				@Override
				public void post(node n) {
				}

				@Override
				public void pre(node n) {
					int ni = nodeIndex.get(n);
					int s = INFINITY, sI = -1;
					if (n.parent == null) {
						for (int cIndex = 0; cIndex < options.length; cIndex++) {
							if (A[ni][cIndex] < s) {
								s = A[ni][cIndex];
								sI = cIndex;
							}
						}
					} else {
						int parentSelectedChar = selectedChar[nodeIndex.get(n.parent)];
						for (int cIndex = 0; cIndex < options.length; cIndex++) {
							int cost = costMatrix[parentSelectedChar][cIndex] + A[ni][cIndex];
							if (cost < s || (cost == s && cIndex == parentSelectedChar)) {
								s = cost;
								sI = cIndex;
							}
						}
					}
					// if (s > 1000) {
					// 	throw new RuntimeException("Something is wrong!");
					// }
					// 
					// if (n.numChilds != 0 && seqIndex == 382-1) {
					// 	System.err.println("S: n=" + n.ID + n.label + " i=" + seqIndex + " s="+s + " sI=" + sI + options[sI] + " A[ni]="+Arrays.toString(A[ni]));
					// }
					if (n.numChilds != 0) {
						if (sI == -1) {
							System.err.println("No minimu found?! s=" + s + " sI=" + sI + " ID="+n.ID + " parent=" + (n.parent != null ? null : n.parent.ID) + " A[ni]="+Arrays.toString(A[ni]));
						}
						selectedChar[ni] = sI;
						n.mSequence += options[sI];
						assert(n.mSequence.length() == seqIndex+1);	
					}
				}
			});
		}

		// recursiveCall(root, new TreeRecursiveUpdate() {
		// 	@Override
		// 	public void pre(node n) {
		// 		if (n.numChilds == 0) {
		// 			System.err.println("" + n.ID + ": " + n.mSequence + "(LEAF)");
		// 		} else {
		// 			System.err.println("" + n.ID + ": " + n.mSequence);
		// 		}
		// 	}

		// 	@Override
		// 	public void post(node n) {
		// 	}
			
		// });
	}

	private void recursiveCall(node curNode, TreeRecursiveUpdate callable) {
		callable.pre(curNode);
		for (int i = 0; i < curNode.numChilds; i++) 
			recursiveCall(curNode.childArray[i], callable);
		callable.post(curNode);
	}

	// private void fillNodeIndex(node curNode, Map<node, Integer> nodeIndex, List<node> indexNode) {
	// 	for (int i = 0; i < curNode.numChilds; i++) 
	// 		fillNodeIndex(curNode.childArray[i], nodeIndex, indexNode);
	// }

	// private void fillCharCount(node curNode, boolean[] charCount) {
	// 	if (curNode.numChilds == 0) {
	// 		//check if no bad character is in the sequences
	// 		boolean[] allowedOptions = new boolean[256];
	// 		if (curNode.type.equals("DNA")) {
	// 			for (char c : "ATCG-".toCharArray()) {
	// 				allowedOptions[(int)c] = true;
	// 			}
	// 		} else {
	// 			for (char c : "ARNDCQEGHILKMFPSTWYV-".toCharArray()) {
	// 				allowedOptions[(int)c] = true;
	// 			}
	// 		}
	// 		for (char c : curNode.mSequence.toCharArray()) {
	// 			if (!allowedOptions[Character.toUpperCase(c)]) {
	// 				throw new RuntimeException("Invalid character in " + curNode.ID + " c='" + c + "' with type " + curNode.type);
	// 			}
	// 		}
	// 		//end of check

	// 		for (char c : curNode.mSequence.toCharArray()) {
	// 			charCount[(int)c] = true;
	// 		}
	// 	}
	// 	for (int i = 0; i < curNode.numChilds; i++) 
	// 		fillCharCount(curNode.childArray[i], charCount);
	// }
}
