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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Stream;

import javax.management.RuntimeErrorException;


interface TreeRecursiveUpdate {

	void pre(node n);
	void post(node n);
}

public class Sankoff {
	
	/**
	 * Fills all node tmpSequence parameter based on tmpSequence of the leaf nodes
	 * @param root
	 * @param seqLength
	 * @param gaps allow gap. If false, we consider gap as unknown.
	 */
	public void ancestralStateReconstruction(node root, int seqLength, boolean gaps, String costMatrixFile) {
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
		if (costMatrixFile == null) {
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
		} else {
			for (int[] row: costMatrix)
    			Arrays.fill(row, INFINITY);
			try (Stream<String> lines = Files.lines(Paths.get(costMatrixFile))) {
				// stream.forEach(System.out::println);
				int index = 0;
				String[] firstRowOptions;
				int[] firstRowIndex = null;
				for (String line : lines.toArray(String[]::new)) {
					if (index == 0) {
						firstRowOptions = line.strip().replaceAll("  *", " ").split(" ");
						firstRowIndex = new int[firstRowOptions.length];
						for (int i=0; i<firstRowOptions.length; i++) {
							firstRowIndex[i] = optionsIndex[firstRowOptions[i].charAt(0)];
						}
					} else {
						String[] x = line.strip().replaceAll("  *", " ").split(" ");
						int rowIndex = optionsIndex[x[0].charAt(0)];
						// System.err.println("R " + rowIndex + " " + x[0].charAt(0));
						for (int j=1; j<x.length; j++) {
							if (rowIndex == -1 || firstRowIndex[j-1] == -1)
								continue;
							costMatrix[rowIndex][firstRowIndex[j-1]] = Integer.parseInt(x[j]);
						}
					}
					index++;
				}
			} catch (IOException e) {
				throw new RuntimeException(e);
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
