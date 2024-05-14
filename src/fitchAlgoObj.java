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

public class fitchAlgoObj {
	
	private void fitchPhaseOne (node curNode, int pos, boolean gaps) {
		if (curNode.numChilds == 0) {
			curNode.intermAnc = new HashSet<String>();
			curNode.fitchType = "child";
			
			if (curNode.type.equals("DNA")) {
					if (curNode.tmpSequence.charAt(pos) != '-') {
						switch (curNode.tmpSequence.charAt(pos)) {
							case 'R': case 'r': curNode.intermAnc.add("A"); curNode.intermAnc.add("G"); break;
							case 'Y': case 'y': curNode.intermAnc.add("C"); curNode.intermAnc.add("T"); break;
							case 'M': case 'm': curNode.intermAnc.add("A"); curNode.intermAnc.add("C"); break;
							case 'K': case 'k': curNode.intermAnc.add("G"); curNode.intermAnc.add("T"); break;
							case 'S': case 's': curNode.intermAnc.add("C"); curNode.intermAnc.add("G"); break;
							case 'W': case 'w': curNode.intermAnc.add("A"); curNode.intermAnc.add("T"); break;
							case 'H': case 'h': curNode.intermAnc.add("A"); curNode.intermAnc.add("C"); curNode.intermAnc.add("T"); break;
							case 'B': case 'b': curNode.intermAnc.add("C"); curNode.intermAnc.add("G"); curNode.intermAnc.add("T"); break;
							case 'V': case 'v': curNode.intermAnc.add("A"); curNode.intermAnc.add("C"); curNode.intermAnc.add("G"); break;
							case 'D': case 'd': curNode.intermAnc.add("A"); curNode.intermAnc.add("G"); curNode.intermAnc.add("T"); break;
							case 'N': case 'n': case '?': curNode.intermAnc.add("A"); curNode.intermAnc.add("C"); curNode.intermAnc.add("G"); curNode.intermAnc.add("T"); if (gaps) curNode.intermAnc.add("-"); break; // if (gaps) curNode.intermAnc.add("-");
							default: curNode.intermAnc.add(""+Character.toUpperCase(curNode.tmpSequence.charAt(pos))); break;
						}
					}
					else {
						if (gaps) {
							curNode.intermAnc.add(""+Character.toUpperCase(curNode.tmpSequence.charAt(pos)));
						}
						else {
							curNode.intermAnc.add("A"); 
							curNode.intermAnc.add("C"); 
							curNode.intermAnc.add("G"); 
							curNode.intermAnc.add("T");
						}
					}
				}
				else {
					if (curNode.tmpSequence.charAt(pos) != '-') {
						switch (curNode.tmpSequence.charAt(pos)) {
							case 'B': case 'b': curNode.intermAnc.add("D"); curNode.intermAnc.add("N"); break;
							case 'Z': case 'z': curNode.intermAnc.add("E"); curNode.intermAnc.add("Q"); break;
							case 'J': case 'j': curNode.intermAnc.add("L"); curNode.intermAnc.add("I"); break;
							case 'X': case 'x': case '?': for (int i = 65; i < 91; i++) { curNode.intermAnc.add("" + (char)i); } break; //if (gaps) curNode.intermAnc.add("-"); 
							default: curNode.intermAnc.add(""+Character.toUpperCase(curNode.tmpSequence.charAt(pos))); break;
						}
					}
					else {
						if (gaps) {
							curNode.intermAnc.add(""+Character.toUpperCase(curNode.tmpSequence.charAt(pos)));
						}
						else {
							for (int i = 65; i < 91; i++) { curNode.intermAnc.add("" + (char)i); }
						}
					}
				}
		}
		else {
			for (int i = 0; i < curNode.numChilds; i++) this.fitchPhaseOne (curNode.childArray [i], pos, gaps);
			if (curNode.numChilds == 2) {
				curNode.intermAnc = new HashSet<String>(curNode.childArray[0].intermAnc);
				// if one is empty build the union
				if (curNode.childArray[0].intermAnc.isEmpty() || curNode.childArray[1].intermAnc.isEmpty()) {
					curNode.intermAnc.addAll (curNode.childArray[1].intermAnc);
					curNode.fitchType = (curNode.childArray[0].intermAnc.isEmpty() ? curNode.childArray[0].fitchType : curNode.childArray[1].fitchType);
				}
				else {
					curNode.intermAnc.retainAll(curNode.childArray[1].intermAnc);	// intersection
					if (curNode.intermAnc.isEmpty()) {
						curNode.intermAnc.addAll(curNode.childArray[0].intermAnc);	// union
						curNode.intermAnc.addAll(curNode.childArray[1].intermAnc);	// union
						curNode.fitchType = "union";
					}
					else {
						curNode.fitchType = "intersection";
					}
				}
			}
			else {	// multifurcating tree => majority voting
				curNode.fitchType = "multi";
				curNode.intermAnc = new HashSet<String>();
				
				int	intermCount [] = new int [27];
				int	max = 1,
					numSet = 0;
				
				for (int i = 0; i < curNode.numChilds; i++) {
					if (!curNode.childArray [i].intermAnc.isEmpty()) {
						for (String s : curNode.childArray[i].intermAnc) {
							if (s.equals("-")) {
									if(gaps) intermCount[intermCount.length-1]++;
							}
							else if (s.equals("?")) {
								
							}
							else {
								intermCount[((int)s.charAt(0)) -65]++;
							}
						}
					}
				}
				
				for(int i = 0; i < intermCount.length; i++) {
					if (intermCount[i] > 0) numSet++;
					if (intermCount [i] > max) {
						curNode.intermAnc = new HashSet<String>();
						curNode.intermAnc.add("" + (char) (i + (i==(intermCount.length-1) ? 19 : 65)));
						max = intermCount [i];
					}
					else if (intermCount [i] == max) curNode.intermAnc.add("" + (char) (i + (i==(intermCount.length-1) ? 19 : 65)));
				}
			}
		}
	}
	
	private Set<String> [] checkDownwards (node curNode) {
		Set<String>	returnArray [] = new Set [2];
		returnArray [0] = new HashSet<String>();
		returnArray [1] = null;
		
		if (curNode.numChilds == 0) {
			returnArray [0].add("child");
		}
		else {
			boolean	isNonEmptyChild [] = new boolean [curNode.numChilds];
			int	sum = 0;
			
			for (int i = 0; i < curNode.numChilds; i++) {
				if (!curNode.childArray[i].intermAnc.isEmpty()) {
					sum++;
					isNonEmptyChild [i] = true;
				}
			}
			
			if (sum == 1) {
				for (int i = 0; i < isNonEmptyChild.length; i++) {
					if (isNonEmptyChild [i]) {
						returnArray = this.checkDownwards(curNode.childArray[i]);
						break;
					}
				}
			}
			else if (sum == 2) {
				for (int i = 0; i < isNonEmptyChild.length; i++) {
					if (isNonEmptyChild[i]) {
						if (returnArray [1] == null) returnArray [1] = new HashSet<String>(curNode.childArray[i].intermAnc);
						else {
							Set<String> tmpSet = new HashSet<String>(returnArray[1]);
							tmpSet.retainAll(curNode.childArray[i].intermAnc);
							if (tmpSet.isEmpty()) returnArray [0].add("union");
							else returnArray [0].add("intersection");
							returnArray [1].addAll(curNode.childArray[i].intermAnc);
							break;
						}
					}
				}
			}
			else if (sum > 2) {
				returnArray[0].add("multi");
			}
		}
		
		return returnArray;
	}
	
	private void fitchPhaseTwo (node curNode, int pos, boolean gaps) {
		if (curNode.numChilds > 0) {
			if (curNode.parent == null) for (int i = 0; i < curNode.numChilds; i++) fitchPhaseTwo (curNode.childArray[i], pos, gaps);
			else if (!curNode.intermAnc.isEmpty()) {
				Set<String>	finalSet = new HashSet<String>();
				for (String sc : curNode.intermAnc) {
					for (String sp : curNode.parent.intermAnc) {
						if (sc.equals(this.toUpper(sp))) {
							finalSet.add(sc);
							break;
						}
					}
				}
				
				if (curNode.numChilds > 2) {
					if (!finalSet.isEmpty()) curNode.intermAnc = new HashSet<String>(finalSet);
				}
				else {
					boolean	fixed = true;
					for (String s : curNode.parent.intermAnc) {
						if (!finalSet.contains(s) && !finalSet.contains(this.toUpper(s))) {
							fixed = false;
							break;
						}
					}
					if (fixed) {
						curNode.intermAnc = new HashSet<String>(finalSet);
					}
					else {
						if (curNode.fitchType.equals("union")) {
							for (String s : curNode.parent.intermAnc) {
								if (!curNode.intermAnc.contains(this.toUpper(s))) {
									curNode.intermAnc.add(this.toLower(s));
								}
							}
						}
						else if (curNode.fitchType.equals("child")) {
							;
						}
						else {
							Set<String>	down [] = checkDownwards (curNode); // introduced
							if (!down[0].contains("child") && !down[0].contains("multi")) {
								for (String s : curNode.parent.intermAnc) {
									if (!curNode.intermAnc.contains(this.toUpper(s)) && down[1].contains(this.toUpper(s))) curNode.intermAnc.add(this.toLower(s));
								}
							}
						}
					}
				}
				
				for (int i = 0; i < curNode.numChilds; i++) fitchPhaseTwo (curNode.childArray[i], pos, gaps);
			}
		}
	}
	
	private void reconstructAncestralStates (node curNode, String strategy, int pos, boolean gaps, Random r) {
		String	finalSet [],
			parentChar;
		
		int	ind = 0;
		
		Set<String> alphabetSet = new HashSet<String> ();
		for (int i = 65; i < 91; i++) {
			alphabetSet.add(("" + (char)i).toUpperCase());
		}
		if (gaps) alphabetSet.add("-");
		
		if (curNode.parent == null) {
			if (curNode.intermAnc.isEmpty()) curNode.tmpSequence += "-";
			else {
				finalSet = new String [curNode.intermAnc.size()];
				for (String s : curNode.intermAnc) finalSet[ind++] = s;
				
				if (finalSet.length == 1) {
					curNode.tmpSequence += this.toUpper(finalSet [0]);
				}
				else {
					// choose randomly
					curNode.tmpSequence += this.toUpper(finalSet [Math.abs(r.nextInt()) % finalSet.length]);
				}
			}
		}
		else if (curNode.numChilds > 0){
			parentChar = this.toUpper("" + curNode.parent.tmpSequence.charAt(curNode.tmpSequence.length()));
			if (curNode.intermAnc.isEmpty()) curNode.tmpSequence += parentChar;
			else {
				finalSet = new String [curNode.intermAnc.size()];
				for (String s : curNode.intermAnc) finalSet[ind++] = s;
				if (finalSet.length == 1) {
					curNode.tmpSequence += this.toUpper(finalSet [0]);
				}
				else {
					if(curNode.intermAnc.contains(parentChar) || curNode.intermAnc.contains(this.toLower(parentChar))) {
						if (curNode.intermAnc.contains(parentChar)) curNode.tmpSequence += parentChar;
						else {
							if (strategy.equals("DelTran")) {
								curNode.tmpSequence += parentChar;
							}
							else {
								curNode.intermAnc.retainAll(alphabetSet);
								ind = 0;
								finalSet = new String [curNode.intermAnc.size()];
								for (String s : curNode.intermAnc) finalSet[ind++] = s;
								
								// choose randomly
								curNode.tmpSequence += this.toUpper(finalSet [Math.abs(r.nextInt()) % finalSet.length]);
							}
						}
					}
					else {
						ind = 0;
						curNode.intermAnc.retainAll(alphabetSet);
						finalSet = new String [curNode.intermAnc.size()];
						for (String s : curNode.intermAnc) {
							finalSet[ind++] = s;
						}
						
						// choose randomly
						curNode.tmpSequence += this.toUpper(finalSet [Math.abs(r.nextInt()) % finalSet.length]);
					}
				}
			}
		}
		
		if (curNode.numChilds > 0) for (int i = 0; i < curNode.numChilds; i++) this.reconstructAncestralStates(curNode.childArray[i], strategy, pos, gaps, r);

		/*
		Ni -> Ni || (Ni) -> Ni is obligatory
		(Ni) -> (Ni) || Ni -> (Ni) is deltran
		(Ni) -> Nj || Ni -> Nj is acctran
		(Ni) -> (Nj) || Ni -> (Nj) is not permitted
		*/
	}
	
	private String toUpper (String s) {
		if (s.charAt(0) == '^') return "-";
		else return s.toUpperCase();
	}
	
	private String toLower (String s) {
		if (s.charAt(0) == '-') return "^";
		else return s.toLowerCase();
	}

	private void assignFixedAncestrals (node curNode, int pos, String ancChar) {
		if (curNode.numChilds != 0) {
			for (int i = 0; i < curNode.numChilds; i++) this.assignFixedAncestrals (curNode.childArray [i], pos, ancChar);
			curNode.tmpSequence += ancChar;
		}
	}
	
	private void copyLeaveSeqs (node curNode) {
		if (curNode.numChilds > 0) {
			for (int i = 0; i < curNode.numChilds; i++) this.copyLeaveSeqs(curNode.childArray[i]);
			curNode.tmpSequence = "";
		}
		else curNode.tmpSequence = curNode.mSequence;
	}
	
	private int [] checkSiteRecursion (node curNode, int pos, int counts [], boolean gaps) {
		int	tmpCounts [] = new int [counts.length];
		System.arraycopy(counts,0,tmpCounts,0,counts.length);
		if (curNode.numChilds > 0) {
			for (int i = 0; i < curNode.numChilds; i++) tmpCounts = checkSiteRecursion (curNode.childArray[i], pos, tmpCounts, gaps);
		}
		else {
			try{
				if (curNode.tmpSequence.charAt(pos) == '-' || curNode.tmpSequence.charAt(pos) == '/') {
					if (gaps) tmpCounts [tmpCounts.length-2]++;
				}
				else {
					if (curNode.tmpSequence.charAt(pos) != '?') tmpCounts [((int)Character.toUpperCase(curNode.tmpSequence.charAt(pos))) - 65]++;
				}
				tmpCounts [tmpCounts.length-1]++;
			}
			catch (Exception e) {
				System.out.println(e);
				System.out.println(curNode.identifier + " " + curNode.tmpSequence + " " + pos);
				System.exit(1);
			}
		}
		
		return tmpCounts;
	}
	
	private String checkSite (node curNode, int pos, boolean gaps) {
		int	counts [] = new int [28],
			sum = 0,
			which = 0;
		for (int i = 0; i < curNode.numChilds; i++) counts = checkSiteRecursion (curNode.childArray[i], pos, counts, gaps);
		for (int i = 0; i < (counts.length-1); i++) {
			if (counts [i] > 0) {
				sum++;
				which = i;
			}
		}
		if (sum > 1) {
			return "multi";
		}
		else {
			return "" + (char)(which+(which == 26 ? 19: 65));
		}
	}
	
	private void setIntermediates (node curNode, int pos) {
		if (curNode.numChilds > 0) {
			curNode.mSequence = (pos == -1 ? new String(curNode.tmpSequence) : curNode.mSequence.subSequence(0, pos) + curNode.tmpSequence + (pos < (curNode.mSequence.length()-1) ? curNode.mSequence.substring(pos+1) : ""));
			curNode.tmpSequence = "";
			curNode.intermAnc = null;
		}
		
		for (int i = 0; i < curNode.numChilds; i++) this.setIntermediates(curNode.childArray [i], pos);
	}
	
	public void ancestralStateReconstruction (node curNode, int seqLength, boolean gaps, String strategy, Long seed) {
		String	site = "";
		
		this.copyLeaveSeqs(curNode);
		Random random = seed == null ? new Random() : new Random(seed);
		// System.err.println("Seed = " + seed);
		
		for (int i = 0; i < seqLength; i++) {
			site = this.checkSite (curNode, i, gaps);
			if (!site.equals("multi")) {
				this.assignFixedAncestrals (curNode, i,site);
			}
			else {
				this.fitchPhaseOne(curNode, i, gaps);
				this.fitchPhaseTwo(curNode, i, gaps);
				this.reconstructAncestralStates(curNode, strategy, i, gaps, random);
			}
		}
		this.setIntermediates (curNode, -1);
	}
}
