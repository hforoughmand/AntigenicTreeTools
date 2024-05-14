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

import java.io.*;
import java.util.*;

public class treeObj {
	
	// -------------------- read data -------------------- //
	
	public String readNewickTree (String fileName, String nexus) {
		// read tree in newick format from file
		// can be in nexus format (tree block) or phylip format (only the newick tree is given)
		try {
			BufferedReader br = new BufferedReader (new FileReader (new File (fileName)));
			String	inpLine = br.readLine(),
				newickString = "";
			while (inpLine != null) {
				if (nexus != null) {
					newickString += inpLine;
					break;
				}
				else {
					int index = inpLine.indexOf("[&U]");
					index = (index == -1 ? inpLine.indexOf("[&R]") : index);
					if (index != -1) {
						index = inpLine.indexOf("(");
						newickString = inpLine.substring(index,inpLine.length());
						break;
					}
				}
				inpLine = br.readLine();
			}
			br.close();
			return newickString;
		}
		catch (IOException e) { 
			System.out.println("Caught the following exception: " + e); 
			System.exit(1);
		}
		return "";
	}
	
	public String[][] readSequences (String fileName) {
		// read sequences from file
		// sequences have to be in fasta format
		// sequences first put onto a stack, taken in reverse order and put into an array
		// each sequence has 2 fields: first contains identifier, second the sequence
		try {
			BufferedReader br = new BufferedReader (new FileReader (new File (fileName)));
			Stack<String[]> tmpData = new Stack<String[]>();
			int numSeqs = 0;
			String	inpLine = br.readLine(),
				tmpObject [] = null;
			while (inpLine != null) {
				if (inpLine.length() == 0) break;
				if (inpLine.charAt(0) == '>') {
					if (numSeqs != 0) {
						tmpData.push(tmpObject);
					}
					tmpObject = new String [2];
					tmpObject [0] = inpLine.substring(1,inpLine.length()).trim();
					tmpObject [1] = "";
					numSeqs++;
				}
				else {
					tmpObject [1] += inpLine;
				}
				inpLine = br.readLine();
			}
			tmpData.push(tmpObject);
			br.close();
			
			String resObject [][] = new String [numSeqs][];
			for (int i = numSeqs; i > 0; i--) {
				resObject [i-1] = (String[])tmpData.pop();
			}
			
			return resObject;
		}
		catch (IOException e) { 
			System.out.println("Caught the following exception: " + e); 
			System.exit(1);
		}
		return null;
	}
	
	public String [][] getSeqSplit (String seqs [][], String pattern) {
		String	resArray [][] = new String [seqs.length][],
				tmpSplit [] = null;
		
		for (int i = 0; i < seqs.length; i++) {
			tmpSplit = seqs [i][1].split(pattern);
			resArray [i] = new String [tmpSplit.length + (pattern == "" ? 0 : 1)];
			resArray [i][0] = seqs [i][0];
			System.arraycopy(tmpSplit, (pattern == "" ? 1 : 0), resArray [i], 1, tmpSplit.length - (pattern == "" ? 1 : 0));
		}
		
		return resArray;
	}
	
	public String [][] readLinkage (String fileName) {
		// read linkage list
		// layout: children	parent (tab seperated)
		// in preorder
		try {
			BufferedReader br = new BufferedReader (new FileReader (new File (fileName)));
			Stack<String[]> tmpData = new Stack<String[]>();
			int numLinks = 0;
			String	inpLine = br.readLine(),
				tmpObject [] = null;
			while (inpLine != null) {
				tmpObject = inpLine.split("\t");
				tmpData.push(tmpObject);
				numLinks++;
				inpLine = br.readLine();
			}
			br.close();
			
			String linkData[][] = new String [numLinks][];
			for (int i = numLinks; i > 0; i--) {
				linkData [i-1] = (String[])tmpData.pop();
			}
			return linkData;
		}
		catch (IOException e) { 
			System.out.println("Caught the following exception: " + e); 
			System.exit(1);
		}
		return null;
	}
	
	// -------------------- build tree -------------------- //
	
	public int createTree (String newick, node root, int counter, String type) {
		// build tree from newick format
		// internal nodes don't have labels (if do so, they are interpreted as bootstrap values)
		int index = 0;
		int colonPos, commaPos, bracePos, semicolonPos, actPos;
		
		String 	splitArray [] = newick.split(",");
		int	maxChilds = splitArray.length;
		node 	curNode = root;
		
		while (index < newick.length()) {
			if (newick.charAt(index) == '(') {	// child
				curNode = new node (counter++,curNode,maxChilds);
				curNode.identifier = "" + (counter-1);
				curNode.parent.childArray[curNode.parent.numChilds++] = curNode;
				curNode.type = type;
				index++;
			}
			else if (newick.charAt(index) == ',') {	// another child
				curNode = new node (counter++,curNode.parent,maxChilds);
				curNode.identifier = "" + (counter-1);
				curNode.parent.childArray[curNode.parent.numChilds++] = curNode;
				curNode.type = type;
				index++;
			}
			else if (newick.charAt(index) == ')') {	// closing brace
				curNode = curNode.parent;
				index++;
				// look for internal labels
				if ((newick.charAt(index) != ')') && (newick.charAt(index) != ',') && (newick.charAt(index) != ':')) {
					commaPos = newick.indexOf(',',index);
					bracePos = newick.indexOf(')',index);
					colonPos = newick.indexOf(':',index);
					actPos = Math.min(Math.min((colonPos == -1 ? newick.length() : colonPos),(commaPos== -1 ? newick.length() : commaPos)),(bracePos == -1 ? newick.length() : bracePos));
					
					curNode.identifier = newick.substring(index,actPos);
					if (curNode.identifier.indexOf(';') > -1) curNode.identifier = curNode.identifier.substring(0,curNode.identifier.indexOf(';'));
					index = actPos;
				}
			}
			else if (newick.charAt(index) == ';') {	// end of tree
				break;
			}
			else if (newick.charAt(index) == ':') {	// distance
				semicolonPos = newick.indexOf(';',index);
				commaPos = newick.indexOf(',',index);
				bracePos = newick.indexOf(')',index);
				actPos = Math.min(Math.min((semicolonPos == -1 ? newick.length() : semicolonPos),(commaPos== -1 ? newick.length() : commaPos)),(bracePos == -1 ? newick.length() : bracePos));
				curNode.distance = Double.parseDouble(newick.substring(index+1,actPos));
				if (curNode.distance < 0.0) curNode.distance = 0.0;
				index = actPos;
			}
			else {	// child values, after '(' or ','
				commaPos = newick.indexOf(',',index);
				bracePos = newick.indexOf(')',index);
				colonPos = newick.indexOf(':',index);
				actPos = Math.min(Math.min((colonPos == -1 ? newick.length() : colonPos),(commaPos== -1 ? newick.length() : commaPos)),(bracePos == -1 ? newick.length() : bracePos));//Math.min((commaPos==-1?newick.length():commaPos),(bracePos==-1?newick.length():bracePos));
				curNode.identifier = newick.substring(index,actPos);
				index = actPos;
			}
		}
		pruneChildArrays (root);
		return counter;
	}
	
	// -------------------- modify tree -------------------- //
	
	public void cleanIdentifier (node curNode) {
		if (curNode.numChilds > 0) {
			for (int i = 0; i < curNode.numChilds; i++) this.cleanIdentifier(curNode.childArray [i]);
			curNode.identifier = "";
		}
	}
	
	private void pruneChildArrays (node curNode) {
		// crop children array
		// replace array by array with only the given children nodes
		// original array is initialised with maximum number of childs (indicated by number of colons ind the newick string)
		node tmpArray [] = new node [curNode.numChilds];
		System.arraycopy(curNode.childArray,0,tmpArray,0,curNode.numChilds);
		curNode.childArray = tmpArray;
		
		for (int i = 0; i < curNode.numChilds; i++) {
			pruneChildArrays(curNode.childArray[i]);
		}
	}

	public void mapIntermediates (node root, String linkage[][]) {
		node curNode = null;
		
		// array of linkages
		// 0:child	1:parent
		// look for child, if can be found, curNode is a leave node, then map identifier to parent
		// if not found, look for parent (if can't be found, semantic error in the list), go through children list, and map identifier to first children which has no label
		// linkage list is ordered, preorder traversal
		
		for (int i = 0; i < linkage.length; i++) {
			curNode = root.getNode(linkage[i][0]);
			if (curNode != null) {
				if (curNode.parent.identifier.equals("")) curNode.parent.identifier = new String(linkage[i][1]);
			}
			else {
				curNode = root.getNode(linkage[i][1]);
				if (curNode == null) {
					curNode = root.getNode(linkage[linkage.length-1][0]);
					if (curNode == null) {
						System.out.println(linkage[i][0]);
						System.out.println(linkage[i][1]);
						System.out.println(linkage[linkage.length-1][0]);
						System.out.println("Semantic error! Node not found while trying to map intermediate identifiers!");
						System.exit(1);
					}
					else {
						if (curNode.parent.identifier.equals("")) {
							curNode.parent.identifier = new String(linkage[i][1]);
							curNode = curNode.parent;
						}
					}
				}
				for (int j = 0; j < curNode.numChilds; j++) {
					if (curNode.childArray[j].identifier.equals("")) {
						curNode.childArray[j].identifier = new String(linkage[i][0]);
						break;
					}
				}
			}
		}
	}
	
	public boolean[] mapTerminalGaps (node curNode, String seq) {
		// go through sequence from the beginning and the end
		// mark all terminal gaps by adding true for a given position
		// stop if once there is no gap
		boolean array[] = new boolean [seq.length()];
		int n = seq.length()-1;
		array[0] = (seq.charAt(0) == '-' || seq.charAt(0) == '?' ? true : false);
		array[n] = (seq.charAt(n) == '-' || seq.charAt(n) == '?' ? true : false);
		if (array[0]) {
			for (int i = 1; i <= n; i++) {
				array[i] = (seq.charAt(i) == '-' || seq.charAt(i) == '?' ? true : false);
				if (!array[i]) break;
			}
		}
		if (array[n]) {
			for (int i = n-1; i >= 0; i--) {
				array[i] = (seq.charAt(i) == '-' || seq.charAt(i) == '?' ? true : false);
				if (!array[i]) break;
			}
		}
		return array;
	}
	
	public void mapSequences (String seqs[][], node root) {
		// go through sequence array
		// add sequence to given node
		// and mark terminal gaps
		node curNode;
		for (int i = 0; i < seqs.length; i++) {
			curNode = root.getNode (seqs[i][0]);
			
			if (curNode == null) {
				System.out.println("No node for identifier '" + seqs[i][0] + "' found!");
				System.out.println("Please check input files!");
				System.exit (1);
			}
			
			curNode.mSequence = seqs [i][1];
			curNode.terminalGap = this.mapTerminalGaps (curNode, seqs[i][1]);
		}
	}
	
	public void setTerminalGaps (node curNode) {
		if (curNode.numChilds > 0) {
			for(int i = 0; i < curNode.numChilds; i++) this.setTerminalGaps(curNode.childArray [i]);
			curNode.terminalGap = this.mapTerminalGaps (curNode, curNode.mSequence);
		}
	}
	
	public void mapMutations_singleNode (node curNode) {
		// map mutations on branches, branch is from child node to parent --> mapped to child
		// root node doesn't have any mutations
		// iterate from position to position and save mutations if neither of the two sequences expose a terminal gap
		char a, b;
		int	c = 0;
		
		curNode.mutations = new StringBuffer ("");
		if (curNode.parent != null) {
			for (int i = 0; i < curNode.mSequence.length(); i++) {
				a = curNode.parent.mSequence.charAt(i);
				b = curNode.mSequence.charAt(i);
				
				if (a != b & !curNode.terminalGap[i] & !curNode.parent.terminalGap[i] & a != '?' & b != '?') {	//  & a != '-' & b != '-'
					if (b == '/' || Character.isLowerCase(b)) {
						curNode.mSequence = curNode.mSequence.substring(0, i) + a + (i == (curNode.mSequence.length()-1) ? "" : curNode.mSequence.substring(i+1));
					}
					else {
						if (curNode.type.equals("DNA") && (a == 'N' || b == 'N')) {
							curNode.mSequence = curNode.mSequence.substring(0, i) + a + (i == (curNode.mSequence.length()-1) ? "" : curNode.mSequence.substring(i+1));
							continue; // skip ambiguous characters
						}
						
						if (curNode.mutations.length() != 0) curNode.mutations.append(" ");	// if it's not the first mutation, add whitespace
						curNode.mutations.append("" + a + (i+1) + b + "");
						c++;
					}
				}
			}
		}
	}
	
	public void mapMutations (node curNode) {
		// traverse tree and map mutations to the branches
		this.mapMutations_singleNode (curNode);
		for (int i = 0; i < curNode.numChilds; i++) {
			this.mapMutations (curNode.childArray[i]);
		}
	}
	
	public void mapLabels (String mapFile, node root, int numNodes, String removed) {
		// map sequence identifiers to leave nodes
		// identical identifiers retrieve additional numbering
		String 	tmpArray [],
				tmpLabel = null,
				removedNodes[] = (removed == null ? null : removed.split(","));
		node	curNode;
		
		String	labelArray [] = new String [numNodes];
		int	numLabels[] = new int [numNodes];
		int	usedLabels = 0;
		int	index = -1;
		
		try {
			BufferedReader br = new BufferedReader (new FileReader (new File (mapFile)));
			String	inpLine = br.readLine(),
					tmpSplit [] = null;
			while (inpLine != null) {
				// parse label
				tmpArray = inpLine.split("\t");
				tmpLabel = (tmpArray[2].equals("") ? (tmpArray[1].equals("") ? tmpArray[6] : tmpArray[1]) : tmpArray[2]); // no serotype

				curNode = root.getNode(tmpArray[0]);
				if (curNode == null) {
					boolean	wasRemoved = false;
					if (removedNodes != null) {
						for (int j = 0; j < removedNodes.length; j++) if (removedNodes[j].equals(tmpArray[0])) {
							wasRemoved = true;
							break;
						}
					}
					
					if (!wasRemoved) {
						System.out.println("Couldn't map sequence name to node! No according node for identifier '" + tmpArray[0] + "' found!");
						System.out.println("Please check input!");
						System.exit(1);
					}
					
					inpLine = br.readLine();
					continue;
				}
				tmpLabel = tmpLabel.replace('\'',' ');	// remove invalid characters from name
				// search for identical names
				// go through array and look if label is already used
				// if so, increase counts for this label and add the according number to this isolate
				// otherwise save label as used
				for (int j = 0; j < usedLabels; j++) {
					if (labelArray[j].equals(tmpLabel)) {
						index = j;
						break;
					}
				}
				if (index > -1) {
					tmpLabel += " ("+(numLabels[index]+1)+")";
					numLabels[index]++;
					index = -1;
				}
				else {
					labelArray[usedLabels] = new String(tmpLabel);
					numLabels[usedLabels++]++;
				}
				curNode.label = tmpLabel;
				curNode.wholeIdentifierString = tmpArray[6];
				curNode.strain = tmpArray [2];
				curNode.year = tmpArray[4];
				curNode.serotype = tmpArray[3];
				curNode.host = tmpArray[5];
				curNode.accession = tmpArray[1];
				if (tmpArray.length > 7) {
					tmpSplit = tmpArray[7].split("\\/");
					curNode.dateYear = Integer.parseInt (tmpSplit [0]);
					curNode.dateMonth = Integer.parseInt (tmpSplit [1]);
					curNode.dateDay = Integer.parseInt (tmpSplit [2]);
				}
				inpLine = br.readLine();
			}
			br.close();
		}
		catch (IOException e) { 
			System.out.println("Caught the following exception: " + e); 
			System.exit(1);
		}
	}
	
	public void deleteSubtree (String id, node root) {
		// delete a given subtree from the whole tree
		// if this results in a node with only one child, delete this node and attach child to the parent node
		node curNode = root.getNode(id);
		if (curNode == null) {
			System.out.println("Can't delete subtree, given identifier '" + id + "' not found!");
			System.exit(1);
		}
		curNode = curNode.parent;
		if (curNode == null) {
			System.out.println ("Sorry, can't delete root node!");
			System.exit(1);
		}
		
		boolean found = false;
		
		// delete child
		for (int i = 0; i < curNode.numChilds; i++) {
			if (!found && curNode.childArray[i].identifier.equals(id)) {
				found = true;
				continue;
			}
			if (found) {
				curNode.childArray[i-1] = curNode.childArray[i];
			}
		}
		curNode.numChilds--;
		
		// check if only one child left
		if (curNode.numChilds < 2) {
			// if we are at the root node we simply put the only child as new root
			if (curNode.parent == null) {
				node	tmpChild = curNode.childArray [0];
				curNode.numChilds = curNode.childArray [0].numChilds;
				curNode.childArray = new node [curNode.numChilds];
				System.arraycopy(tmpChild.childArray,0,curNode.childArray,0,curNode.numChilds);
				//System.arraycopy(curNode.childArray [0].childArray,0,curNode.childArray,0,curNode.numChilds);
				for (int i = 0; i < curNode.numChilds; i++) curNode.childArray [i].parent = curNode;
			}
			else {
				for (int i = 0; i < curNode.parent.numChilds; i++) {
					if (curNode.parent.childArray[i].identifier.equals(curNode.identifier)) {
						// update distance
						curNode.childArray[0].distance += curNode.distance;
						// attach child to new parent
						curNode.parent.childArray[i] = curNode.childArray[0];
						// update parent pointer
						curNode.childArray[0].parent = curNode.parent;
						// update mutations on the given branch
						this.mapMutations_singleNode (curNode.childArray[0]);
						break;
					}
				}
			}
		}
	}
	
	public void pruneTree (String nodeString, node curNode) {
		// delete all given nodes from tree (and all its children)
		String nodeArray[] = nodeString.split(",");
		
		for (int i = 0; i < nodeArray.length; i++) {
			this.deleteSubtree (nodeArray[i], curNode);
		}
	}
	
	public void collapseBranches (node curNode, double threshold,node root) {
		// only collapse if not the root or a leave node
		boolean delete = false;
		int tmpNumChilds;
		if (curNode.parent != null && curNode.numChilds > 0) {
			if (curNode.distance < threshold) {
				delete = true;
				node tmpArray [] = new node [curNode.numChilds + curNode.parent.numChilds];
				// update pointer
				for (int i = 0; i < curNode.numChilds; i++) {
					curNode.childArray[i].parent = curNode.parent;
					curNode.childArray[i].distance += curNode.distance;
				}
				// copy children
				System.arraycopy (curNode.parent.childArray,0,tmpArray,0,curNode.parent.numChilds);
				System.arraycopy (curNode.childArray,0,tmpArray,curNode.parent.numChilds,curNode.numChilds);
				// update parents childArray
				curNode.parent.childArray = tmpArray;
				curNode.parent.numChilds += curNode.numChilds - 1;
				// delete curNode & update numChilds
				boolean found = false;
				for (int i = 0; i < (curNode.parent.numChilds+1); i++) {
					if (!found && curNode.parent.childArray[i].identifier.equals(curNode.identifier)) {
						found = true;
						continue;
					}
					if (found) curNode.parent.childArray[i-1] = curNode.parent.childArray[i];
				}
			}
		}
		if (!delete) {
			for (int i = 0; i < curNode.numChilds; i++) {
				tmpNumChilds = curNode.numChilds;
				this.collapseBranches (curNode.childArray[i],threshold,root);
				if (tmpNumChilds < curNode.numChilds) i--;
			}
		}
	}
	
	public void reorderTree (node curNode, boolean distance) {	// by number of nodes in the subtree
		node	tmpNode = null;
		double	tmpSum [] = null,
			tmpVal = 0;
		
		if (curNode.numChilds > 0) {
			tmpSum = new double [curNode.numChilds];
			for (int i = 0; i < curNode.numChilds; i++) {
				if (distance) {
					tmpSum [i] = curNode.childArray [i].distance;
				}
				else {
					tmpSum [i] = 0.0;
					for (int j = 0; j < curNode.childArray [i].subHosts.length; j++) {
						tmpSum [i] += curNode.childArray [i].subHosts [j];
					}
				}
			}
			
			for (int i = 0; i < (curNode.numChilds-1); i++) {
				for (int j = curNode.numChilds -1; j > i; j--) {
					if (tmpSum [j] > tmpSum [j-1]) {
						// switch values
						tmpVal = tmpSum [j-1];
						tmpSum [j-1] = tmpSum [j];
						tmpSum [j] = tmpVal;
						// switch nodes
						tmpNode = curNode.childArray [j-1];
						curNode.childArray [j-1] = curNode.childArray [j];
						curNode.childArray [j] = tmpNode;
					}
				}
			}
			
			for (int i = 0; i < curNode.numChilds; i++) this.reorderTree (curNode.childArray [i], distance);
		}
	}
	
	public void countSubHosts (node curNode) {
		if (curNode.numChilds == 0) {
			if (curNode.host.equals("human")) { curNode.subHosts [0] = 1; }
			else if (curNode.host.equals("swine")) { curNode.subHosts [1] = 1; }
			else { curNode.subHosts [2] = 1; }
		}
		else {
			for (int i = 0; i < curNode.numChilds; i++) this.countSubHosts (curNode.childArray[i]);
			for (int i = 0; i < curNode.numChilds; i++) {
				for (int j = 0; j < curNode.subHosts.length; j++) curNode.subHosts [j] += curNode.childArray[i].subHosts [j];
			}
		}
	}
	
	// -------------------- tree output -------------------- //
	
	public StringBuffer getNewick (node curNode, boolean map) {
		// return tree in newick format for given node
		// use label if a map is given, identifier otherwise
		if (curNode.numChilds == 0) {	// leave node
			StringBuffer f = new StringBuffer("");
			f.append("'" + (map ? curNode.label : curNode.identifier) + "'");
			f.append("[");
			f.append("&Mutations=\"");
			f.append(curNode.mutations.toString());
			f.append("\",internalNode=\"\"");
			if (curNode.branchIndex_up > 0) {
				f.append(",branchIndex=");
				f.append(curNode.branchIndex_up);
			}
			if (curNode.distLabel != null) {
				f.append(",distanceLabel=");
				f.append(curNode.distLabel);
			}
			f.append("]");
			if (!Double.isNaN(curNode.distance)) {
				f.append(":");
				f.append(curNode.distance);
			}
			return f;
		}
		else {	// internal node
			StringBuffer help = new StringBuffer ("");
			help.append("(");
			for (int i = 0; i < curNode.numChilds; i++) {
				help.append(getNewick(curNode.childArray[i],map));
				if (i < (curNode.numChilds-1)) {
					help.append(",");
				}
			}
			help.append(")");
			if (curNode.parent != null) {
				help.append("[&Mutations=\"");
				help.append(curNode.mutations.toString());
				help.append("\",internalNode=\"");
				help.append(curNode.identifier);
				help.append("\"");
				if (curNode.branchIndex_up > 0) {
					help.append(",branchIndex=");
					help.append(curNode.branchIndex_up);
				}
				if (curNode.distLabel != null) {
					help.append(",distanceLabel=");
					help.append(curNode.distLabel);
				}
				help.append("]");
				if (!Double.isNaN(curNode.distance)) {
					help.append(":");
					help.append(curNode.distance);
				}
			}
			//if (curNode.parent == null) help.append(";");
			return help;
		}
	}
	
	public void writeNexusTree (String newick, String fileName, String seqs [][], node root, int numNodes, boolean map) {
		// write a nexus tree file to be read by figtree
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter (new File (fileName)));
			String newLine = System.getProperty("line.separator");
			bw.write("#NEXUS");
			bw.write(newLine+newLine);
			bw.write("begin trees;"+newLine);
			bw.write("\t tree [&u] tree1 = ");
			bw.write(newick);
			bw.write(newLine+"end;"+newLine);
			bw.close();
		}
		catch (IOException e) { 
			System.out.println("Caught the following exception: " + e); 
			System.exit (1);
		}
	}
	
	public node rerootTree (node curNode, String id) {
		node	tmpRoot = curNode.getNode(id),
				tmpNode = null,
				tmpArray [] = null;
		double	tmpDistance = 0.0,
				oldDistance = 0.0,
				newDistance = 0.0;
		boolean	found = false;
		
		if (tmpRoot == null) {
			System.out.println("New root not found!");
			System.exit(1);
		}
		
		if (tmpRoot.parent == null || (tmpRoot.numChilds == 0 && tmpRoot.parent.parent == null)) {
			return curNode;
		}	// we are fine, nothing to do
		else {
			if (tmpRoot.numChilds == 0) {
				tmpRoot = tmpRoot.parent;
				tmpNode = tmpRoot.parent;
				curNode = tmpRoot;
			}
			else {
				tmpNode = tmpRoot.parent;
				curNode = tmpRoot;
			}
			
			tmpRoot.numChilds++;
			tmpArray = new node [tmpRoot.numChilds];
			System.arraycopy(tmpRoot.childArray, 0, tmpArray, 0, tmpRoot.numChilds-1);
			tmpArray [tmpRoot.numChilds-1] = tmpRoot.parent;
			tmpRoot.childArray = tmpArray;
			tmpRoot.parent = null;
			
			curNode = tmpRoot;
			tmpNode = curNode.childArray [curNode.numChilds-1];
			tmpDistance = newDistance = curNode.distance;
			
			while (tmpNode.parent != null) {
				for (int i = 0; i < tmpNode.numChilds; i++) {
					if (tmpNode.childArray [i] == curNode) {
						oldDistance = tmpNode.distance;
						tmpNode.distance = newDistance;
						newDistance = oldDistance;
						tmpNode.childArray [i] = tmpNode.parent;
						tmpNode.parent = curNode;
						curNode = tmpNode;
						tmpNode = tmpNode.childArray [i];
						break;
					}
				}
			}
			
			tmpNode.parent = curNode;
			
			tmpArray = new node [tmpNode.numChilds - 1];
			found = false;
			for (int i = 0; i < tmpNode.numChilds; i++) {
				if (tmpNode.childArray [i] == curNode) {
					found = true;
					continue;
				}
				if (found) tmpArray [i-1] = tmpNode.childArray [i];
				else tmpArray [i] = tmpNode.childArray [i];
			}
			tmpNode.childArray = tmpArray;
			tmpNode.numChilds--;
			tmpNode.distance = tmpDistance;
			
			// only one child left --> remove tmpNode
			if (tmpNode.numChilds == 1) {
				for (int i = 0; i < tmpNode.parent.numChilds; i++) {
					if (tmpNode.parent.childArray[i].identifier.equals(tmpNode.identifier)) {
						tmpNode.parent.childArray[i] = tmpNode.childArray[0];
						tmpNode.childArray[0].distance += tmpNode.distance;
						tmpNode.childArray[0].parent = tmpNode.parent;
					}
				}
			}
		}
		
		tmpRoot.distance = 0.0;
		return tmpRoot;
	}
	
	// -------------------- additional functions -------------------- //
	
	public String sequenceType (String inpSeq) {
		int	baseCount [] = new int [27],
			sum = 0,
			num = 0;
		String	seq = inpSeq.toUpperCase().replace('/', '-');
		
		for (int i = 0; i < seq.length(); i++) {
			if (seq.charAt(i) == '-') {
				baseCount[(baseCount.length-1)]++;
			}
			else {
				baseCount[((int)seq.charAt(i)) -65]++;
				num++;
			}
		}
		
		sum = baseCount [(int)'A' -65] + baseCount [(int)'C' -65] + baseCount [(int)'G' -65] + baseCount [(int)'T' -65];
		if ((((float)sum)/((float)num)) > 0.5) return "DNA";
		else return "PROTEIN";
	}
	
}
