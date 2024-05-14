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

public class node {
	public	String	identifier,
			wholeIdentifierString,
			mSequence,
			tmpSequence,
			label,
			strain,
			year,
			date,
			serotype,
			host,
			accession,
			type,
			fitchType,
			distLabel;
	public	StringBuffer	mutations;
	public	node	parent,
			childArray[],
			subtreeLeafs [];
	public	int	numChilds,
			ID,
			subHosts [],
			dateYear,
			dateMonth,
			dateDay,
			branchIndex_up,
			branchIndex_down;
	public	double	distance;
	public	boolean	terminalGap [],
			negative;
	public	Set<String>	intermAnc;
	
	public node (int id) {
		this.accession = "";
		this.branchIndex_down = -1;
		this.branchIndex_up = -1;
		this.childArray = null;
		this.date = "";
		this.dateDay = 0;
		this.dateMonth = 0;
		this.dateYear = 0;
		this.distance = Double.NaN;
		this.distLabel = null;
		this.fitchType = "";
		this.host = "";
		this.ID = id;
		this.identifier = "";
		this.intermAnc = null;
		this.label = null;
		this.mSequence = "";
		this.mutations = new StringBuffer("");
		this.negative = false;
		this.numChilds = 0;
		this.parent = null;
		this.serotype = "";
		this.strain = "";
		this.subHosts = new int [3]; // human, swine, avian (everything else)
		this.subtreeLeafs = null;
		this.terminalGap = null;
		this.tmpSequence = "";
		this.type = "";
		this.wholeIdentifierString = "";
		this.year = "";
	}
	
	public node (int id, node newParent) {
		this(id);
		this.parent = newParent;
	}
	
	public node (int id, int limitChilds) {
		this(id);
		this.childArray = new node [limitChilds];
	}
	
	public node (int id, node newParent, int limitChilds) {
		this(id,newParent);
		this.childArray = new node [limitChilds];
	}
	
	/*public String toString () {
		String nodeString = new String("");
		
		nodeString += this.ID;
		nodeString += "\t" + this.identifier;
		nodeString += "\t" + this.accession;
		nodeString += "\t" + this.label;
		nodeString += "\t" + this.strain;
		nodeString += "\t" + this.year;
		nodeString += "\t" + this.serotype;
		nodeString += "\t" + this.host;
		nodeString += "\t" + this.bootstrapValue;
		nodeString += "\t" + this.mutations.toString();
		nodeString += "\t" + this.mutations_plus.toString();
		nodeString += "\t" + this.distance;
		for (int i = 0; i < this.subHosts.length; i++) nodeString += "\t" + this.subHosts[i];
		nodeString += "\t" + this.numChilds;
		nodeString += "\t" + this.mSequence;
		nodeString += "\t" + this.cSequence;
		
		return nodeString;
	}*/
	
	public node getNode (String id) {
		// is current node the right one, return this node
		// otherwise go through children array and look for the requested node
		// if it can't be found return null
		if (this.identifier.equals(id)) {
			return this;
		}
		else {
			node tmpNode = null;
			for (int i = 0; i < this.numChilds; i++) {
				tmpNode = this.childArray[i].getNode(id);
				if (tmpNode != null) return tmpNode;
			}
			return tmpNode;
		}
	}
	
}
