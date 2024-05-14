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

public class HIMat {
	public double	table [][];
	public String	names [];
	public int	dim, nSera, nAg;
	public boolean	isThreshold [][],
			isSerum [],
			isAntigen [];
	private	double	ref [];
	
	public HIMat () {
		this.table = null;
		this.names = null;
		this.dim = this.nSera = this.nAg = 0;
		this.isThreshold = null;
		this.isSerum = null;
		this.isAntigen = null;
		this.ref = null;
	}
	
	private void getHIMat (String fileName) {
		try {
			BufferedReader	br;
			String	inpLine = new String (),
				tmpSplit [] = null,
				tmpMat [][] = null,
				interMat [][] = null,
				row_names [] = null,
				col_names [] = null,
				tmp_name = null;
			int	n = 0,
				row_indices [] = null,
				col_indices [] = null,
				index = 0,
				matched = 0,
				double_matched = 0;
			Stack<String>	tmpStack = new Stack<String> ();
			boolean	isMatched = false,
				multiple = false;
			double	tmp_ref[] = null;
			
			br = new BufferedReader (new FileReader (new File (fileName)));
			// read HI table
			while ((inpLine = br.readLine()) != null) {
				tmpStack.push (inpLine);
				n++;
			}
			// br.close();
			
			row_names = new String [n-2];
			tmpMat = new String [n-2][];
			
			while (!tmpStack.isEmpty()) {
				inpLine = tmpStack.pop();
				tmpSplit = inpLine.split("\t",-1);
				if (n==1) {
					col_names = new String [tmpSplit.length-1];
					System.arraycopy(tmpSplit,1,col_names,0,col_names.length);
				}
				else if (n == 2) {
					if (tmpSplit[1] == "") continue;
					tmp_ref = new double [tmpSplit.length-1];
					for (int i = 1; i < tmpSplit.length; i++) tmp_ref[i-1] = (tmpSplit[i].equals("NaN") ? Double.NaN : Double.parseDouble(tmpSplit[i]));
				}
				else {
					row_names [n-3] = tmpSplit [0];
					tmpMat [n-3] = new String [tmpSplit.length-1];
					System.arraycopy(tmpSplit,1,tmpMat [n-3],0,tmpSplit.length-1);
				}
				n--;
			}
			br.close();
			
			// reorder columns and rows
			row_indices = new int [row_names.length];
			col_indices = new int [col_names.length];
			index = 0;
			
			for (int i = 0; i < col_names.length; i++) {
				if (col_names [i].indexOf("_") != -1) {
					tmp_name = col_names [i].substring(0,col_names[i].indexOf("_"));
					multiple = true;
				}
				else {
					tmp_name = col_names [i];
					multiple = false;
				}
				
				isMatched = false;
				for (int j = 0; j < row_names.length; j++) {
					if (tmp_name.equals(row_names[j])) {
						col_indices [i] = ++index;
						if (row_indices [j] == 0 && (multiple ? col_names[i].substring(col_names [i].indexOf("_")+1,col_names [i].length()).equals("0") : true)) row_indices [j] = index;
						else double_matched++;
						matched++;
						isMatched = true;
						break;
					}
				}
				
				if (!isMatched) {
					if (col_names [i].indexOf("_") != -1) tmp_name = col_names [i].substring(0,col_names[i].indexOf("_")) + ".0";
					else tmp_name = col_names [i] + ".0";
					for (int j = 0; j < row_names.length; j++) {
						if (tmp_name.equals(row_names[j])) {
							col_indices [i] = ++index;
							if (row_indices [j] == 0) row_indices [j] = index;
							else double_matched++;
							matched++;
							isMatched = true;
							if (col_names [i].indexOf("_") == -1) col_names [i] = tmp_name;
							break;
						}
					}
				}
			}
			
			for (int i = 0; i < col_names.length; i++) {
				if (col_indices [i] == 0) col_indices [i] = ++index;
			}
			for (int i = 0; i < row_names.length; i++) {
				if (row_indices [i] == 0) row_indices [i] = ++index;
			}
			
			this.dim = row_names.length + col_names.length - matched + double_matched;
			
			interMat = new String [dim][dim];
			if (tmp_ref != null) this.ref = new double[tmp_ref.length];
			
			for (int i = 0; i < col_names.length; i++) {
				for (int j = 0; j < row_names.length; j++) {
					interMat [row_indices [j]-1][col_indices [i]-1] = tmpMat [j][i];
				}
				
				if (tmp_ref != null) {
					this.ref[col_indices[i]-1] = tmp_ref[i];
				}
			}
			
			// set names
			this.names = new String [this.dim];
			this.isSerum = new boolean [this.dim];
			this.isAntigen = new boolean [this.dim];
			index = 0;
				// columns
			for (int i = 0; i < col_names.length; i++) {
				for (int j = 0; j < col_names.length; j++) {
					if (col_indices [j] == (i+1)) {
						this.isSerum [index] = true;
						this.names [index++] = col_names [j];
						break;
					}
				}
			}
				// rows
			for (int i = col_names.length; i < this.dim; i++) {
				for (int j = 0; j < row_names.length; j++) {
					if (row_indices [j] == (i+1)) {
						this.names [index++] = row_names [j];
						break;
					}
				}
			}
			
			// set isAntigen
			index = 0;
			for (int i = 0; i < row_names.length; i++) {
				for (int j = 0; j < row_names.length; j++) {
					if (row_indices [j] == (i+1)) {
						this.isAntigen [index++] = true;
						break;
					}
				}
			}
			
			// set matrix
			this.table = new double [this.dim][this.dim];
			this.isThreshold = new boolean [this.dim][this.dim];
			
			// convert first
			for (int i = 0; i < dim; i++) {
				for (int j = 0; j < dim; j++) {
					if (interMat [i][j] == null || interMat [i][j].length () == 0 || interMat [i][j].equals ("*") || interMat [i][j].equals ("nd")) {
						this.table [i][j] = Double.NaN;
					}
					else if (interMat [i][j].charAt(0) == '<') {
						// this.table [i][j] = Double.NaN;
						this.table [i][j] = Double.parseDouble(interMat [i][j].substring(1,interMat [i][j].length()));
						this.isThreshold [i][j] = true;
					}
					else {
						this.table [i][j] = Double.parseDouble(interMat [i][j]);
					}
				}
			}
		}
		catch (IOException e) {
			System.out.println("Caught the following IOException: " + e);
		}
	}
	
	private void getDistances () {
		double	maxValues [] = null;
		
		if (this.ref == null) {
			maxValues  = new double [this.dim];
			for (int i = 0; i < this.dim; i++) maxValues [i] = 0.0;
			
			for (int i = 0; i < this.dim; i++) {
				for (int j = 0; j < this.dim; j++) {
					if (!Double.isNaN(this.table [i][j]) && this.table [i][j] > maxValues [j]) maxValues [j] = this.table [i][j];
				}
			}
		}
		else maxValues = this.ref;
		
		// normalize by maximum value
		for (int i = 0; i < this.dim; i++) {
			for (int j = 0; j < this.dim; j++) {
				if (!Double.isNaN(this.table [i][j])) {
					this.table [i][j] = maxValues [j] / this.table [i][j];
				}
			}
		}
	}
	
	private double logNormalize (double a) {
		return Math.round(Math.log(a)/Math.log(2)*10000.0)/10000.0;
	}
	
	private void log2 () {
		for (int i = 0; i < this.dim; i++) {
			for (int j = 0; j < this.dim; j++) {
				if (!Double.isNaN (this.table [i][j])) this.table [i][j] = this.logNormalize (this.table [i][j]);
			}
		}
	}
	
	public void writeHIMat (String fileName, boolean useThresholds){
		try {
			BufferedWriter	bw = new BufferedWriter (new FileWriter (new File (fileName)));
			String	newLine = System.getProperty("line.separator");
			
			// write names
			for (int i = 0; i < this.dim; i++) bw.write("\t" + this.names[i]);
			bw.write(newLine);
			
			// write table
			for (int i = 0; i < this.dim; i++) {
				bw.write(this.names[i] + "\t");
				for (int j = 0; j < this.dim; j++) {
					if (!this.isSerum [j]) break;
					if (useThresholds) bw.write ((this.isThreshold [i][j] ? "<" : "") + this.table [i][j] + "\t");
					else bw.write ("" + (this.isThreshold [i][j] ? Double.NaN : this.table [i][j]) + "\t");
				}
				bw.write (newLine);
			}
			bw.close();
		}
		catch (IOException e) {
			System.out.println("Caught the following IOException: " + e);
		}
	}
	
	private void countN () {
		// antigens
		this.nAg = this.dim;
		// antisera
		for (int i = 0; i < this.dim; i++) if (this.isSerum[i]) this.nSera++;
	}
	
	private void mergeColumns () {
		int	tmp_index = 0,
			index = 0,
			numSameSera = 0,
			tmp_indices [] = null,
			indices [] = null,
			middle_index = 0,
			numUse = this.dim;
		double	new_table_values [] = null,
			tmp_values [] = null,
			help_value = 0.0;
		boolean	new_isThreshold_values [] = null,
			tmp_isThreshold = false,
			use [] = new boolean [this.dim];
		
		for (int i = 0; i < this.dim; i++) use [i] = true;
		
		for (int i = 0; i < this.dim; i++) {
			tmp_index = this.names [i].indexOf("_");
			if (tmp_index != -1 && Integer.parseInt (this.names [i].substring(tmp_index+1,this.names [i].length())) == 0) {
				this.names [i] = this.names [i].substring(0,tmp_index);
				// collect sera
				tmp_indices = new int [this.dim - i];
				tmp_indices [0] = i;
				numSameSera = 1;
				for (int j = 0; j < this.dim; j++) {
					tmp_index = this.names [j].indexOf("_");
					if (tmp_index != -1 && this.names [j].substring(0,tmp_index).equals(this.names[i]) && Integer.parseInt (this.names [j].substring(tmp_index+1,this.names [j].length())) != 0) {
						tmp_indices [numSameSera++] = j;
					}
				}
				
				indices = new int [numSameSera];
				System.arraycopy(tmp_indices,0,indices,0,numSameSera);
				
				new_table_values = new double [this.dim];
				new_isThreshold_values = new boolean [this.dim];
				// similar to mergeDuplicates, but look at rows instead
				for (int k = 0; k < this.dim; k++) {
					tmp_index = 0;
					tmp_values = new double [indices.length];
					tmp_isThreshold = false;
					// ignore threshold values
					for (int l = 0; l < indices.length; l++) {
						if (!Double.isNaN(this.table [k][indices [l]]) && !this.isThreshold[k][indices[l]]) {
							tmp_values [tmp_index++] = this.table [k][indices[l]];
						}
					}
					// if nothing was found check for threshold values
					if (tmp_index == 0) {
						for (int l = 0; l < indices.length; l++) {
							if (!Double.isNaN(this.table [k][indices [l]])) {
								tmp_values [tmp_index++] = this.table [k][indices[l]];
							}
						}
						if (tmp_index != 0) tmp_isThreshold = true;
					}
					
					new_isThreshold_values [k] = tmp_isThreshold;
					
					if (tmp_index == 0) {
						new_table_values [k] = Double.NaN;
					}
					else if (tmp_index == 1) {
						new_table_values [k] = tmp_values [0];
					}
					else {
						// set new values
						// sort values
						for (int m = 0; m < tmp_index; m++) {
							for (int n = tmp_index-1; n > m; n--) {
								if (tmp_values [n] < tmp_values [n-1]) {
									help_value = tmp_values [n];
									tmp_values [n] = tmp_values [n-1];
									tmp_values [n-1] = help_value;
								}
							}
						}
						
						middle_index = tmp_index / 2;
						
						if ((tmp_index % 2) == 1) {
							new_table_values [k] = tmp_values [middle_index];
						}
						else {
							new_table_values [k] = Math.round((tmp_values [middle_index] + tmp_values[middle_index-1]) / 2.0 * 1000.0) / 1000.0;
						}
					}
				}
				
				for (int k = 0; k < this.dim; k++) {
					this.table [k][indices [0]] = new_table_values [k];
					this.isThreshold [k][indices [0]] = new_isThreshold_values [k];
				}
				
				for (int k = 1; k < indices.length; k++) {
					use [indices [k]] = false;
					numUse--;
				}
			}
		}
		
		
		// remove unused stuff
		double	new_table [][] = new double [numUse][numUse];
		boolean	new_isThreshold [][] = new boolean [numUse][numUse],
			new_isSerum [] = new boolean [numUse],
			new_isAntigen [] = new boolean [numUse];
		String	new_names [] = new String [numUse];
		
		index = 0;
		for (int i = 0; i < this.dim; i++) {
			if (!use [i]) continue;
			tmp_index = 0;
			for (int j = 0; j < this.dim; j++) {
				if (!use[j]) continue;
				new_table [index][tmp_index] = this.table [i][j];
				new_isThreshold [index][tmp_index] = this.isThreshold [i][j];
				tmp_index++;
			}
			new_isSerum [index] = this.isSerum [i];
			new_isAntigen [index] = this.isAntigen [i];
			new_names [index] = this.names [i];
			index++;
		}
		
		this.table = new_table;
		this.isThreshold = new_isThreshold;
		this.names = new_names;
		this.isSerum = new_isSerum;
		this.isAntigen = new_isAntigen;
		this.dim = numUse;
	}
	
	private void mergeRows () {
		int	indices [] = null,
			tmp_indices [] = null,
			index = 0,
			tmp_index = 0,
			middle_index = 0,
			numUse = this.dim,
			numSameAg = 0,
			regular_i = 0;
		double	tmp_values [] = null,
			new_table_values [] = null,
			help_value = 0.0;
		boolean	new_isThreshold_values [] = null,
			tmp_isThreshold = false,
			use [] = new boolean [this.dim];
		String	tmp_name = null;
		
		for (int i = 0; i < this.dim; i++) use [i] = true;
		
		for (int i = 0; i < this.dim; i++) {
			if (!use [i]) continue;
			
			tmp_index = this.names [i].indexOf(".");
			if (tmp_index != -1) {
				regular_i = -1;
				if (Integer.parseInt (this.names [i].substring(tmp_index+1,this.names [i].length())) != 0) {
					// check if zero exists
					tmp_name = this.names [i].substring(0,tmp_index) + ".0";
					for (int j = 0; j < this.dim; j++) {
						if (!use [i]) continue;
						
						if (this.names [j].equals(tmp_name)) {
							regular_i = i;
							tmp_name = tmp_name.substring(0,tmp_name.length()-2);
							i = j;
							break;
						}
					}
					
					if (regular_i == -1) {
						tmp_name = this.names [i].substring(0,tmp_index);
						for (int j = 0; j < this.dim; j++) {
							if (!use [i]) continue;
							
							if (this.names [j].equals(tmp_name)) {
								regular_i = i;
								i = j;
								break;
							}
						}
					}
					
					this.names [i] = tmp_name;
				}
				else {
					this.names [i] = this.names [i].substring(0,tmp_index);
				}
				
				// collect ag
				tmp_indices = new int [this.dim - i];
				tmp_indices [0] = i;
				numSameAg = 1;
				for (int j = 0; j < this.dim; j++) {
					if (!use [j]) continue;
					tmp_index = this.names [j].indexOf(".");
					if (tmp_index != -1 && this.names [j].substring(0,tmp_index).equals(this.names[i]) && Integer.parseInt (this.names [j].substring(tmp_index+1,this.names [j].length())) != 0) {
						tmp_indices [numSameAg++] = j;
					}
				}
				
				indices = new int [numSameAg];
				System.arraycopy(tmp_indices,0,indices,0,numSameAg);
				
				new_table_values = new double [this.dim];
				new_isThreshold_values = new boolean [this.dim];
				// for each column look at the values and decide what to take
				// first check if there are non threshold values?
				// if only threshold values look at these, otherwise ignore threshold values
				for (int k = 0; k < this.dim; k++) {
					tmp_index = 0;
					tmp_values = new double [indices.length];
					tmp_isThreshold = false;
					// ignore threshold values
					for (int l = 0; l < indices.length; l++) {
						if (!Double.isNaN(this.table [indices [l]][k]) && !this.isThreshold[indices[l]][k]) {
							tmp_values [tmp_index++] = this.table [indices[l]][k];
						}
					}
					// if nothing was found check for threshold values
					if (tmp_index == 0) {
						for (int l = 0; l < indices.length; l++) {
							if (!Double.isNaN(this.table [indices [l]][k])) {
								tmp_values [tmp_index++] = this.table [indices[l]][k];
							}
						}
						if (tmp_index != 0) tmp_isThreshold = true;
					}
					
					new_isThreshold_values [k] = tmp_isThreshold;
					
					if (tmp_index == 0) {
						new_table_values [k] = Double.NaN;
					}
					else if (tmp_index == 1) {
						new_table_values [k] = tmp_values [0];
					}
					else {
						// set new values
						// sort values
						for (int m = 0; m < tmp_index; m++) {
							for (int n = tmp_index-1; n > m; n--) {
								if (tmp_values [n] < tmp_values [n-1]) {
									help_value = tmp_values [n];
									tmp_values [n] = tmp_values [n-1];
									tmp_values [n-1] = help_value;
								}
							}
						}
						
						middle_index = tmp_index / 2;
						
						if ((tmp_index % 2) == 1) {
							new_table_values [k] = tmp_values [middle_index];
						}
						else {
							new_table_values [k] = Math.round((tmp_values [middle_index] + tmp_values[middle_index-1]) / 2.0 * 1000.0) / 1000.0;
						}
					}
				}
				
				for (int k = 0; k < this.dim; k++) {
					this.table [indices [0]][k] = new_table_values [k];
					this.isThreshold [indices [0]][k] = new_isThreshold_values [k];
				}
				
				for (int k = 1; k < indices.length; k++) {
					use [indices [k]] = false;
					numUse--;
				}
			
				if (regular_i != -1) {
					i = regular_i;
				}
			}
		}
		
		// remove unused stuff
		double	new_table [][] = new double [numUse][numUse];
		boolean	new_isThreshold [][] = new boolean [numUse][numUse],
			new_isSerum [] = new boolean [numUse],
			new_isAntigen [] = new boolean [numUse];
		String	new_names [] = new String [numUse];
		
		index = 0;
		for (int i = 0; i < this.dim; i++) {
			if (!use [i]) continue;
			tmp_index = 0;
			for (int j = 0; j < this.dim; j++) {
				if (!use[j]) continue;
				new_table [index][tmp_index] = this.table [i][j];
				new_isThreshold [index][tmp_index] = this.isThreshold [i][j];
				tmp_index++;
			}
			new_isSerum [index] = this.isSerum [i];
			new_isAntigen [index] = this.isAntigen [i];
			new_names [index] = this.names [i];
			index++;
		}
		
		this.table = new_table;
		this.isThreshold = new_isThreshold;
		this.names = new_names;
		this.isSerum = new_isSerum;
		this.isAntigen = new_isAntigen;
		this.dim = numUse;
	}
	
	public void getLog2Table (String infile_name, boolean isDistance) {
		//System.out.print("> read data\t\t\t");
		this.getHIMat (infile_name);
		//System.out.println("done");
		
		if (!isDistance) {
			//System.out.print("> retrieve distances\t\t");
			this.getDistances ();
			this.log2 ();
			//System.out.println("done");
		}
		
		//System.out.print("> merge sera\t\t\t");
		this.mergeColumns ();
		//System.out.println("done");
		
		//System.out.print("> merge ag\t\t\t");
		this.mergeRows ();
		//System.out.println("done");
		
		this.countN ();
	}
}
