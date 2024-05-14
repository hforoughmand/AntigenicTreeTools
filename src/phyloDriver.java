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

public class phyloDriver {
	
	public static void printHelp () {
		String helpString = 	"Use the following options:	(options indicated with [] are optional)" +
					"\n" +
					"TREE INPUT\n" +
					"[-f strategy	-- infer intermediate sequences (AccTran/DelTran/Sankoff)]\n" +
					"[-g			-- count gaps as changes (when ancestral states are reconstructed)]\n" +
					"[-i file		-- input file with intermediate sequences in fasta format]\n" +
					"[-l file		-- file with node linkage]\n" +
					"[-m file		-- file with leaf node mapping]\n" +
					" -n file		-- file with tree in newick format\n" +
					"[-o name		-- output name]\n" +
					"[-p			-- given tree is in phylip format (default is nexus)]\n" +
					" -t file		-- input file with leave sequences in fasta format\n" +
					"\n" +
					"TREE MANIPULATION\n" +
					"[-col			-- permit branch collapsing]\n" +
					"[-not list		-- comma separated list of nodes to be pruned]\n" +
					"[-r name		-- reroot tree at leaf 'name']\n" +
					"\n" +
					"NNLS FIT\n" +
					" -ls file		-- input matrix for least squares fit\n" +
					"[-d			-- HI input matrix contains already log2 normalized distance values]\n" +
					"[-loo			-- do loo for fit?]\n" + 
					"[-cv \"x,y\"	-- do x-fold cross validation y times for fit?]\n" +
					"[-seed seed	-- seed for ancestoral reconstruction]\n";
		System.out.println(helpString);
	}
	
	public static boolean isSet (String arg) {
		return (arg == null ? false : true);
	}

	final static int ARGUMENT_SEED_INDEX = 17;
	final static int ARGUMENT_COST_MATRIX_FILE = 18;
	
	public static String [] readArgs (String [] args) {
		if (args.length == 0) {
			printHelp ();
			System.exit (1);
		}
		String arguments [] = new String [19];
		for (int i = 0; i < arguments.length; i++) arguments [i] = null;
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h")) { 
				printHelp ();
				System.exit (1);
			}
			else if (args[i].equals("-n")) { arguments [0] = args[++i]; }	//	newick tree
			else if (args[i].equals("-t")) { arguments [1] = args[++i]; }	//	leave sequences
			else if (args[i].equals("-i")) { arguments [2] = args[++i]; }	//	intermediate sequences
			else if (args[i].equals("-l")) { arguments [3] = args[++i]; }	//	linkage
			else if (args[i].equals("-p")) { arguments [4] = "p"; }		//	trees in nexus format or phylip (nexus is default)?
			else if (args[i].equals("-m")) { arguments [5] = args[++i]; }	//	leave labels?
			else if (args[i].equals("-not")) { arguments [6] = args[++i]; }	//	don't use subtrees?
			else if (args[i].equals("-o")) { arguments [7] = args[++i]; }	//	output name
			else if (args[i].equals("-col")) { arguments [8] = "col"; }	//	permit collapse branches
			else if (args[i].equals("-f")) { arguments [9] = args[++i]; }	//	infer main intermediates
			else if (args[i].equals("-r")) { arguments [10] = args[++i]; }	//	re-root tree
			else if (args[i].equals("-ls")) { arguments [11] = args[++i]; }	//	solve least squares problem
			else if (args[i].equals("-g")) { arguments [12] = "gaps"; }	//	count gaps as change?
			else if (args[i].equals("-loo")) { arguments [14] = "loo"; }	//	do leave-one-out cross validation?
			else if (args[i].equals("-cv")) { arguments [15] = args[++i]; }	//	do cross validation?
			else if (args[i].equals("-d")) { arguments [16] = "isDistance"; }	//	HI values already distances?
			else if (args[i].equals("-seed")) { arguments [ARGUMENT_SEED_INDEX] = args[++i]; }	//	HI values already distances?
			else if (args[i].equals("-cost")) { arguments [ARGUMENT_COST_MATRIX_FILE] = args[++i]; }	//	HI values already distances?
			else { 
				// System.out.println ("Unknown option '" + args[i] + "'! Use '-h' for help."); 
				throw new RuntimeException("Unknown option '" + args[i] + "'! Use '-h' for help.");
			}
		}
		if (!isSet(arguments [0])) {
			System.out.println("No tree file given!");
			System.exit(1);
		}
		if (!isSet(arguments [1])) {
			System.out.println("No sequences for leave nodes given!");
			System.exit(1);
		}
		if (!isSet(arguments [2]) && arguments [9] == null) {
			System.out.println("No sequences for intermediate nodes given!");
			System.exit(1);
		}
		if (!isSet(arguments [3]) && arguments [9] == null) {
			System.out.println("Node linkage for sequence mapping not given!");
			System.exit(1);
		}
		if (!isSet(arguments [7])) {
			System.out.println("No output name given!");
			System.exit(1);
		}
		if (!isSet(arguments [11])) {
			System.out.println("Distance matrix is missing!");
			System.exit(1);
		}
		if (isSet(arguments [14]) && isSet(arguments [15])) {
			System.out.println("Leve-one-out and cross validation specified! Use only one of the two options.");
			System.exit(1);
		}
		
		return arguments;
	}
	
	public static void main(String[] args) {
		System.out.println("\nAntigenic Tree Tools");
		System.out.println("\nPlease cite: L. Steinbrueck and A. C. McHardy (2012). Inference of Genotype-Phenotype Relationships in the Antigenic Evolution of Human Influenza A (H3N2) Viruses. PLoS Comp Biol 8(4): e1002492.\n");
		
		// read arguments
		String argv [] = readArgs(args);
		
		treeObj	tree = new treeObj();
		fitchAlgoObj	fitch = new fitchAlgoObj ();
		
		// read tree and sequences
		System.out.print ("> read tree and sequences ...\t");
		String	newick = tree.readNewickTree(argv[0],argv[4]),
			leaveSeqs[][] = null,
			interSeqs[][] = null;
		
		// 
		leaveSeqs = tree.readSequences(argv[1]);
		if (!isSet (argv [9]))	{
			interSeqs = tree.readSequences(argv[2]);
		}
		String	type = tree.sequenceType(leaveSeqs[0][1]);
		
		System.out.println ("\tdone");

		node	root = new node(0,newick.split(",").length);
		root.identifier = "0";
		
		// build tree from newick format
		System.out.print ("> build tree ...\t\t");
		int numNodes = tree.createTree(newick,root,1,type);
		if (isSet(argv [10])) root = tree.rerootTree(root,argv [10]);
		System.out.println ("\tdone");
		
		// read linkage
		if(!isSet (argv [9])) {
			System.out.print ("> map linkage ...\t\t");
			String	linkage [][] = tree.readLinkage(argv [3]);
			// map linkage to tree
			tree.cleanIdentifier(root);
			tree.mapIntermediates(root,linkage);
			System.out.println ("\tdone");
		}	
		
		// map sequences to each node, leave nodes and intermediates separated
		System.out.print ("> map sequences to nodes ...\t");
		tree.mapSequences(leaveSeqs, root);
		if (interSeqs != null && !isSet(argv[9])) tree.mapSequences(interSeqs, root);
		System.out.println ("\tdone");
		
		// collapse branches
		if (isSet(argv [8])) {
			System.out.print ("> collapse branches ...\t\t");
			tree.collapseBranches (root, 0.0000001,root);
			System.out.println ("\tdone");
		}
		
		if (isSet(argv[9])){
			System.out.print ("> infer ancestral character states ...");
			if (argv[9].equals("Sankoff")) {
				Sankoff sankoff = new Sankoff();
				sankoff.ancestralStateReconstruction(root, leaveSeqs[0][1].length(), isSet(argv[12]), argv[ARGUMENT_COST_MATRIX_FILE]);
			} else {
				fitch.ancestralStateReconstruction(root, leaveSeqs[0][1].length(), isSet(argv[12]), argv[9], isSet(argv[ARGUMENT_SEED_INDEX]) ? Long.parseLong(argv[ARGUMENT_SEED_INDEX]) : null);
			}
			tree.setTerminalGaps(root);
			System.out.println ("\tdone");
		}
		
		// map mutations to branches
		System.out.print ("> map mutations ...\t\t");
		tree.mapMutations(root);
		System.out.println ("\tdone");
		
		boolean	map = (isSet(argv [5]));

		// map sequence labels to leave nodes?
		if (isSet(argv [5])) {
			System.out.print ("> map leaf labels ...\t\t");
			tree.mapLabels(argv [5],root,numNodes, argv[6]);
			System.out.println ("\tdone");
		}
		tree.countSubHosts (root);
		
		// delete subtrees from whole tree?
		if (isSet(argv [6])) {
			System.out.print ("> delete subtree(s) ...\t\t");
			tree.pruneTree (argv [6], root);
			/*if (isSet(argv [5]))*/ tree.countSubHosts (root);
			System.out.println ("\tdone");
		}
		
		// reorder tree
		tree.reorderTree (root,false);
		
		// read HI matrix
		System.out.print ("> compute HI distances ...\t");
		HIMat	him = new HIMat();
		him.getLog2Table(argv[11], isSet(argv[16]));
		System.out.println ("\tdone");
		
		// solve least squares problem
		System.out.print ("> solve least squares problem ...\t");
		NNLSsolver	nnls = new NNLSsolver ();
		nnls.computeBranchWeights(root, him, argv [7]+(isSet(argv [6]) ? ".pruned" : "") + ".leastSquares", argv [14], argv [15], map, leaveSeqs, numNodes);
		System.out.println ("done");
		
		System.out.println("");
	}
}
