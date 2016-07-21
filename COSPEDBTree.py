#!/usr/bin/env python


#---------------------------------------------
""" 
this program is used to generate a supertree (consensus) from a set of constituent trees
the input is multiple source trees
each of the trees need to be decomposed to couplets 
then the couplets need to be joined
there may be conflicts among the input tree - we have to select the consensus

Author: Sourya Bhattacharyya
Dept of CSE, IIT Kharagpur
V1.0 - 15.01.2014 - public release
V1.1 - 15.03.2014 - rewritten for lowering time complexity
V1.2 - 28.02.2015 - added binary suprtree and level based scoring options
V1.3 - 30.03.2015 - added fast input tree processing
V1.4 - 17.06.2015 - binary refinement and github release
V1.5 - 30.09.2015 - weighted taxa cluster and frequency count 
V1.6 - 06.04.2016 - MPP module update
V1.7 - 09.04.2016 - Graph based cluster connectivity model incorporated
V1.8 - 05.06.2016 - Improved graph based cluster connectivity model with simplified code
"""

# Copyright 2013, 2014, 2015, 2016 
# Sourya Bhattacharyya and Jayanta Mukherjee.
# All rights reserved.
#
# See "LICENSE.txt" for terms and conditions of usage.
#
#---------------------------------------------

import Header
from Header import *
import Cost_Update
from Cost_Update import *
import Process_Queues
from Process_Queues import *
import ReachGraph_Update
from ReachGraph_Update import *
import UtilFunc
from UtilFunc import *
import RefineTree
from RefineTree import *
import Conflict_Detect
from Conflict_Detect import *
import Cluster_Manage
from Cluster_Manage import *

#-----------------------------------------------------
"""
this function is useful to parse various options for input data processing
"""
def parse_options():  
	parser = OptionParser()

	parser.add_option("-I", "--INPFILE", \
				type="string", \
				action="store", \
				dest="INP_FILENAME", \
				default="", \
				help="name of the input file (along with complete or relative path) containing input trees")

	parser.add_option("-O", "--OUTFILE", \
				type="string", \
				action="store", \
				dest="OUT_FILENAME", \
				default="", \
				help="name of the output file which will contain the output supertree")  

	parser.add_option("-p", "--inpform", \
				type="int", \
				action="store", \
				dest="inp_file_format", \
				default=1, \
				help="1 - input file format is NEWICK (default) \
				2 - input file format is NEXUS")

	parser.add_option("-b", "--binary", \
				action="store_true", \
				dest="binary_suptr", \
				default=False, \
				help="if TRUE, it produces a strictly binary supertree. \
				Otherwise, the tree can be non-binary. Default FALSE.")
					
	parser.add_option("-u", "--underscore", \
				action="store_false", \
				dest="preserve_underscores", \
				default=True, \
				help="this is a boolean flag option \
				using this option toggles the existing configuration (Default TRUE) \
				if TRUE, then this option preserves the underscores of the names of taxa \
				so, enabling this option do not preserve the underscores")  
					
	parser.add_option("-n", "--njrule", \
				type="int", \
				action="store", \
				dest="NJ_type", \
				default=2, \
				help="valid only if binary supertree is produced \
				1 - classical NJ method \
				2 - Normalized couplet statistic for agglomeration (Default)")     

	parser.add_option("-w", "--weighttaxa", \
				action="store_false", \
				dest="weight_taxa_subset", \
				default=True, \
				help="this is a boolean flag option \
				using this option toggles the existing configuration (Default TRUE) \
				if TRUE, then this option weighs couplet statistics according \
				to the size of taxa subset underlying MRCA of that couplet")  

	parser.add_option("-d", "--distmat", \
				type="int", \
				action="store", \
				dest="dist_mat_type", \
				default=1, \
				help="1 - Mean of XL \
				2 - Mean of internode count \
				3 - Mean of R1 reln Frequency")     

	opts, args = parser.parse_args()
	return opts, args

#-----------------------------------------------------
""" 
this is the main function of the package
"""
def main():  
	opts, args = parse_options()

	ROOTED_TREE = False #opts.default_rooted
	PRESERVE_UNDERSCORE = opts.preserve_underscores
	if (opts.inp_file_format == 1):
		INPUT_FILE_FORMAT = 'newick'
	else:
		INPUT_FILE_FORMAT = 'nexus'
	INPUT_FILENAME = opts.INP_FILENAME
	OUTPUT_FILENAME = opts.OUT_FILENAME
	BINARY_SUPERTREE_OPTION = opts.binary_suptr
	NJ_RULE_USED = opts.NJ_type 
	WEIGHT_TAXA_SUBSET = opts.weight_taxa_subset
	DIST_MAT_TYPE = opts.dist_mat_type
	
	if (INPUT_FILENAME == ""):
		print '******** THERE IS NO INPUT FILE SPECIFIED - RETURN **********'
		return
	else:
		print 'input filename: ', INPUT_FILENAME
	
	"""
	according to the location of input filename
	adjust the locations of the output files as well
	"""
	k = INPUT_FILENAME.rfind("/")
	if (k == -1):
		dir_of_inp_file = './'
	else:
		dir_of_inp_file = INPUT_FILENAME[:(k+1)]
	input_file_name = INPUT_FILENAME[(k+1):]
	if (DEBUG_LEVEL > 1):
		print 'dir_of_inp_file: ', dir_of_inp_file  
		
	#-------------------------
	if (OUTPUT_FILENAME == ""):
		"""
		first create the output directory containing the results
		the directory name should also reflect the settings of input parameters
		"""
		dir_of_curr_exec = dir_of_inp_file + 'COSPEDBTree_Try'
		"""
		according to the settings of input parameters, 
		output directory name is customized
		"""
		if (BINARY_SUPERTREE_OPTION == True):
			dir_of_curr_exec = dir_of_curr_exec + '_B_1'
		else:
			dir_of_curr_exec = dir_of_curr_exec + '_B_0'
		if (WEIGHT_TAXA_SUBSET == True):
			dir_of_curr_exec = dir_of_curr_exec + '_W_1'
		else:
			dir_of_curr_exec = dir_of_curr_exec + '_W_0'
		
		dir_of_curr_exec = dir_of_curr_exec + '_D_' + str(DIST_MAT_TYPE)
		
		if (BINARY_SUPERTREE_OPTION == True):
			dir_of_curr_exec = dir_of_curr_exec + '_N_' + str(NJ_RULE_USED)

		""" 
		append the current output directory in the text file
		"""
		Output_Text_File = dir_of_curr_exec + '/' + 'COSPEDBTree_Complete_Desription.txt'
		""" 
		create the directory
		"""
		if (os.path.isdir(dir_of_curr_exec) == False):
			mkdr_cmd = 'mkdir ' + dir_of_curr_exec
			os.system(mkdr_cmd)               
	
	else:
		k = OUTPUT_FILENAME.rfind("/")
		if (k == -1):
			dir_of_curr_exec = './'
		else:
			dir_of_curr_exec = OUTPUT_FILENAME[:(k+1)]
		Output_Text_File = dir_of_curr_exec + input_file_name + '_COSPEDBTree_Complete_Desription.txt'

	if (DEBUG_LEVEL > 1):
		print 'dir_of_curr_exec: ', dir_of_curr_exec  
		print 'Output_Text_File: ', Output_Text_File      

	#-------------------------
	# open the output text file
	fp = open(Output_Text_File, 'w')    
	
	fp.write('\n ================ status of options ================= (1 means ON)')
	fp.write('\n ROOTED_TREE: ' + str(ROOTED_TREE))
	fp.write('\n PRESERVE_UNDERSCORE: ' + str(PRESERVE_UNDERSCORE))
	fp.write('\n BINARY SUPERTREE OPTION: ' + str(BINARY_SUPERTREE_OPTION))
	fp.write('\n ===>>>  processing the file now ======== ')

	"""
	this file contains the execution timing information of the code
	"""
	timing_file = dir_of_curr_exec + '/' + 'timing_info.txt'
	fp_time = open(timing_file, 'w')

	# note the program beginning time 
	start_timestamp = time.time()
	#-------------------------------------  
	""" 
	read the source trees collection and store it in a treelist
	individual elements of this collection is a source tree 
	"""
	Input_Treelist = Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME)  
	#-------------------------------------  
	""" 
	from the input trees, note the number of taxa (total)
	and also define the class instances corresponding to single taxa
	"""
	for tr_idx in range(len(Input_Treelist)):
		taxa_labels_curr_tree = Input_Treelist[tr_idx].infer_taxa().labels()
		if (DEBUG_LEVEL > 1):
			fp.write('\n Tree no : ' + str(tr_idx+1) +  'no of leaf nodes: ' + str(len(taxa_labels_curr_tree)))
		if (DEBUG_LEVEL > 2):
			fp.write('\n taxa set belonging to current tree: ' + str(taxa_labels_curr_tree))
		for i in range(len(taxa_labels_curr_tree)):
			if taxa_labels_curr_tree[i] not in COMPLETE_INPUT_TAXA_LIST:
				COMPLETE_INPUT_TAXA_LIST.append(taxa_labels_curr_tree[i])

	"""
	we also define one structure "Taxa_Info_Dict" marked by a taxa
	this defines a single taxon, and becomes an instance of the class "Single_Taxa"
	The key of individual element is the index of the taxon in the COMPLETE_INPUT_TAXA_LIST
	"""
	for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
		Taxa_Info_Dict.setdefault(i, Single_Taxa())

	if (DEBUG_LEVEL >= 0):
		fp.write('\n  total no of taxa: ' + str(len(COMPLETE_INPUT_TAXA_LIST)))
		fp.write('\n  total no of trees: ' + str(len(Input_Treelist)))

	data_read_timestamp1 = time.time()	# note the timestamp

	if (DEBUG_LEVEL >= 0):
		fp_time.write('\n Time taken to read the taxa information: ' + str(data_read_timestamp1 - start_timestamp))
		
	#---------------------------
	"""
	here we add for individual taxon, the input trees supporting this taxon
	"""
	for tr_idx in range(len(Input_Treelist)):
		taxa_labels_curr_tree = Input_Treelist[tr_idx].infer_taxa().labels()
		for t in taxa_labels_curr_tree:
			taxa_idx = COMPLETE_INPUT_TAXA_LIST.index(t)
			Taxa_Info_Dict[taxa_idx]._AddSupportTree(tr_idx)

	#---------------------------
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n ***** Printing information of all taxa and their support tree count **** \n\n\n')
		for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
			fp.write('\n\n *** Taxa idx: ' + str(i) + '  Label: ' + str(COMPLETE_INPUT_TAXA_LIST[i]) + \
				'  Supporting tree list: ' + str(Taxa_Info_Dict[i]._GetSupportTreeList()))
		fp.close()
	#---------------------------
	"""
	for individual couplets, this function lists all the taxa 
	lying within their LCA node for all the input trees
	"""
	if (WEIGHT_TAXA_SUBSET == True):
		for tr_idx in range(len(Input_Treelist)):
			FindCoupletUnderlyingTaxon(Input_Treelist[tr_idx])

	data_read_timestamp2 = time.time()	# note the timestamp
	
	if (DEBUG_LEVEL >= 0):
		fp_time.write('\n\n Time taken to compute the underlying taxon for weight computation: ' \
			+ str(data_read_timestamp2 - start_timestamp))
	#---------------------------
	"""
	now process individual trees to find the couplet relations within those trees
	"""
	for tr_idx in range(len(Input_Treelist)):
		DeriveCoupletRelations(Input_Treelist[tr_idx], WEIGHT_TAXA_SUBSET)

	# note the timestamp
	data_read_timestamp3 = time.time()
	
	if (DEBUG_LEVEL >= 0):
		fp_time.write('\n Time taken to read the couplet relations: ' \
			+ str(data_read_timestamp3 - data_read_timestamp2))
	
	if (DEBUG_LEVEL > 2):
		fp.write('\n len Taxa_Info_Dict: ' + str(len(Taxa_Info_Dict)))
		fp.write('\n len COMPLETE_INPUT_TAXA_LIST: ' + str(COMPLETE_INPUT_TAXA_LIST))
		fp.write('\n len TaxaPair_Reln_Dict : ' + str(len(TaxaPair_Reln_Dict)))
	
	fp.close()
		
	# note the timestamp
	data_read_timestamp = time.time()
	#------------------------------------------------------------
	""" 
	Here we allocate the list / storage space of taxa clusters
	each taxa cluster is supposed to store the taxa subsets 
	related via relation r3 (simultaneous speciation)
	initially all the clusters contain one taxa
	each of the cluster has the index of the corresponding 
	taxa in the COMPLETE_INPUT_TAXA_LIST 
	"""
	for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
		Create_Cluster_Taxa_Label(i, COMPLETE_INPUT_TAXA_LIST[i])

	#------------------------------------------------------------
	"""
	here we analyze individual couplets and fill the support score queues
	with the supported relations
	"""
	Fill_Support_Score_Queues_Couplet_Based()
	
	#-----------------------------------------------------------
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n ***** Printing information of all the couplets **** \n\n\n')
		fp.close()
		for l in TaxaPair_Reln_Dict:
			TaxaPair_Reln_Dict[l]._PrintRelnInfo(l, Output_Text_File, WEIGHT_TAXA_SUBSET)
	
	#-------------------------------------------------------------------------------
	"""
	first process the queues containing the R3 relations 
	as a consensus or as a single occurring relation between a couplet
	to process corresponding couplets
	"""
	if (DEBUG_LEVEL > 0):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n *****  Now processing the following couplets \n 1) Only R3 relation supported, \
			2) R3 relation is majority consensus **** \n\n')
		fp.close()

	if (len(Queue_Score_R3_SingleReln) > 0):
		Sort_Priority_Queue(Queue_Score_R3_SingleReln)
		Proc_Queue_Couplet_Reln_R3(Output_Text_File, 1)

	if (len(Queue_Score_R3_MajCons) > 0):
		Sort_Priority_Queue(Queue_Score_R3_MajCons)
		Proc_Queue_Couplet_Reln_R3(Output_Text_File, 2)

	# note the timestamp
	data_initialize_timestamp = time.time()
	
	#------------------------------------------------------------
	"""
	this is the total number of taxa clusters
	individually containing one or more taxon
	which are themselves related via R3 relation
	"""
	nClusters = len(CURRENT_CLUST_IDX_LIST)
	
	#------------------------------------------------------------
	""" 
	we initialize the Reachability_Graph_Mat
	dimension: N X N where N = number of taxa clusters
	this is a numpy 2D array 
	values Mat[x][y] = 1 means x->y
	Mat[x][y] = Mat[y][x] = 2 means x and y are connected via relation r4
	"""
	Reachability_Graph_Mat = numpy.zeros((nClusters, nClusters), dtype=numpy.int)
	
	#------------------------------------------------------------
	if (DEBUG_LEVEL > 0):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n *****  Now processing the pair of clusters for their connectivity **** \n\n')
		fp.close()
	
	"""
	now we check individual cluster pairs, and check the allowed relations between them
	"""
	Initialize_Cluster_Pair(Output_Text_File)

	# print the cluster information 
	if (DEBUG_LEVEL > 0):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n ========== Initially formed taxa clusters (Initialize_Cluster_Pair) =============\n\n')
		fp.close()
		for i in Cluster_Info_Dict:
			#print 'printing the information for cluster node: ', i
			Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File, False)

	#----------------------------------------------------------
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n ***** Printing information of all the cluster pairs **** \n\n\n')
		fp.close()
		for l in Cluster_Pair_Info_Dict:
			r1_freq = Cluster_Pair_Info_Dict[l]._GetFreq(RELATION_R1)
			r2_freq = Cluster_Pair_Info_Dict[l]._GetFreq(RELATION_R2)
			r3_freq = Cluster_Pair_Info_Dict[l]._GetFreq(RELATION_R3)
			r4_freq = Cluster_Pair_Info_Dict[l]._GetFreq(RELATION_R4)
			if (r1_freq == 0) and (r2_freq == 0) and (r3_freq == 0) and (r4_freq > 0):
				continue
			else:
				Cluster_Pair_Info_Dict[l]._PrintClusterPairInfo(l, Output_Text_File)
	#----------------------------------------------------------
	"""
	form the support score queue for the cluster pairs having either R1 or R2 (or both)
	as their allowed relation
	"""
	Form_Conflict_Support_Score_Queue()
	#-------------------------------------------------------
	"""
	sorting the queue containing the support scores for individual cluster pairs
	then we process this support score queue and relate individual pair of clusters 
	"""
	if (len(Queue_Score_Cluster_Pair_NonConflict) > 0):
		Sort_Priority_Queue(Queue_Score_Cluster_Pair_NonConflict)
		Reachability_Graph_Mat = Proc_Queue_Clust(Reachability_Graph_Mat, Output_Text_File, 1)

	if (len(Queue_Score_Cluster_Pair) > 0):
		Sort_Priority_Queue(Queue_Score_Cluster_Pair)
		Reachability_Graph_Mat = Proc_Queue_Clust(Reachability_Graph_Mat, Output_Text_File, 2)

	#------------------------------------------------------------
	# print the cluster information 
	if (DEBUG_LEVEL > 0):
		fp = open(Output_Text_File, 'a')
		fp.write('\n **** total number of clusters: ' + str(len(CURRENT_CLUST_IDX_LIST)))
		fp.write('\n CURRENT_CLUST_IDX_LIST contents: ')
		fp.write(str(CURRENT_CLUST_IDX_LIST))
		fp.write('\n\n\n ========== cluster information after processing Queue_Score_Cluster_Pair =============\n\n')
		fp.close()
		for i in Cluster_Info_Dict:
			#print 'printing the information for cluster node: ', i
			Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)

	# note the timestamp
	reachability_graph_form_timestamp = time.time()  

	#------------------------------------------------------------
	""" 
	now perform the transitive reduction of the closure formed by 
	connection of the cluster of nodes in the above operation
	this is required to handle the following scenario:
	suppose, there exists a case such that A->C, B->C and A->B
	then in the final graph, only A->B and B->C information needs to be preserved
	in order to form the DAG 
	"""
	CompressDirectedGraph(Reachability_Graph_Mat, Output_Text_File)
	
	# print the cluster information 
	if (DEBUG_LEVEL > 0):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n ========== cluster information after transitive reduction =============\n\n')
		fp.close()
		for i in Cluster_Info_Dict:
			#print 'printing the information for cluster node: ', i
			Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)
	
	#------------------------------------------------------------
	"""
	now instead of arbitrary assignment of the parent node for individual clusters 
	we assign parent node according to the source tree relationships
	this will solve the multiple parent problem C2 as discussed in the manuscript 
	this is a new addition and marked under the DFS based parent refinement option 
	"""
	SelectUniqueParent_Directed(Output_Text_File, DIST_MAT_TYPE)
	
	# print the cluster information 
	if (DEBUG_LEVEL > 0):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n ========== cluster information after solving Transitive Reduction + MPP (for directed edges) =============\n\n')
		fp.close()
		for i in Cluster_Info_Dict:
			#print 'printing the information for cluster node: ', i
			Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)
	
	# note the timestamp
	cluster_of_node_refine_species_timestamp1 = time.time()  

	#----------------------------------------------
	""" 
	construct supertree from the generated DAG 
	scheme: repeatedly extract the nodes (taxa clusters) with minimum indegree
	print taxa information constituent within this cluster
	maintain appropriate paranthesis according to the newick format
	"""
	""" 
	the variable no_of_components is used in association with the used depth first technique 
	if no_of_components > 1, it signifies that a forest has been created
	in other words, the problem C3 (no parent problem exists)
	"""
	no_of_components = 0	# for forest
	while (1):
		root_clust_node_idx = Extract_Node_Min_Indeg(Output_Text_File, Cluster_Info_Dict)	# Cluster_Dict_Backup)
		if (root_clust_node_idx == -1):
			break
		Tree_Str = PrintNewick(root_clust_node_idx)
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n\n *** Tree_Str: ' + Tree_Str + '   *** \n\n')
			fp.close()
		no_of_components = no_of_components + 1
		if (no_of_components == 1):	# first component
			Final_Supertree_Str = Tree_Str
		else:
			Final_Supertree_Str = Final_Supertree_Str + ',' + Tree_Str
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n\n *** Final_Supertree_Str: ' + Final_Supertree_Str + '   *** \n\n')
			fp.close()

	""" 
	with the final tree string, finally generate the tree result 
	this is also required to tackle the creation of a forest during the DFS (no parent problem)
	this procedure adheres to the basic COSPEDTree mechanism
	"""
	Final_Supertree_Str = '(' + Final_Supertree_Str + ')'

	#---------------------------------------------------------
	fp = open(Output_Text_File, 'a')
	fp.write('\n\n\n\n **** original supertree as newick string --- ' + Final_Supertree_Str) 

	"""
	now read this super string in a supertree containing all the input taxa
	Note: This is an important step, since the parsed tree from the string has the correct 
	tree structure, with respect to dendropy format
	"""
	Supertree_Final = dendropy.Tree.get_from_string(Final_Supertree_Str, schema="newick")	# #preserve_underscores=PRESERVE_UNDERSCORE, default_as_rooted=ROOTED_TREE)
	if 0:
		Supertree_Final.print_plot()  
		
	# note the timestamp
	newick_str_formation_timestamp = time.time()  

	fp = open(Output_Text_File, 'a')
	#fp.write('\n initial supertree: ' + Supertree_without_branch_len.as_newick_string())    
	Supertree_Final.update_splits(delete_outdegree_one=True)
	fp.write('\n\n\n ***** after update splits: supertree: ' + Supertree_Final.as_newick_string())    
	fp.close()

	if (BINARY_SUPERTREE_OPTION == True):    
		""" 
		this function removes all multifurcating clusters and produces binary tree 
		it also solves the problem C3, as mentioned in the manuscript
		"""
		Supertree_Final = Refine_Supertree_Binary_Form(Input_Treelist, Supertree_Final, Output_Text_File, NJ_RULE_USED, DIST_MAT_TYPE)
		fp = open(Output_Text_File, 'a')
		fp.write('\n --- user provided option for producing strict binary supertree')
		fp.write('\n --- after binary refinement --- output tree without branch length (in newick format): ' + Supertree_Final.as_newick_string())    
		fp.close()
	else:
		fp = open(Output_Text_File, 'a')
		fp.write('\n --- user did not provide option for producing strict binary supertree - so output tree can be non-binary')
		fp.close()
		
	"""
	write this tree on a separate text file
	"""
	if (OUTPUT_FILENAME == ""):
		out_treefilename = dir_of_curr_exec + '/' + 'cospedbtree_newick.tre'
	else:
		out_treefilename = OUTPUT_FILENAME

	Supertree_Final.write_to_path(out_treefilename, 'newick')
	#---------------------------------
	"""
	Note: This seems to be a redundant operation, but an important one
	Somehow, the structure "Supertree_Final" fails to read the underscores properly
	so we read the tree from the file itself, and store it again in the "Supertree_Final" structure
	"""
	Supertree_Final = dendropy.Tree.get_from_path(out_treefilename, schema='newick', preserve_underscores=PRESERVE_UNDERSCORE)

	# final timestamp
	data_process_timestamp = time.time()      

	#---------------------------------
	"""
	timing information write
	"""
	fp_time.write('\n \n\n ===============>>>>>>>>>>>>>>> TIME COMPLEXITY OF THE METHOD (in seconds) ')
	fp_time.write('\n \n reading the data: ' + str(data_read_timestamp - start_timestamp) + \
	'\n initialization of the structure: ' + str(data_initialize_timestamp - data_read_timestamp) + \
	'\n formation of the reachability graph (cluster) (after loop): ' + \
				str(reachability_graph_form_timestamp - data_initialize_timestamp) + \
	'\n multiple parent (related) problem: ' + \
				str(cluster_of_node_refine_species_timestamp1 - reachability_graph_form_timestamp) + \
	'\n newick string formation: ' + str(newick_str_formation_timestamp - cluster_of_node_refine_species_timestamp1) + \
		'\n binary tree construction: ' + str(data_process_timestamp - newick_str_formation_timestamp))

	fp_time.write('\n \n Total time taken (in seconds) : ' + str(data_process_timestamp - start_timestamp))  
	fp_time.close()

	#----------------------------------------------
	# Performance metric code
	#----------------------------------------------
	# open the output text file
	outtextfile = dir_of_curr_exec + '/' + 'FP_FN_RF_Perf.txt'
	fp = open(outtextfile, 'w')

	# examine each of the source trees and find the FP, FN and RF distance with respect to the generated supertree  
	sumFP = sumFN = sumRF = 0  
	sumLenSrcTree = 0
	sum_symmetric_diff = 0
	fp.write('\n \n\n total edges of supertree: ' + str(len(Supertree_Final.get_edge_set())))  


	# add - sourya - debug
	Input_Treelist = Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME)  

	for inp_tree_idx in range(len(Input_Treelist)):
		Curr_src_tree = Input_Treelist[inp_tree_idx]
		curr_src_tree_taxa = Curr_src_tree.infer_taxa().labels()
		curr_src_tree_no_of_taxa = len(curr_src_tree_taxa)
		
		"""
		according to the taxa set of the current source tree, 
		prune the supertree to get the tree portion containing only this taxa set
		"""
		pruned_tree = dendropy.Tree(Supertree_Final)
		pruned_tree.retain_taxa_with_labels(curr_src_tree_taxa)
		
		"""
		source tree number of edges calculation
		it is used to compute normalized RF metric values
		"""
		lenSrcTree = len(Curr_src_tree.get_edge_set())
		sumLenSrcTree = sumLenSrcTree + lenSrcTree
		fp.write('\n\n\n src tree : ' + str(Curr_src_tree))
		fp.write('\n\n pruned supertree: ' + str(pruned_tree))
		fp.write('\n\n src tree len: ' + str(lenSrcTree) + ' pruned supertree len: ' + str(len(pruned_tree.get_edge_set())))
		
		"""
		determine the false positives and the false negatives 
		"""
		tup = Curr_src_tree.false_positives_and_negatives(pruned_tree)
		fp.write('   FP_int: ' + str(tup[0]) + '  FN_int:  ' + str(tup[1]))
		sumFP = sumFP + tup[0]
		sumFN = sumFN + tup[1]
		sumRF = sumRF + ((tup[0] + tup[1]) / 2.0)
		
		symm_diff = Curr_src_tree.symmetric_difference(pruned_tree)
		fp.write('   Symmetric difference: ' + str(symm_diff))
		sum_symmetric_diff = sum_symmetric_diff + symm_diff

	"""
	final normalized sumFP's are computed by dividing with the number of trees
	"""
	normsumFP = (sumFP * 1.0) / sumLenSrcTree
	normsumFN = (sumFN * 1.0) / sumLenSrcTree
	normsumRF = (sumRF * 1.0) / sumLenSrcTree
	norm_symm_diff = sum_symmetric_diff / (2.0 * sumLenSrcTree)
	
	"""
	print the final result
	"""
	fp.write('\n\n\n ===============>>>>>>>>>>>>>>> FINAL RESULTS \n \n')
	fp.write('\n ******* absolute sumFP: ' + str(sumFP) + \
		'\n ******* absolute sumFN: ' + str(sumFN) + \
		'\n ******* absolute sumRF: ' + str(sumRF) + \
		'\n ******* absolute Symmetric difference: ' + str(sum_symmetric_diff))

	fp.write('\n ===============>>>>>>>>>>>>>>> IN TERMS OF NORMALIZED \
		(DIVIDED BY THE SUM OF INTERNAL EDGES OF THE SOURCE TREES) ''')  
	fp.write('\n normsumFP: ' + str(normsumFP) + '\n normsumFN: ' + str(normsumFN) + \
		'\n normsumRF: ' + str(normsumRF) + '\n norm Symmetric Diff: ' + str(norm_symm_diff)) 

	fp.close()
	#----------------------------------------------
	# end Performance metric code
	#----------------------------------------------
	#--------------------------------------------------------------  
	# delete the storage variables associated with the current execution 

	# clear the dictionaries
	Cluster_Info_Dict.clear()
	Taxa_Info_Dict.clear()
	TaxaPair_Reln_Dict.clear()
	Cluster_Pair_Info_Dict.clear()
	
	# clear the lists associated
	if (len(Queue_Score_R3_SingleReln) > 0):
		Queue_Score_R3_SingleReln[:] = []
	if (len(Queue_Score_R3_MajCons) > 0):
		Queue_Score_R3_MajCons[:] = []
	if (len(Queue_Score_Cluster_Pair) > 0):
		Queue_Score_Cluster_Pair[:] = []
	if (len(Queue_Score_Cluster_Pair_NonConflict) > 0):
		Queue_Score_Cluster_Pair_NonConflict[:] = []
	if (len(COMPLETE_INPUT_TAXA_LIST) > 0):
		COMPLETE_INPUT_TAXA_LIST[:] = []
	if (len(CURRENT_CLUST_IDX_LIST) > 0):
		CURRENT_CLUST_IDX_LIST[:] = []

	# free the reachability graph (numpy array)
	del Reachability_Graph_Mat
	#del Clust_Possible_R1_Mat
	
	# free other arrays
	Input_Treelist = []
  
#-----------------------------------------------------
if __name__ == "__main__":
	main() 
  
