#!/usr/bin/env python


##---------------------------------------------
''' 
this program is used to generate a supertree (consensus) from a set of constituent trees
the input is multiple source trees
each of the trees need to be decomposed to couplets 
then the couplets need to be joined
there may be conflicts among the input tree - we have to select the consensus

Author: Sourya Bhattacharyya
Dept of CSE, IIT Kharagpur
V1.0 - 15.01.2014 - public release
V2.0 - 15.03.2014 - rewritten for lowering time complexity
V3.0 - 28.02.2015 - added binary suprtree and level based scoring options
V4.0 - 30.03.2015 - added fast input tree processing
V5.0 - 17.06.2015 - binary refinement and github release
''' 

## Copyright 2013 Sourya Bhattacharyya and Jayanta Mukherjee.
## All rights reserved.
##
## See "LICENSE.txt" for terms and conditions of usage.
##
##---------------------------------------------

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

##-----------------------------------------------------
# this function is useful to parse various options for input data processing
def parse_options():  
  parser = OptionParser()
  
  parser.add_option("-I", "--INPFILE", \
			  type="string", \
			  action="store", \
			  dest="INP_FILENAME", \
			  default="", \
			  help="name of the input file containing candidate source trees")
  
  parser.add_option("-O", "--OUTFILE", \
			  type="string", \
			  action="store", \
			  dest="OUT_FILENAME", \
			  default="", \
			  help="name of the output file which will contain the target supertree")  
  
  parser.add_option("-p", "--inpform", \
			  type="int", \
			  action="store", \
			  dest="inp_file_format", \
			  default=1, \
			  help="1 - input file format is NEWICK (default) \
			  2 - input file format is NEXUS")
  
  parser.add_option("-q", "--queues", \
			  type="int", \
			  action="store", \
			  dest="no_of_queues", \
			  default=2, \
			  help="1 - only a single max priority queue is used for storing the score metrics \
			  2 - two separate queues are used to store the conflicting and non conflicting taxa pairs and corresponding score metrics (default)")

  parser.add_option("-t", "--treewt", \
			  type="int", \
			  action="store", \
			  dest="tree_weight", \
			  default=0, \
			  help="0 - there is no separate weight assigned to individual phylogenetic trees \
			  1 - input phylogenetic trees are assigned separate weights (default)")

  parser.add_option("-d", "--dynscore", \
			  action="store_true", \
			  dest="dynamic_score", \
			  default=False, \
			  help="if true, then this option dynamically updates support score values - Default FALSE")

  parser.add_option("-b", "--binary", \
			  action="store_false", \
			  dest="binary_suptr", \
			  default=True, \
			  help="if TRUE, it produces a strictly binary supertree. Otherwise, the tree can be non-binary. Default FALSE.")
			    
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
			  1 - classical NJ method (Default) \
			  2 - Normalized couplet statistic for agglomeration")     
			    
  parser.add_option("-s", "--supportcoupletrule", \
			  action="store_false", \
			  dest="Support_Couplet_Check", \
			  default=True, \
			  help="this is a boolean flag option \
			  valid only if binary supertree is produced \
			  using this option toggles the existing configuration (Default TRUE) \
			  if TRUE, then this option checks minimum XL entry for only supported couplets")  
			    
  parser.add_option("-m", "--metric", \
			  type="int", \
			  action="store", \
			  dest="dist_metric", \
			  default=1, \
			  help="valid only if binary supertree is produced \
			  1 - sum of extra taxa (XT) with respect to individual input trees \
			  2 - Accumulated branch count measure (employed in NJ_st method) \
			  3 - product of branch count and excess taxa \
			  4 - product of LCA rank and excess taxa")
			    
  parser.add_option("-c", "--classavg", \
			  type="int", \
			  action="store", \
			  dest="class_of_metric", \
			  default=2, \
			  help="valid only if binary supertree is produced \
			  1 - absolute sum of metric value (either XT or Branch count) between couplets\
			  2 - simple average of  metric value (either XT or Branch count) between couplets \
			  3 - mode based average of metric value (either XT or Branch count) between couplets")     
			
  opts, args = parser.parse_args()
  return opts, args
  
  
##-----------------------------------------------------
''' main function '''
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
  NO_OF_QUEUES = opts.no_of_queues
  DFS_PARENT_REFINE = True #opts.dfs_parent_refine
  BINARY_SUPERTREE_OPTION = opts.binary_suptr
  NJ_RULE_USED = opts.NJ_type
  SUPPORT_COUPLET_CHECK = opts.Support_Couplet_Check
  DIST_METRIC = opts.dist_metric
  CLASS_OF_METRIC = opts.class_of_metric
  TREE_WEIGHT_ASSIGNED = opts.tree_weight
  DYNAMIC_SCORE_UPDATE = opts.dynamic_score
  
  if (INPUT_FILENAME == ""):
    print '******** THERE IS NO INPUT FILE SPECIFIED - RETURN **********'
    return
  else:
    print 'input filename: ', INPUT_FILENAME
        
  # according to the location of input filename
  # adjust the locations of the output files as well
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
    # first create the output directory containing the results
    dir_of_curr_exec = dir_of_inp_file + 'COSPEDBTree'
    # update dir_of_curr_exec according to the settigs of input parameters
    dir_of_curr_exec = dir_of_curr_exec + '_Q_' + str(NO_OF_QUEUES)
    if (BINARY_SUPERTREE_OPTION == True):
      dir_of_curr_exec = dir_of_curr_exec + '_B_1'
    else:
      dir_of_curr_exec = dir_of_curr_exec + '_B_0'

    if (TREE_WEIGHT_ASSIGNED == True):
      dir_of_curr_exec = dir_of_curr_exec + '_T_1'
    else:
      dir_of_curr_exec = dir_of_curr_exec + '_T_0'
    
    if (DYNAMIC_SCORE_UPDATE == True):
      dir_of_curr_exec = dir_of_curr_exec + '_D_1'
    else:
      dir_of_curr_exec = dir_of_curr_exec + '_D_0'
        
    # following options are used if binary supertree is produced
    if (BINARY_SUPERTREE_OPTION == True):
      dir_of_curr_exec = dir_of_curr_exec + '_N_' + str(NJ_RULE_USED)
      if (SUPPORT_COUPLET_CHECK == True):
	dir_of_curr_exec = dir_of_curr_exec + '_S_1'
      else:
	dir_of_curr_exec = dir_of_curr_exec + '_S_0'
      dir_of_curr_exec = dir_of_curr_exec + '_M_' + str(DIST_METRIC)
      dir_of_curr_exec = dir_of_curr_exec + '_C_' + str(CLASS_OF_METRIC)
    
    # append the current output directory in the text file
    Output_Text_File = dir_of_curr_exec + '/' + 'COSPEDBTree_Complete_Desription.txt'
    # create the directory
    if (os.path.isdir(dir_of_curr_exec) == False):
      mkdr_cmd = 'mkdir ' + dir_of_curr_exec
      os.system(mkdr_cmd)               
    if (DEBUG_LEVEL > 1):
      print 'Output_Text_File: ', Output_Text_File      
  else:
    k = OUTPUT_FILENAME.rfind("/")
    if (k == -1):
      dir_of_curr_exec = './'
    else:
      dir_of_curr_exec = OUTPUT_FILENAME[:(k+1)]
    Output_Text_File = dir_of_curr_exec + input_file_name + '_COSPEDBTree_Complete_Desription.txt'
  
  #-------------------------
  fp = open(Output_Text_File, 'w')    
  #fp.write('\n ================ status of options ================= (1 means ON)')
  #fp.write('\n ROOTED_TREE: ' + str(ROOTED_TREE))
  #fp.write('\n PRESERVE_UNDERSCORE: ' + str(PRESERVE_UNDERSCORE))
  fp.write('\n NO_OF_QUEUES: ' + str(NO_OF_QUEUES))
  fp.write('\n MULTIPLE PARENT PROBLEM C2 solve: ' + str(DFS_PARENT_REFINE))
  fp.write('\n BINARY SUPERTREE OPTION: ' + str(BINARY_SUPERTREE_OPTION))
  fp.write('\n ===>>>  processing the file now ======== ')
  
  # this variable notes the count of input source trees  
  tree_count = 0
    
  # note the program beginning time 
  start_timestamp = time.time()
    
  #-------------------------------------  
  """ read the source trees collection and store it in a tree collection structure
  individual elements of this collection is thus a source tree """
  Input_Treelist = Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME)  
  
  #-------------------------------------  
  # from the input source trees, note the number of taxa (total)
  # and also define the class instances corresponding to single taxa
  for tr_idx in range(len(Input_Treelist)):
    tree_count = tree_count + 1
    taxa_labels_curr_tree = Input_Treelist[tr_idx].infer_taxa().labels()
    if (DEBUG_LEVEL > 1):
      fp.write('\n Tree no : ' + str(tree_count) +  'no of leaf nodes: ' + str(len(taxa_labels_curr_tree)))
    if (DEBUG_LEVEL > 2):
      fp.write('\n taxa set belonging to current tree: ' + str(taxa_labels_curr_tree))
    for i in range(len(taxa_labels_curr_tree)):
      if taxa_labels_curr_tree[i] not in COMPLETE_INPUT_TAXA_LIST:
	COMPLETE_INPUT_TAXA_LIST.append(taxa_labels_curr_tree[i])
    # add - sourya
    if (TREE_WEIGHT_ASSIGNED == 1):
      for i in range(len(taxa_labels_curr_tree) - 1):
	for j in range(i+1, len(taxa_labels_curr_tree)):
	  key1 = (taxa_labels_curr_tree[i], taxa_labels_curr_tree[j])
	  key2 = (taxa_labels_curr_tree[j], taxa_labels_curr_tree[i])
	  if key1 in TaxaPair_Reln_Dict:
	    TaxaPair_Reln_Dict[key1]._AddSupportingTree()
	  elif key2 in TaxaPair_Reln_Dict:
	    TaxaPair_Reln_Dict[key2]._AddSupportingTree()
	  else:
	    TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
	    TaxaPair_Reln_Dict[key1]._AddSupportingTree()
    # end add - sourya
  
  # we also define one structure "Taxa_Info_Dict" marked by a taxa
  for label in COMPLETE_INPUT_TAXA_LIST:
    Taxa_Info_Dict.setdefault(label, Single_Taxa())

  # close the text file
  fp.close()
  #---------------------------------------------
  if (TREE_WEIGHT_ASSIGNED == 1):
    AssignMatrixWeights(Input_Treelist)
    PrintMatrixWeights(Output_Text_File)
    
  # now process individual trees to find the couplet relations of those trees
  for tr_idx in range(len(Input_Treelist)):
    DeriveCoupletRelations(Input_Treelist[tr_idx], tr_idx, TREE_WEIGHT_ASSIGNED, DYNAMIC_SCORE_UPDATE)
  
  fp = open(Output_Text_File, 'a')    
  if (DEBUG_LEVEL > 0):
    fp.write('\n  total no of taxa: ' + str(len(COMPLETE_INPUT_TAXA_LIST)))
  if (DEBUG_LEVEL > 1):
    fp.write('\n len Taxa_Info_Dict: ' + str(len(Taxa_Info_Dict)))
    fp.write('\n len COMPLETE_INPUT_TAXA_LIST: ' + str(COMPLETE_INPUT_TAXA_LIST))
    fp.write('\n len TaxaPair_Reln_Dict : ' + str(len(TaxaPair_Reln_Dict)))
  fp.close()
    
  data_read_timestamp = time.time()	# note the timestamp
    
  #------------------------------------------------------------
  ''' we also calculate the connection value between each pair of nodes in the output tree
  the value defines the majority of the edge type that is between those 2 nodes '''
  for l in TaxaPair_Reln_Dict:
    """
    calculate the consensus and priority measures associated with each couplet for different relations
    single_edge_exist means that the couplet is non-conflicting
    only one type of relation exists between them in the input trees
    in such a case, include only that relation in the priority queue
    """
    single_edge_exist_list = TaxaPair_Reln_Dict[l]._SetConnPrVal(True)
    single_edge_exist = single_edge_exist_list[0]
    edge_type = single_edge_exist_list[1]

    """ 
    we calculate the support score value between individual couplets and for each relations
    previously the support score value was equal to the priority metric
    now we change it to make it a product of the frequency and the priority measures 
    """
    TaxaPair_Reln_Dict[l]._SetCostMetric(l[0], l[1])
  
    #------------------------------------------------------------
    """ also update the cost value for individual elements in the list "Cost_List_Taxa_Pair_Multi_Reln or Cost_List_Taxa_Pair_Single_Reln"
    each list element contains the following values:
    1) taxa1 and taxa2    
    3) edge type (one at a time - so there will be 4 entries for each node pair)
    4) edge cost (the cost associated with one particular edge type """
    if (single_edge_exist == 0):
      for edge_type in range(4):
	# now we add only those relations between the taxa pair which are supported by at least one source tree
	if (TaxaPair_Reln_Dict[l]._GetEdgeWeight(edge_type) > 0):
	  sublist = [l[0], l[1], edge_type, TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(edge_type)]
	  Cost_List_Taxa_Pair_Multi_Reln.append(sublist)
    else:
      # for single edge type, assign that particular connection
      sublist = [l[0], l[1], edge_type, TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(edge_type)]
      """ this connection is only possible between the current examined taxa pairs
      if SINGLE_EDGE_TYPE_CONN_PRIORITY is TRUE then the information is placed in Cost_List_Taxa_Pair_Single_Reln
      otherwise it is placed in Cost_List_Taxa_Pair_Multi_Reln """
      if (NO_OF_QUEUES == 2):
	Cost_List_Taxa_Pair_Single_Reln.append(sublist)
      else:
	Cost_List_Taxa_Pair_Multi_Reln.append(sublist)
    #------------------------------------------------------------
  
  # we print the original connection status for all the tree nodes
  if (DEBUG_LEVEL > 2):
    for l in TaxaPair_Reln_Dict:
      #print 'printing info for the TaxaPair_Reln_Dict key: ', l
      TaxaPair_Reln_Dict[l]._PrintRelnInfo(l, Output_Text_File)
    
  #------------------------------------------------------------
  """ here we allocate the list of clusters
  initially all the clusters contain one taxa
  each of the cluster has the index of the corresponding taxa in the COMPLETE_INPUT_TAXA_LIST """
  for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
    Create_Cluster_Taxa_Label(i, COMPLETE_INPUT_TAXA_LIST[i])
  
  #------------------------------------------------------------
  """ we initialize the Reachability_Graph_Mat
  for all the clusters of nodes possible, this indicate the edges between a cluster pair
  initially the number of clusters are thought of equal as the number of taxa
  gradually, as the taxa subsets are merged in a cluster, corresponding count is decreased 
  convention -  out edge from the first node to the second node 
  this is a numpy 2D array 
  values Mat[x][y] = 1 means x->y
  Mat[x][y] = Mat[y][x] = 2 means x and y are connected via NO EDGE """
  Reachability_Graph_Mat = numpy.zeros((len(COMPLETE_INPUT_TAXA_LIST), len(COMPLETE_INPUT_TAXA_LIST)), dtype=numpy.int)
  
  fp = open(Output_Text_File, 'a')
  if (DEBUG_LEVEL > 0):
    fp.write('\n shape of Reachability_Graph_Mat: ' + str(numpy.shape(Reachability_Graph_Mat)))
  
  # this information is printed to know the maximum possible iterations that the while loops will undergo
  if (DEBUG_LEVEL > 1):
    fp.write('\n =========== max connection pair ============= : ' + str((len(Cost_List_Taxa_Pair_Single_Reln) + len(Cost_List_Taxa_Pair_Multi_Reln))))      
    fp.write('\n len Cost_List_Taxa_Pair_Single_Reln: ' + str(len(Cost_List_Taxa_Pair_Single_Reln)))
    fp.write('\n len Cost_List_Taxa_Pair_Multi_Reln : ' + str(len(Cost_List_Taxa_Pair_Multi_Reln)))
  
  fp.close()
  #------------------------------------------------------------
  """ now we have to sort the Cost_List_Node_Pair according to the edge cost value 
  that is the 4th field of the sublist 
  we use custom sorting operation """
  Sort_Priority_Queue(Cost_List_Taxa_Pair_Multi_Reln)
  # if there is provision to include single connectivity edges then we sort that list also
  if (NO_OF_QUEUES == 2):
    Sort_Priority_Queue(Cost_List_Taxa_Pair_Single_Reln)
  
  data_initialize_timestamp = time.time()	# note the timestamp
  
  #------------------------------------------------------------
  # print the queue storing the scores of individual relations
  if (DEBUG_LEVEL > 2):
    if (NO_OF_QUEUES == 2):
      fp = open(Output_Text_File, 'a')
      fp.write('\n printing contents for the non-conflicting queue (couplet relation score)')
      fp.close()
      PrintQueueInfo(Cost_List_Taxa_Pair_Single_Reln, Output_Text_File)
      
    fp = open(Output_Text_File, 'a')
    fp.write('\n printing contents for the conflicting queue (couplet relation score)')
    fp.close()
    PrintQueueInfo(Cost_List_Taxa_Pair_Multi_Reln, Output_Text_File)
    
  #------------------------------------------------------------
  if (NO_OF_QUEUES == 2):
    Reachability_Graph_Mat = Proc_Queue(Reachability_Graph_Mat, 1, Output_Text_File, DYNAMIC_SCORE_UPDATE)
  
  """ then we process the queue containing the taxa pairs depicting multiple relation instance """
  Reachability_Graph_Mat = Proc_Queue(Reachability_Graph_Mat, 0, Output_Text_File, DYNAMIC_SCORE_UPDATE)
  
  # we print the final connection status for all the tree nodes
  if (DEBUG_LEVEL > 2):
    for l in Taxa_Info_Dict:
      #print 'printing information for the Taxa ', l
      Taxa_Info_Dict[l]._PrintFinalTaxaInfo(l, Output_Text_File) 
  
  # print the cluster information 
  if (DEBUG_LEVEL > 0):
    fp = open(Output_Text_File, 'a')
    fp.write('\n **** total number of clusters: ' + str(len(CURRENT_CLUST_IDX_LIST)))
    fp.write('\n CURRENT_CLUST_IDX_LIST contents: ')
    fp.write(str(CURRENT_CLUST_IDX_LIST))
    fp.write('\n ========== cluster information after reachability graph generation =============')
    fp.close()
    for i in Cluster_Info_Dict:
      #print 'printing the information for cluster node: ', i
      Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)
  
  # note the timestamp
  reachability_graph_form_timestamp = time.time()  
  #------------------------------------------------------------
  """ now perform the transitive reduction of the closure formed by connection of the cluster of nodes in the above operation
  this is required to handle the following scenario:
  suppose, there exists a case such that A->C, B->C and A->B
  then in the final graph, only A->B and B->C information needs to be preserved
  in order to form the DAG """
  CompressDirectedGraph(Reachability_Graph_Mat)
    
  #------------------------------------------------------------
  # now instead of arbitrary assignment of the parent node for individual clusters 
  # we assign parent node according to the source tree relationships
  # this will solve the multiple parent problem C2 as discussed in the manuscript 
  # this is a new addition and marked under the DFS based parent refinement option 
  if (DFS_PARENT_REFINE == True):
    SolveMultipleParentC2Problem(Output_Text_File)
        
  # print the cluster information 
  if (DEBUG_LEVEL > 0):
    fp = open(Output_Text_File, 'a')
    fp.write('\n **** total number of clusters: ' + str(len(CURRENT_CLUST_IDX_LIST)))
    fp.write('\n CURRENT_CLUST_IDX_LIST contents: ')
    fp.write(str(CURRENT_CLUST_IDX_LIST))    
    fp.write('\n ========== cluster information after transitive reduction =============')
    fp.close()
    for i in Cluster_Info_Dict:
      #print 'printing the information for cluster node: ', i
      Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)
      
  # note the timestamp
  cluster_of_node_refine_species_timestamp1 = time.time()  
  
  ##----------------------------------------------
  ## add - sourya
  #if (BINARY_SUPERTREE_OPTION == True):
    #fp = open(Output_Text_File, 'a')
    #fp.write('\n --- We first refine at cluster level - clusters containing species more than 2 will be divided')
    #fp.close()
  
    #Refine_Clusters_Binary(Output_Text_File)

    ## print the cluster information 
    #if (DEBUG_LEVEL > 0):
      #fp = open(Output_Text_File, 'a')
      #fp.write('\n **** after modification of the clusters \n total number of clusters: ' + str(len(CURRENT_CLUST_IDX_LIST)))
      #fp.write('\n CURRENT_CLUST_IDX_LIST contents: ')
      #fp.write(str(CURRENT_CLUST_IDX_LIST))    
      #fp.write('\n ========== cluster information after division of constituent taxa =============')
      #fp.close()
      #for i in Cluster_Info_Dict:
	##print 'printing the information for cluster node: ', i
	#Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)
  
  ## end add - sourya
  ##----------------------------------------------
  ''' now this section constructs the supertree from the generated DAG 
  this is performed by repeatedly extracting the nodes with minimum indegree
  basically we first form a string which represents the supertree '''
  no_of_components = 0	# for forest
  while (1):
    root_clust_node_idx = Extract_Node_Min_Indeg(len(CURRENT_CLUST_IDX_LIST))
    if (root_clust_node_idx == -1):
      break
    Tree_Str = PrintNewick(root_clust_node_idx)	#, Reachability_Graph_Mat, len(CURRENT_CLUST_IDX_LIST))
    no_of_components = no_of_components + 1
    if (no_of_components == 1):	# first component
      Final_Supertree_Str = Tree_Str
    else:
      Final_Supertree_Str = Final_Supertree_Str + ',' + Tree_Str
  
  # with the final tree string, finally generate the tree result 
  Final_Supertree_Str = '(' + Final_Supertree_Str + ')'
    
  fp = open(Output_Text_File, 'a')
  fp.write('\n --- original supertree as newick string --- ' + Final_Supertree_Str) 
  
  Final_Supertree_Str = Remove_Extra_Paranthesis(Final_Supertree_Str)  
  fp.write('\n --- after removing extra paranthesis -- supertree as newick string --- ' + Final_Supertree_Str) 
  fp.close()
  
  # now read this super string in a supertree containing all the input taxa
  Supertree_Final = dendropy.Tree.get_from_string(Final_Supertree_Str, schema="newick")	#preserve_underscores=PRESERVE_UNDERSCORE, default_as_rooted=ROOTED_TREE)
  if 0:
    Supertree_Final.print_plot()  
    
  # note the timestamp
  newick_str_formation_timestamp = time.time()  
  
  # add - sourya  
  if (BINARY_SUPERTREE_OPTION == True):
    # ADD - SOURYA
    # we apply the update splits routine so as to remove any internal node with outdegree 1
    # this is required since during DAG based supertree formation, many internal nodes with one outdegree is created
    # this causes problem when binary refinement is employed
    #Supertree_Final.update_splits(delete_outdegree_one=True)
    #fp.write('\n --- after update splits --- output tree without branch length (in newick format): ' + Supertree_Final.as_newick_string())    
    
    ## add - sourya
    ## we add the function to refine binary of supertree
    #Refine_Latest_Supertree_Binary(Supertree_Final, Output_Text_File)
    ## end add - sourya
        
    # this function removes all multifurcating clusters and produces binary tree (except problem C3)
    Refine_Supertree_Binary_Form(Supertree_Final, NJ_RULE_USED, SUPPORT_COUPLET_CHECK, DIST_METRIC, CLASS_OF_METRIC, Output_Text_File)
    fp = open(Output_Text_File, 'a')
    fp.write('\n --- user provided option for producing strict binary supertree')
    fp.write('\n --- after binary refinement --- output tree without branch length (in newick format): ' + Supertree_Final.as_newick_string())    
    fp.close()
  else:
    fp = open(Output_Text_File, 'a')
    fp.write('\n --- user did not provide option for producing strict binary supertree - so output tree can be non-binary')
    fp.close()
    
  # write this tree on a separate text file
  if (OUTPUT_FILENAME == ""):
    out_treefilename = dir_of_curr_exec + '/' + 'cospedbtree_newick.tre'
  else:
    out_treefilename = OUTPUT_FILENAME
  
  ## comment - sourya
  #outfile = open(out_treefilename, 'w')
  #outfile.write(Supertree_Final.as_newick_string())
  ##outfile.write('\n \n final tree \n \n')
  ##outfile.write(Supertree_Final.as_ascii_plot())
  #outfile.close()
  ## end comment - sourya

  # add - sourya
  Supertree_Final.write_to_path(out_treefilename, 'newick')
  # end add - sourya
  
  # read the tree from the file itself
  Supertree_Final = dendropy.Tree.get_from_path(out_treefilename, schema='newick', preserve_underscores=PRESERVE_UNDERSCORE)
  
  # final timestamp
  data_process_timestamp = time.time()      
    
  #----------------------------------------------
  # Performance metric code
  #----------------------------------------------
  if 1:
    # open the output text file
    fp = open(Output_Text_File, 'a')
    
    # examine each of the source trees and find the FP, FN and RF distance with respect to the generated supertree  
    sumFP = sumFN = sumRF = 0  
    sumLenSrcTree = 0
    sum_symmetric_diff = 0
    fp.write('\n \n\n total edges of supertree: ' + str(len(Supertree_Final.get_edge_set())))  
    
    #print 'taxon set of supertree: ', Supertree_Final.infer_taxa()
    for inp_tree_idx in range(len(Input_Treelist)):
      Curr_src_tree = Input_Treelist[inp_tree_idx]
      curr_src_tree_taxa = Curr_src_tree.infer_taxa().labels()
      curr_src_tree_no_of_taxa = len(curr_src_tree_taxa)
      
      # according to the taxa set of the current source tree, 
      # prune the supertree to get the tree portion containing only this taxa set
      pruned_tree = dendropy.Tree(Supertree_Final)
      pruned_tree.retain_taxa_with_labels(curr_src_tree_taxa)
      
      # source tree number of edges calculation
      # it is used to compute normalized RF metric values
      lenSrcTree = len(Curr_src_tree.get_edge_set())
      sumLenSrcTree = sumLenSrcTree + lenSrcTree
      fp.write('\n src tree : ' + str(Curr_src_tree))
      fp.write('\n pruned supertree: ' + str(pruned_tree))
      fp.write('\n src tree len: ' + str(lenSrcTree) + ' pruned supertree len: ' + str(len(pruned_tree.get_edge_set())))
      
      # determine the false positives and the false negatives 
      tup = Curr_src_tree.false_positives_and_negatives(pruned_tree)
      fp.write('   FP_int: ' + str(tup[0]) + '  FN_int:  ' + str(tup[1]))
      sumFP = sumFP + tup[0]
      sumFN = sumFN + tup[1]
      sumRF = sumRF + ((tup[0] + tup[1]) / 2.0)
      
      symm_diff = Curr_src_tree.symmetric_difference(pruned_tree)
      fp.write('   Symmetric difference: ' + str(symm_diff))
      sum_symmetric_diff = sum_symmetric_diff + symm_diff
	  
    # final normalized sumFP's are computed by dividing with the number of trees
    normsumFP = (sumFP * 1.0) / sumLenSrcTree
    normsumFN = (sumFN * 1.0) / sumLenSrcTree
    normsumRF = (sumRF * 1.0) / sumLenSrcTree
    norm_symm_diff = sum_symmetric_diff / (2.0 * sumLenSrcTree)
      
    # print the final result
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
  fp = open(Output_Text_File, 'a')  
  fp.write('\n \n\n ===============>>>>>>>>>>>>>>> TIME COMPLEXITY OF THE METHOD (in seconds) ')
  fp.write('\n \n reading the data: ' + str(data_read_timestamp - start_timestamp) + \
	'\n initialization of the structure: ' + str(data_initialize_timestamp - data_read_timestamp) + \
	'\n formation of the reachability graph (cluster) (after loop): ' + \
	      str(reachability_graph_form_timestamp - data_initialize_timestamp) + \
	'\n multiple parent (related) problem: ' + \
	      str(cluster_of_node_refine_species_timestamp1 - reachability_graph_form_timestamp) + \
	'\n newick string formation: ' + str(newick_str_formation_timestamp - cluster_of_node_refine_species_timestamp1) + \
	  '\n binary tree construction: ' + str(data_process_timestamp - newick_str_formation_timestamp))
	
  fp.write('\n \n Total time taken (in seconds) : ' + str(data_process_timestamp - start_timestamp))
  
  fp.close()
  #--------------------------------------------------------------  
  # delete the storage variables associated with the current execution 
  
  # clear the dictionaries
  Cluster_Info_Dict.clear()
  Taxa_Info_Dict.clear()
  TaxaPair_Reln_Dict.clear()
  
  # clear the lists associated
  if (len(Cost_List_Taxa_Pair_Multi_Reln) > 0):
    Cost_List_Taxa_Pair_Multi_Reln[:] = []
  if (len(Cost_List_Taxa_Pair_Single_Reln) > 0):
    Cost_List_Taxa_Pair_Single_Reln[:] = []
  if (len(COMPLETE_INPUT_TAXA_LIST) > 0):
    COMPLETE_INPUT_TAXA_LIST[:] = []
  if (len(CURRENT_CLUST_IDX_LIST) > 0):
    CURRENT_CLUST_IDX_LIST[:] = []
  
  # free the reachability graph (numpy array)
  del Reachability_Graph_Mat
  
##-----------------------------------------------------

if __name__ == "__main__":
    main() 
  
