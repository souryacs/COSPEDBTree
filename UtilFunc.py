#!/usr/bin/env python

import Header
from Header import *
    
##-----------------------------------------------------
# this function computes the score (ancestor relation) from clust1 to clust2
def ComputeScore(clust1, clust2, Output_Text_File):
  
  if (DEBUG_LEVEL > 2):
    fp = open(Output_Text_File, 'a')
    fp.write('\n ===>>> Score --- cluster 1 - taxa set: ' + str(Cluster_Info_Dict[clust1]._GetSpeciesList())\
      + ' cluster 2 taxa set ' + str(Cluster_Info_Dict[clust2]._GetSpeciesList()))
    
  score_val = 0
  for t1 in Cluster_Info_Dict[clust1]._GetSpeciesList():
    for t2 in Cluster_Info_Dict[clust2]._GetSpeciesList():
      key1 = (t1, t2)
      key2 = (t2, t1)
      if key1 in TaxaPair_Reln_Dict:
	score_val = score_val + TaxaPair_Reln_Dict[key1]._GetConnPrVal(DIRECTED_OUT_EDGE)
      elif key2 in TaxaPair_Reln_Dict:
	score_val = score_val + TaxaPair_Reln_Dict[key2]._GetConnPrVal(DIRECTED_IN_EDGE)
      else:
	if (DEBUG_LEVEL > 2):
	  fp.write('\n score compute -- key pair ' + str(t1) + ',' + str(t2) + ' does not exist ')

  if (DEBUG_LEVEL > 2):
    fp.write('\n pairwise score of this cluster pair is : ' + str(score_val))
    fp.close()
    
  return score_val

##-----------------------------------------------------    
""" this function solves multiple parenting problem (C2)
by uniquely selecting one particular parent
the selection is carried out using a scoring mechanism """
def SolveMultipleParentC2Problem(Output_Text_File):
  for cx in Cluster_Info_Dict:
    
    if (DEBUG_LEVEL > 2):
      fp = open(Output_Text_File, 'a')
      fp.write('\n ***** Examining cluster -- ')
      fp.close()      
      Cluster_Info_Dict[cx]._PrintClusterInfo(cx, Output_Text_File)    
    
    if (Cluster_Info_Dict[cx]._Get_Indegree() > 1):
      # at first form the list to contain the score values for all the child nodes
      # we'll define the score value later
      if (DEBUG_LEVEL > 2):
	fp = open(Output_Text_File, 'a')
	fp.write('\n ***** Examining cluster with more than one indegree -- before in edge list fixing: ')
	fp.close()      
	Cluster_Info_Dict[cx]._PrintClusterInfo(cx, Output_Text_File)
      
      # initialize one dictionary keyed by cluster indices
      scoring_dict = dict()
      for cz in Cluster_Info_Dict[cx]._GetInEdgeList():
	scoring_dict.setdefault(cz, 0)
      # now for each of the in clusters, compute the score 
      for cz in Cluster_Info_Dict[cx]._GetInEdgeList():
	scoring_dict[cz] = ComputeScore(cz, cx, Output_Text_File)
      
      # open the output file
      if (DEBUG_LEVEL > 2):
	fp = open(Output_Text_File, 'a')      
      
      # after computing all such scores for all the in edge clusters
      # we store the values in a list and sort it
      Scoring_List = []
      for cz in Cluster_Info_Dict[cx]._GetInEdgeList():
	if (DEBUG_LEVEL > 2):
	  fp.write('\n scoring dict elem: ' + str(cz) + ' score: ' + str(scoring_dict[cz]))
	temp_subl = [cz, scoring_dict[cz]]
	Scoring_List.append(temp_subl)
      
      # after obtaining scores for different taxa set belonging under individual child 
      # nodes of the current node under study, we decide about their ancestral / descendant relationships
      if (DEBUG_LEVEL > 2):
	fp.write('\n --- before sorting the scoring list --- ')
	for i in range(len(Scoring_List)):
	  fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(Scoring_List[i][0]) + ' score: ' + str(Scoring_List[i][1]))
	
      # sort the scoring list
      Scoring_List.sort(key=lambda x: x[1])
      if (DEBUG_LEVEL > 2):
	fp.write('\n --- after sorting the scoring list --- ')
	for i in range(len(Scoring_List)):
	  fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(Scoring_List[i][0]) + ' score: ' + str(Scoring_List[i][1]))
   
      # close the output file
      if (DEBUG_LEVEL > 2):
	fp.close()
	
      # after sorting the scoring list, now remove all except the last element from the 
      # cx cluster in edge lists
      for i in range(len(Scoring_List) - 1):
	target_delete_clust_idx = Scoring_List[i][0]
	Cluster_Info_Dict[cx]._RemoveInEdge(target_delete_clust_idx)
	Cluster_Info_Dict[target_delete_clust_idx]._RemoveOutEdge(cx)
    
      if (DEBUG_LEVEL > 2):
	fp = open(Output_Text_File, 'a')
	fp.write('\n ***** Examining cluster with more than one indegree -- after in edge list fixing: ')
	fp.close()      
	Cluster_Info_Dict[cx]._PrintClusterInfo(cx, Output_Text_File)
    
##-----------------------------------------------------    
# following code removes extra paranthesis (thus producing insignificant edges) 
# from the tree expression contained in the string
def Remove_Extra_Paranthesis(Final_Supertree_Str):
  L = list(Final_Supertree_Str)
  SL = []	# stack list
  first_bracket_dict = dict()
  
  for i in range(len(L)):
    if (L[i] == '('):
      SL.append(i)
    elif (L[i] == ')'):
      first_bracket_idx = SL.pop()
      first_bracket_dict.setdefault(first_bracket_idx, i)
 
  if 0:	#(DEBUG_LEVEL > 2):
    print 'L : ', L
    print 'first_bracket_dict: ', first_bracket_dict
  
  for i in range(len(L) - 1):
    if (L[i] == '(') and (L[i+1] == '('):
      sb1 = first_bracket_dict[i]
      sb2 = first_bracket_dict[i+1]
      if (sb1 - sb2 == 1):
	L.pop(i)
	L.insert(i, '$')
	L.pop(sb1)
	L.insert(sb1, '$')
  
  if 0:	#(DEBUG_LEVEL > 2):
    print 'L : ', L
  
  while (1):
    if '$' in L:
      L.remove('$')
    else:
      break
  
  if 0:	#(DEBUG_LEVEL > 2):
    print 'L : ', L

  # construct the string containing final supertree
  outstr = ''.join(L)
  return outstr
    
##-----------------------------------------------------        
''' this function returns the root node for the final supertree 
for a depth first forest, multiple root nodes can be possible - so it returns the node with 0 indegree '''
def Extract_Node_Min_Indeg(no_of_clusters):
  min_indeg_node_idx = -1
  valid_node_found = 0
  for i in Cluster_Info_Dict:
    if (Cluster_Info_Dict[i]._GetExploredStatus() == 0):
      if (valid_node_found == 0):
	min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
	min_indeg_node_idx = i
	valid_node_found = 1
      elif (valid_node_found == 1) and (Cluster_Info_Dict[i]._Get_Indegree() < min_indeg):
	min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
	min_indeg_node_idx = i
      elif (valid_node_found == 1) and (Cluster_Info_Dict[i]._Get_Indegree() == min_indeg)\
	    and (Cluster_Info_Dict[i]._Get_Outdegree() > Cluster_Info_Dict[min_indeg_node_idx]._Get_Outdegree()):    
	min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
	min_indeg_node_idx = i
    
  return min_indeg_node_idx
  
##-----------------------------------------------------
# this function prints the tree in Newick format
def PrintNewick(root_clust_node_idx):
  if 0:
    print 'in function printnewick:   root_clust_node_idx: ', root_clust_node_idx
    print 'taxa set: ', Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()  
    print 'out clust list: ', Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList()
  
  Tree_Str_List = ''
  # process the node provided it has not been explored yet
  if (Cluster_Info_Dict[root_clust_node_idx]._GetExploredStatus() == 0):  
    # set the explored status of the current node to true
    Cluster_Info_Dict[root_clust_node_idx]._SetExploredStatus()
    # get the out edge list of the current node which are not explored yet 
    outnodes = []
    for l in Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList():
      if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
	outnodes.append(l)
    if (len(outnodes) == 0):
      spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
      if (len(spec_list) > 1):
	Tree_Str_List = Tree_Str_List + '('
      Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
      if (len(spec_list) > 1):
	Tree_Str_List = Tree_Str_List + ')'
    else:
      Tree_Str_List = Tree_Str_List + '('
      Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList())
      Tree_Str_List = Tree_Str_List + ','    
      Tree_Str_List = Tree_Str_List + '('
      for i in range(len(outnodes)):
	if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):  
	  Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i])
	  if (i < (len(outnodes) - 1)):
	    # we check whether any subsequent node belonging to the outnodes list
	    # is left for traverse
	    j = i + 1
	    while (j < len(outnodes)):
	      if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
		break
	      j = j + 1	      
	    # in this case, we append one comma
	    if (j < len(outnodes)):
	      Tree_Str_List = Tree_Str_List + ','
      
      Tree_Str_List = Tree_Str_List + ')'
      Tree_Str_List = Tree_Str_List + ')'
    
  return Tree_Str_List    
    
##-----------------------------------------------------  
''' this function performs transitive reduction of a graph (transitive closure) and subsequently modifies the cluster of nodes
in terms of the edge connectivity, to make it free of redunant edges '''
def CompressDirectedGraph(Reachability_Graph_Mat):  
  no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
  # transitive reduction
  for j in range(no_of_clusters):
    for i in range(no_of_clusters):
      if (Reachability_Graph_Mat[i][j] == 1):
	for k in range(no_of_clusters):
	  if (Reachability_Graph_Mat[j][k] == 1) and (Reachability_Graph_Mat[i][k] == 1):
	    Reachability_Graph_Mat[i][k] = 0
	    # remove the edge from the cluster node directory
	    clust_i = CURRENT_CLUST_IDX_LIST[i]
	    clust_k = CURRENT_CLUST_IDX_LIST[k]
	    Cluster_Info_Dict[clust_i]._RemoveOutEdge(clust_k)
	    Cluster_Info_Dict[clust_k]._RemoveInEdge(clust_i)
    
##-----------------------------------------------------
""" this function creates one new cluster with the given index value
also, it inserts one specified taxa in that cluster """
def Create_Cluster_Taxa_Label(target_clust_idx, target_taxa_label):
  # create the cluster
  Cluster_Info_Dict.setdefault(target_clust_idx, Cluster_node(target_taxa_label))
  # include the cluster idx in the global list CURRENT_CLUST_IDX_LIST
  CURRENT_CLUST_IDX_LIST.append(target_clust_idx)
  # mention the cluster index in the taxa information
  Taxa_Info_Dict[target_taxa_label]._Set_Clust_Idx_taxa_Part(target_clust_idx)
  
##-----------------------------------------------------
""" this function appends one specified taxon on a given cluster """
def Append_Cluster_Taxa_Label(target_clust_idx, target_taxa_label):
  if target_taxa_label not in Cluster_Info_Dict[target_clust_idx]._GetSpeciesList():
    Cluster_Info_Dict[target_clust_idx]._Append_taxa(target_taxa_label)
    # mention the cluster index in the taxa information
    Taxa_Info_Dict[target_taxa_label]._Set_Clust_Idx_taxa_Part(target_clust_idx)      
        
#--------------------------------------------------------
# this function defines relationship between a pair of nodes in a tree
# the relationship is either ancestor / descendant, or siblings, or no relationship 
def DefineLeafPairReln(xl_val, node1, node2, edge_type):
  key1 = (node1.taxon.label, node2.taxon.label)
  key2 = (node2.taxon.label, node1.taxon.label)
  if key1 in TaxaPair_Reln_Dict:
    TaxaPair_Reln_Dict[key1]._AddSupportingTree()
    TaxaPair_Reln_Dict[key1]._AddLineage(xl_val)
    TaxaPair_Reln_Dict[key1]._AddEdgeCount(edge_type)
  elif key2 in TaxaPair_Reln_Dict:
    TaxaPair_Reln_Dict[key2]._AddSupportingTree()
    TaxaPair_Reln_Dict[key2]._AddLineage(xl_val)
    if (edge_type == BI_DIRECTED_EDGE) or (edge_type == NO_EDGE):
      TaxaPair_Reln_Dict[key2]._AddEdgeCount(edge_type)
    elif (edge_type == DIRECTED_OUT_EDGE):
      TaxaPair_Reln_Dict[key2]._AddEdgeCount(DIRECTED_IN_EDGE)
    else:
      TaxaPair_Reln_Dict[key2]._AddEdgeCount(DIRECTED_OUT_EDGE)
  else:
    TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
    TaxaPair_Reln_Dict[key1]._AddSupportingTree()
    TaxaPair_Reln_Dict[key1]._AddLineage(xl_val)
    TaxaPair_Reln_Dict[key1]._AddEdgeCount(edge_type)
      
  return

#--------------------------------------------------------
# this function derives coupket relations belonging to one tree
# that is provided as an input argument to this function
def DeriveCoupletRelations(Curr_tree):
  # traverse the internal nodes of the tree in postorder fashion
  for curr_node in Curr_tree.postorder_internal_node_iter():
    # compute the XL value associated with this node
    xl_val = (len(curr_node.leaf_nodes()) - 2)
    
    # list the leaf and internal children of the current node
    curr_node_child_leaf_nodes = []
    curr_node_child_internal_nodes = []
    for x in curr_node.child_nodes():
      if (x.is_leaf() == True):
	curr_node_child_leaf_nodes.append(x)
      else:
	curr_node_child_internal_nodes.append(x)
    
    # pair of leaf nodes will be related by sibling relations
    if (len(curr_node_child_leaf_nodes) > 1):
      for i in range(len(curr_node_child_leaf_nodes) - 1):
	for j in range(i+1, len(curr_node_child_leaf_nodes)):
	  DefineLeafPairReln(xl_val, curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j], BI_DIRECTED_EDGE)
    
    # one leaf node (direct descendant) and another leaf node (under one internal node)
    # will be related by ancestor / descendant relations
    if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
      for p in curr_node_child_leaf_nodes:
	for q in curr_node_child_internal_nodes:
	  for r in q.leaf_nodes():
	    DefineLeafPairReln(xl_val, p, r, DIRECTED_OUT_EDGE)
    
    # finally a pair of leaf nodes which are descendant of internal nodes will be related by NO_EDGE relation
    if (len(curr_node_child_internal_nodes) > 1):
      for i in range(len(curr_node_child_internal_nodes) - 1):
	for j in range(i+1, len(curr_node_child_internal_nodes)):
	  for p in curr_node_child_internal_nodes[i].leaf_nodes():
	    for q in curr_node_child_internal_nodes[j].leaf_nodes():
	      DefineLeafPairReln(xl_val, p, q, NO_EDGE)


#--------------------------------------------------------
""" this function prints the elements of the queue (which stores the couplet scores 
for individual relations """
def PrintQueueInfo(inp_queue, Output_Text_File):
  fp = open(Output_Text_File, 'a')
  for elem in inp_queue:
    fp.write(' ' + str(elem))
  fp.close()
