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
	#score_val = score_val + TaxaPair_Reln_Dict[key1]._GetConnPrVal(RELATION_R1)
	score_val = score_val + TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R1)
      elif key2 in TaxaPair_Reln_Dict:
	#score_val = score_val + TaxaPair_Reln_Dict[key2]._GetConnPrVal(RELATION_R2)
	score_val = score_val + TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R2)
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
""" 
following code removes extra paranthesis (thus producing insignificant edges) 
from the supertree expression contained in the string Final_Supertree_Str
this is important since the original string (and thus the tree) may contain 
internal nodes of outdegree 1, thus can be removed 
employing the update_splits routine in the dendropy often causes trouble
it is better to manually track the number of taxa and the occurrences of 
enclosing brackets within the string expression
"""
def Remove_Extra_Paranthesis(Final_Supertree_Str):
  # first we convert the input string expression in a list format
  # so as to traverse it in position specific manner
  L = list(Final_Supertree_Str)
  # this is a stack (list structure) which contains the positions of the first brackets encountered
  # within the input string expression
  # this is required to keep track of the correspondence between the opening and the closing brackets
  SL = []
  # this is a dictionary storing the positions of respective opening and closing brackets
  # the position of opening bracket is the dictionary key
  # the value for one key of the dictionary contains the position of the closing bracket
  first_bracket_dict = dict()

  # scan through the tree expression string
  for i in range(len(L)):
    # append the position of an opening bracket in the stack
    if (L[i] == '('):
      SL.append(i)
    # for a closing bracket, retrieve the corresponding opening bracket position
    # and enter those positions in the dictionary
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
	
    # add - sourya
    # if enclosing brackets contain only one species 
    # then those brackets need to be removed
    if (L[i] == '('):
      sb = first_bracket_dict[i]
      f = 0
      for j in range(i,sb):
	if (L[j] == ','):
	  f = 1
	  break
      if (f == 0):
	L.pop(i)
	L.insert(i, '$')
	L.pop(sb)
	L.insert(sb, '$')
    # end add - sourya
  
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
    # comment - sourya
    if (len(outnodes) == 0):
    # add - sourya
    #if (len(outnodes) <= 1):
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
def DefineLeafPairReln(lca_node_rank, node1_rank, node2_rank, xl_val, sum_of_branch_count, \
  node1, node2, edge_type, curr_tree_taxa, taxa_under_curr_node):
  
  key1 = (node1.taxon.label, node2.taxon.label)
  key2 = (node2.taxon.label, node1.taxon.label)
  
  # derive the accumulated rank information
  sum_acc_rank = 0
  if (node1_rank < lca_node_rank):
    sum_acc_rank = sum_acc_rank + (((lca_node_rank - node1_rank) * (lca_node_rank + node1_rank - 1)) / 2)
  if (node2_rank < lca_node_rank):
    sum_acc_rank = sum_acc_rank + (((lca_node_rank - node2_rank) * (lca_node_rank + node2_rank - 1)) / 2)
  # add the rank of the LCA node
  sum_acc_rank = sum_acc_rank + lca_node_rank  
  
  if key1 in TaxaPair_Reln_Dict:
    #intersect_taxa_count = len(set(TaxaPair_Reln_Dict[key1]._GetUnderlyingTaxonList()) & set(curr_tree_taxa)) - 2
    if (len(taxa_under_curr_node) > 0):
      intersect_ratio = (len(set(TaxaPair_Reln_Dict[key1]._GetUnderlyingTaxonList()) & set(curr_tree_taxa)) * 1.0) / len(set(TaxaPair_Reln_Dict[key1]._GetUnderlyingTaxonList()))
    else:
      intersect_ratio = 1
    intersect_ratio1 = 1	#(len(set(taxa_under_curr_node) & set(TaxaPair_Reln_Dict[key1]._GetUnderlyingTaxonList())) * 1.0) / len(set(TaxaPair_Reln_Dict[key1]._GetUnderlyingTaxonList()))    
    TaxaPair_Reln_Dict[key1]._AddSupportingTree()
    TaxaPair_Reln_Dict[key1]._AddLineage(xl_val)
    TaxaPair_Reln_Dict[key1]._AddEdgeCount(edge_type, intersect_ratio)
    TaxaPair_Reln_Dict[key1]._AddLevel(sum_of_branch_count)
    TaxaPair_Reln_Dict[key1]._AddAccumulatedRank(sum_acc_rank)
    TaxaPair_Reln_Dict[key1]._AddLCARank(lca_node_rank)
  elif key2 in TaxaPair_Reln_Dict:
    #intersect_taxa_count = len(set(TaxaPair_Reln_Dict[key2]._GetUnderlyingTaxonList()) & set(curr_tree_taxa)) - 2
    if (len(taxa_under_curr_node) > 0):
      intersect_ratio = (len(set(TaxaPair_Reln_Dict[key2]._GetUnderlyingTaxonList()) & set(curr_tree_taxa)) * 1.0)  / len(set(TaxaPair_Reln_Dict[key2]._GetUnderlyingTaxonList()))
    else:
      intersect_ratio = 1
    intersect_ratio1 = 1	#(len(set(taxa_under_curr_node) & set(TaxaPair_Reln_Dict[key2]._GetUnderlyingTaxonList())) * 1.0) / len(set(TaxaPair_Reln_Dict[key2]._GetUnderlyingTaxonList()))    
    TaxaPair_Reln_Dict[key2]._AddSupportingTree()
    TaxaPair_Reln_Dict[key2]._AddLineage(xl_val)
    TaxaPair_Reln_Dict[key2]._AddEdgeCount(Complementary_Reln(edge_type), intersect_ratio)
    TaxaPair_Reln_Dict[key2]._AddLevel(sum_of_branch_count)
    TaxaPair_Reln_Dict[key2]._AddAccumulatedRank(sum_acc_rank)
    TaxaPair_Reln_Dict[key2]._AddLCARank(lca_node_rank)
  else:
    TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
    #intersect_taxa_count = len(set(TaxaPair_Reln_Dict[key1]._GetUnderlyingTaxonList()) & set(curr_tree_taxa)) - 2
    if (len(taxa_under_curr_node) > 0):
      intersect_ratio = (len(set(TaxaPair_Reln_Dict[key1]._GetUnderlyingTaxonList()) & set(curr_tree_taxa)) * 1.0)  / len(set(TaxaPair_Reln_Dict[key1]._GetUnderlyingTaxonList()))
    else:
      intersect_ratio = 1
    intersect_ratio1 = 1	#(len(set(taxa_under_curr_node) & set(TaxaPair_Reln_Dict[key2]._GetUnderlyingTaxonList())) * 1.0) / len(set(TaxaPair_Reln_Dict[key2]._GetUnderlyingTaxonList()))    
    TaxaPair_Reln_Dict[key1]._AddSupportingTree()
    TaxaPair_Reln_Dict[key1]._AddLineage(xl_val)
    TaxaPair_Reln_Dict[key1]._AddEdgeCount(edge_type, intersect_ratio)
    TaxaPair_Reln_Dict[key1]._AddLevel(sum_of_branch_count)
    TaxaPair_Reln_Dict[key1]._AddAccumulatedRank(sum_acc_rank)
    TaxaPair_Reln_Dict[key1]._AddLCARank(lca_node_rank)
      
  return

#--------------------------------------------------------
# this function derives coupket relations belonging to one tree
# that is provided as an input argument to this function
def DeriveCoupletRelations(Curr_tree, WEIGHT_TAXA_SUBSET):
  
  curr_tree_taxa = Curr_tree.infer_taxa().labels()
  no_of_taxa = len(curr_tree_taxa)  
  
  # traverse the internal nodes of the tree in postorder fashion
  for curr_node in Curr_tree.postorder_internal_node_iter():
    if (WEIGHT_TAXA_SUBSET == True):
      taxa_under_curr_node = GetTaxaUnderInternalNode(curr_node)
    else:
      taxa_under_curr_node = []
    
    # compute the XL value associated with this node    
    # comment - sourya
    # this is the normalized value of extra lineage with respect to the subtree rooted under current node
    # when we divide this absolute value with respect to the number of taxa of the current tree
    ## comment - sourya
    #xl_val = ((len(curr_node.leaf_nodes()) - 2) * 1.0) / no_of_taxa
    # add - sourya
    xl_val = (len(curr_node.leaf_nodes()) - 2)
    curr_node_level = curr_node.level()
    curr_node_rank = no_of_taxa - curr_node_level
    
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
	  node1_rank = no_of_taxa - curr_node_child_leaf_nodes[i].parent_node.level()
	  node2_rank = no_of_taxa - curr_node_child_leaf_nodes[j].parent_node.level()	  
	  sum_of_branch_count = ((curr_node_child_leaf_nodes[i].level() - curr_node_level) + (curr_node_child_leaf_nodes[j].level() - curr_node_level))
	  DefineLeafPairReln(curr_node_rank, node1_rank, node2_rank, xl_val, sum_of_branch_count, \
	    curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j], RELATION_R3, curr_tree_taxa, taxa_under_curr_node)
    
    # one leaf node (direct descendant) and another leaf node (under one internal node)
    # will be related by ancestor / descendant relations
    if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
      for p in curr_node_child_leaf_nodes:
	for q in curr_node_child_internal_nodes:
	  for r in q.leaf_nodes():
	    node1_rank = no_of_taxa - p.parent_node.level()
	    node2_rank = no_of_taxa - r.parent_node.level()	    
	    sum_of_branch_count = ((p.level() - curr_node_level) + (r.level() - curr_node_level))
	    DefineLeafPairReln(curr_node_rank, node1_rank, node2_rank, xl_val, sum_of_branch_count, p, r, RELATION_R1, curr_tree_taxa, taxa_under_curr_node)
    
    # finally a pair of leaf nodes which are descendant of internal nodes will be related by RELATION_R4 relation
    if (len(curr_node_child_internal_nodes) > 1):
      for i in range(len(curr_node_child_internal_nodes) - 1):
	for j in range(i+1, len(curr_node_child_internal_nodes)):
	  for p in curr_node_child_internal_nodes[i].leaf_nodes():
	    for q in curr_node_child_internal_nodes[j].leaf_nodes():
	      node1_rank = no_of_taxa - p.parent_node.level()
	      node2_rank = no_of_taxa - q.parent_node.level()      
	      sum_of_branch_count = ((p.level() - curr_node_level) + (q.level() - curr_node_level))
	      DefineLeafPairReln(curr_node_rank, node1_rank, node2_rank, xl_val, sum_of_branch_count, p, q, RELATION_R4, curr_tree_taxa, taxa_under_curr_node)

#--------------------------------------------------------
# this function derives coupket relations belonging to one tree
# that is provided as an input argument to this function
def FindCoupletUnderlyingTaxon(Curr_tree):
  # traverse the internal nodes of the tree in postorder fashion
  for curr_node in Curr_tree.postorder_internal_node_iter():
    taxa_under_curr_node = GetTaxaUnderInternalNode(curr_node)
    
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
	  node1 = curr_node_child_leaf_nodes[i]
	  node2 = curr_node_child_leaf_nodes[j]
	  key1 = (node1.taxon.label, node2.taxon.label)
	  key2 = (node2.taxon.label, node1.taxon.label)
	  if key1 in TaxaPair_Reln_Dict:
	    TaxaPair_Reln_Dict[key1]._AppendUnderlyingTaxonList(taxa_under_curr_node)
	  elif key2 in TaxaPair_Reln_Dict:
	    TaxaPair_Reln_Dict[key2]._AppendUnderlyingTaxonList(taxa_under_curr_node)
	  else:
	    TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
	    TaxaPair_Reln_Dict[key1]._AppendUnderlyingTaxonList(taxa_under_curr_node)

    
    # one leaf node (direct descendant) and another leaf node (under one internal node)
    # will be related by ancestor / descendant relations
    if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
      for p in curr_node_child_leaf_nodes:
	for q in curr_node_child_internal_nodes:
	  for r in q.leaf_nodes():
	    node1 = p
	    node2 = r
	    key1 = (node1.taxon.label, node2.taxon.label)
	    key2 = (node2.taxon.label, node1.taxon.label)
	    if key1 in TaxaPair_Reln_Dict:
	      TaxaPair_Reln_Dict[key1]._AppendUnderlyingTaxonList(taxa_under_curr_node)
	    elif key2 in TaxaPair_Reln_Dict:
	      TaxaPair_Reln_Dict[key2]._AppendUnderlyingTaxonList(taxa_under_curr_node)
	    else:
	      TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
	      TaxaPair_Reln_Dict[key1]._AppendUnderlyingTaxonList(taxa_under_curr_node)    
    
    # finally a pair of leaf nodes which are descendant of internal nodes will be related by RELATION_R4 relation
    if (len(curr_node_child_internal_nodes) > 1):
      for i in range(len(curr_node_child_internal_nodes) - 1):
	for j in range(i+1, len(curr_node_child_internal_nodes)):
	  for p in curr_node_child_internal_nodes[i].leaf_nodes():
	    for q in curr_node_child_internal_nodes[j].leaf_nodes():
	      node1 = p
	      node2 = q
	      key1 = (node1.taxon.label, node2.taxon.label)
	      key2 = (node2.taxon.label, node1.taxon.label)
	      if key1 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key1]._AppendUnderlyingTaxonList(taxa_under_curr_node)
	      elif key2 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key2]._AppendUnderlyingTaxonList(taxa_under_curr_node)
	      else:
		TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
		TaxaPair_Reln_Dict[key1]._AppendUnderlyingTaxonList(taxa_under_curr_node)

#--------------------------------------------------------
""" this function prints the elements of the queue (which stores the couplet scores 
for individual relations """
def PrintQueueInfo(inp_queue, Output_Text_File):
  fp = open(Output_Text_File, 'a')
  for elem in inp_queue:
    fp.write(' ' + str(elem))
  fp.close()

##-----------------------------------------------------
# this function reads the input tree list file
# parameters: ROOTED_TREE - whether the treelist to be read as rooted format
# PRESERVE_UNDERSCORE: whether underscores of the taxa name will be preserved or not
# INPUT_FILE_FORMAT: data is read from the file according to NEWICK or NEXUS format
# INPUT_FILENAME: file containing the input treelist
def Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME):
  Inp_TreeList = dendropy.TreeList.get_from_path(INPUT_FILENAME, schema=INPUT_FILE_FORMAT, \
						  preserve_underscores=PRESERVE_UNDERSCORE, \
						  default_as_rooted=ROOTED_TREE)
  
  return Inp_TreeList

#--------------------------------------------------
# this function returns the label of an internal or a leaf node 
# in terms of newick representation
def Node_Label(inp_node):
  return str(inp_node.as_newick_string(suppress_edge_lengths=True))

#-----------------------------------------------------
# this is the taxa list generated from current internal node
def GetTaxaUnderInternalNode(curr_node):
  taxa_list_from_curr_internal_node = []
  for n in curr_node.leaf_nodes():
    taxa_list_from_curr_internal_node.append(n.taxon.label)
  return taxa_list_from_curr_internal_node

#----------------------------------------
def Complementary_Reln(inp_reln):
  if (inp_reln == RELATION_R3) or (inp_reln == RELATION_R4):
    return inp_reln
  elif (inp_reln == RELATION_R1):
    return RELATION_R2
  else:
    return RELATION_R1
