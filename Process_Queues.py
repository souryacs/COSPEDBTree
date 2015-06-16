#!/usr/bin/env python

import Header
from Header import *
import Cost_Update
from Cost_Update import *
import ReachGraph_Update
from ReachGraph_Update import *
import UtilFunc
from UtilFunc import *


#-------------------------------------------------------
""" this function checks whether the current relation (given as an input) is feasible to the input supertree configuration
if the current relation is already existing in the supertree, then this function returns 2
if the current relation produces a cycle to the existing configuration of the final supertree, then it returns 1
otherwise the function returns 0 """
def Possible_Conflict_Curr_Reln(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, target_edge_type):
  if (target_edge_type == DIRECTED_IN_EDGE):
    src_taxa_clust_idx = Taxa_Info_Dict[dest_taxa_label]._Get_Taxa_Part_Clust_Idx()
    dest_taxa_clust_idx = Taxa_Info_Dict[src_taxa_label]._Get_Taxa_Part_Clust_Idx()
  else:
    src_taxa_clust_idx = Taxa_Info_Dict[src_taxa_label]._Get_Taxa_Part_Clust_Idx()
    dest_taxa_clust_idx = Taxa_Info_Dict[dest_taxa_label]._Get_Taxa_Part_Clust_Idx()

  #--------------------------------------------------------------------
  # check whether the taxa pairs are already resolved
  
  # if the taxa pair is already part of the same cluster
  # that is, they are already related via equivalence partition 
  # then return
  if (src_taxa_clust_idx == dest_taxa_clust_idx):
    return 1
    
  # we find the indices of the reachability matrix corresponding to individual cluster indices 
  reach_mat_src_taxa_clust_idx = CURRENT_CLUST_IDX_LIST.index(src_taxa_clust_idx)
  reach_mat_dest_taxa_clust_idx = CURRENT_CLUST_IDX_LIST.index(dest_taxa_clust_idx)
  
  # case A - if the clusters are already related (depicted in the reachability matrix) then return 1
  if (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][reach_mat_src_taxa_clust_idx] > 0) or \
      (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][reach_mat_dest_taxa_clust_idx] > 0):
    return 1

  #--------------------------------------------------------------------
  # check for the transitive conflict
  if (target_edge_type == DIRECTED_IN_EDGE) or (target_edge_type == DIRECTED_OUT_EDGE):
    # if A->B is to be established
    # and there exists D->A 
    # then if B->D or B><D then return a conflict
    for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetInEdgeList():
      if (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] > 0):
	return 1
    # if A->B is to be established
    # and there exists B->E 
    # then if E->A or E><A then return a conflict  
    for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
      if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] > 0):
	return 1
    # if A->B is to be established
    # and there exists D->A and B->E 
    # then if E->D or E=D or E><D then return a conflict  
    for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetInEdgeList():
      for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
	if (x == y) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] > 0):
	  return 1
	      
    # if A->B is to be established
    # and there exists D><A
    # then it would be D><B afterwards 
    # so if there is D->B, B->D or D=B then return a conflict  
    for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetNoEdgeList():
      if (dest_taxa_clust_idx == x) \
	or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] == 1) \
	or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1):
	return 1
    
    # if A->B is to be established
    # and there exists D><A
    # then for all B->E
    # it would be D><E afterwards 
    # so if D->E or E->D or D=E then return a conflict  
    for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetNoEdgeList():
      for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
	if (x == y) \
	  or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
	  or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 1):
	  return 1

  elif (target_edge_type == BI_DIRECTED_EDGE):
    # if A=B is to be established
    # and there exists B->D 
    # then if D->A or D><A or D=A then return a conflict    
    for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
      if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] > 0):
	return 1
    # if A=B is to be established
    # and there exists D->B 
    # then if A->D or A><D or A=D then return a conflict    	
    for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetInEdgeList():
      if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] > 0):
	return 1
    
    # if A=B is to be established
    # and there exists B><D 
    # then if A->D or D->A or A=D then return a conflict    	
    for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetNoEdgeList():
      if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
	  or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] == 1):
	return 1
	
    # if A=B is to be established
    # and there exists A->D 
    # then if D->B or D><B or D=B then return a conflict    	
    for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetOutEdgeList():
      if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] > 0):
	return 1
    # if A=B is to be established
    # and there exists D->A 
    # then if B->D or B><D or B=D then return a conflict    	
    for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetInEdgeList():
      if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] > 0):
	return 1
    
    # if A=B is to be established
    # and there exists A><D 
    # then if D->B or B->D or B=D then return a conflict    	
    for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetNoEdgeList():
      if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
	  or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] == 1):
	return 1
	
  else:	#if (target_edge_type == NO_EDGE):
    # construct the out neighborhood of src_cluster
    # it will contain the cluster itself and all the other clusters connected via out edges from this cluster
    src_clust_out_neighb = []
    src_clust_out_neighb.append(src_taxa_clust_idx)
    src_clust_out_neighb.extend(Cluster_Info_Dict[src_taxa_clust_idx]._GetOutEdgeList())
    # construct the out neighborhood of dest cluster
    # it will contain the cluster itself and all the other clusters connected via out edges from this cluster    
    dest_clust_out_neighb = []
    dest_clust_out_neighb.append(dest_taxa_clust_idx)
    dest_clust_out_neighb.extend(Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList())
    
    # if mutually any pair of taxa belonging to respective out edge clusters
    # are themselves related via any relationships other than NO EDGE then this connection is not possible
    for x in src_clust_out_neighb:
      for y in dest_clust_out_neighb:
	if (x == y) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
	  or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 1):
	  return 1
    
  return 0
  
#-------------------------------------------------------
""" this function processes individual queues (score lists)
to construct the final supertree 
Single_Reln_no_Conflict_Queue_process: if 1, then points to the single relation queue of score metrics
	      otherwise, points to the multi relation queue of score metrics
it can be collection of taxa pairs exhibiting single relation instance
or can be taxa pairs exhibiting multi relation instance """
def Proc_Queue(Reachability_Graph_Mat, Single_Reln_no_Conflict_Queue_process, Output_Text_File):
  if (Single_Reln_no_Conflict_Queue_process == 1):
    Inp_Queue = Cost_List_Taxa_Pair_Single_Reln
  else:
    Inp_Queue = Cost_List_Taxa_Pair_Multi_Reln
  
  while (0 < len(Inp_Queue)):
    """ extract the 1st element of "Inp_Queue" 
    since it is sorted to have max cost at the beginning """
    outlist = Heap_Extract_Max(Inp_Queue)
    
    src_taxa_label = outlist[0]
    dest_taxa_label = outlist[1]
    src_to_dest_edge_type = outlist[2]
    conn_score = TaxaPair_Reln_Dict[(src_taxa_label, dest_taxa_label)]._GetEdgeCost_ConnReln(src_to_dest_edge_type)
    
    if (DEBUG_LEVEL >= 2):
      fp = open(Output_Text_File, 'a')
      fp.write('\n ===>> current extracted max element: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
	' edge type: ' + str(src_to_dest_edge_type) + ' conn score: ' + str(conn_score))
      fp.close()
    
    """ if the current score metric based relation does not induce a cycle to the existing configuration of the final supertree
    then include the connection in it """
    conflict_detection = Possible_Conflict_Curr_Reln(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, src_to_dest_edge_type)
    
    if (conflict_detection == 0):      
      """ current element does not create a cycle / conflict
      not that it is already present in the supertree 
      valid connection is found - append this connection to the final formed tree """
      if (DEBUG_LEVEL > 0):
	fp = open(Output_Text_File, 'a')    	
	fp.write('\n ==>>>>>>>>> NEW CONN --- nodes to be connected: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
	      ' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) + ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
	      ' edge type: ' + str(src_to_dest_edge_type) + ' conn score: ' + str(conn_score))
	fp.close()
      
      # also update the reachability graph information
      Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, src_to_dest_edge_type)
                  
  return Reachability_Graph_Mat
