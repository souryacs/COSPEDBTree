#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#-----------------------------------------------------
def FindNonProcessedTaxon(Reachability_Graph_Mat, inp_taxon, reln_type):
  reach_mat_inp_taxon_idx = CURRENT_CLUST_IDX_LIST.index(Taxa_Info_Dict[inp_taxon]._Get_Taxa_Part_Clust_Idx())
  ll = list(Taxa_Info_Dict[inp_taxon]._GetRelnList(reln_type))
  for t in ll:
    reach_mat_t_idx = CURRENT_CLUST_IDX_LIST.index(Taxa_Info_Dict[t]._Get_Taxa_Part_Clust_Idx())
    if (Reachability_Graph_Mat[reach_mat_t_idx][reach_mat_inp_taxon_idx] != 0):
      Taxa_Info_Dict[inp_taxon]._RemSelectedReln(reln_type, inp_taxon)
  
  ll = []
  ll.extend(Taxa_Info_Dict[inp_taxon]._GetRelnList(reln_type))
  return ll

##-----------------------------------------------------
"""
this function updates the support score corresponding to the relation 'reln_type'
between the couplet taxonA and taxonB 
by an amount of delta_score
and also updates the queue containing the support score
"""
def UpdateScore(taxonA, taxonB, delta_score, reln_type, Output_Text_File):
  # queue element
  curr_score = TaxaPair_Reln_Dict[(taxonA, taxonB)]._GetEdgeCost_ConnReln(reln_type)
  orig_list_key = [taxonA, taxonB, reln_type, curr_score]
  
  # update the couplet score first
  TaxaPair_Reln_Dict[(taxonA, taxonB)]._IncrEdgeCost_ConnReln(reln_type, delta_score)
  
  fp = open(Output_Text_File, 'a')
  fp.write('\n couplet ' + str(taxonA) + ' and ' + str(taxonB) + ' reln_type : ' + str(reln_type) + ' orig score ' + str(curr_score) + ' curr score ' + str(curr_score + delta_score))
  fp.close()
  
  if (orig_list_key in Cost_List_Taxa_Pair_Single_Reln):
    idx = Cost_List_Taxa_Pair_Single_Reln.index(orig_list_key)
    if (delta_score > 0):
      Queue_Increase_Score(Cost_List_Taxa_Pair_Single_Reln, idx, (curr_score + delta_score))
    else:
      Queue_Decrease_Score(Cost_List_Taxa_Pair_Single_Reln, idx, (curr_score + delta_score))
    
  elif (orig_list_key in Cost_List_Taxa_Pair_Multi_Reln):
    idx = Cost_List_Taxa_Pair_Multi_Reln.index(orig_list_key)
    if (delta_score > 0):
      Queue_Increase_Score(Cost_List_Taxa_Pair_Multi_Reln, idx, (curr_score + delta_score))
    else:
      Queue_Decrease_Score(Cost_List_Taxa_Pair_Multi_Reln, idx, (curr_score + delta_score))

#-----------------------------------------------------
""" 
this function modules the UpdateEdgeCost_Conn_Reln function
arguments: taxonA --- whose neighborhood (x) needs to be reviewed
	   taxonB --- relationship between this node and neighborhood (x) of taxonA needs to be examined
	   neighborhood_reln_type - type of neighborhood of taxonA that is to be examined
	   target_reln_type: species the cost of this type of edge between taxonB and x that needs to be updated
	   node_A_B_reln_type: edge type between taxonA and taxonB
"""
def ModuleUpdCost(taxonA, taxonB, neighborhood_reln_type, target_reln_type, node_A_B_reln_type, Reachability_Graph_Mat, Output_Text_File):    
  taxonA_neighb = FindNonProcessedTaxon(Reachability_Graph_Mat, taxonA, neighborhood_reln_type)
  if taxonB in taxonA_neighb:
    taxonA_neighb.remove(taxonB)	# dont consider taxonB
    
  if (len(taxonA_neighb) > 0):
    for l in taxonA_neighb:    
      # now we calculate the score of the derived edge connection
      # (in terms of deviation from the best possible edge between these 2 nodes) 
      key1 = (taxonB, l)
      key2 = (l, taxonB)
      if key1 in TaxaPair_Reln_Dict:
	delta_score = TaxaPair_Reln_Dict[key1]._GetConnPrVal(target_reln_type)
	if (target_reln_type != RELATION_R4):
	  # check whether originally between these 2 nodes, RELATION_R4 connection was predominant 
	  if (TaxaPair_Reln_Dict[key1]._GetConnPrVal(RELATION_R4) > 0):
	    delta_score = delta_score + ORIG_NO_EDGE_BECOME_CONNECTED         
      elif key2 in TaxaPair_Reln_Dict:
	delta_score = TaxaPair_Reln_Dict[key2]._GetConnPrVal(FindComplementaryReln(target_reln_type))
	if (target_reln_type != RELATION_R4):
	  # check whether originally between these 2 nodes, NO_EDGE connection was predominant 
	  if (TaxaPair_Reln_Dict[key2]._GetConnPrVal(RELATION_R4) > 0):
	    delta_score = delta_score + ORIG_NO_EDGE_BECOME_CONNECTED
      else:
	if (target_reln_type == RELATION_R4):
	  delta_score = ORIG_NO_EDGE_REMAIN_NO_EDGE
	else:
	  delta_score = ORIG_DIFF_SRC_TREE_CONNECT_SCORE

      # update the score of relation between taxonA and l 
      # (where l belongs to the neighborhood of taxonA)
      if (delta_score != 0):
	kt1 = (taxonA, l)
	kt2 = (l, taxonA)
	if kt1 in TaxaPair_Reln_Dict:
	  UpdateScore(taxonA, l, delta_score, neighborhood_reln_type, Output_Text_File)
	elif kt2 in TaxaPair_Reln_Dict:
	  UpdateScore(l, taxonA, delta_score, FindComplementaryReln(neighborhood_reln_type), Output_Text_File)  

#-----------------------------------------------------
"""
this function is called to update the relation score of all non-processed couplets
which are transitively related to a just connected couplet (taxonA, taxonB), 
via the relation 'reln_type'

we check the NON-PROCESSED neighborhood of taxonA and taxonB

1)(say x belongs to non processed neighborhood of taxonB, with ? as the edge type) --- 
  check the possible connection of x to taxonA
  the score of this connection from x to taxonA is the priority associated with corresponding edge type
  this score is added with the score of edge type ? between taxonB to x
2)(say y belongs to non processed neighborhood of taxonA, with ? as the edge type) --- 
  check the possible connection of y to taxonB
  the score of this connection from y to taxonB is the priority associated with corresponding edge type
  this score is added with the score of edge type ? between taxonA to y 
"""
  
def UpdateEdgeCost_Conn_Reln(taxonA, taxonB, reln_type, Reachability_Graph_Mat, Output_Text_File):
  
  # case 1 - A->B connection 
  if (reln_type == RELATION_R1):        
    ##-----------------------------------------
    ## check non processed neighborhood of taxonA
    ##-----------------------------------------
    # check the non processed neighbors of taxonA that can be connected in future via edge C->A
    # check the possible connection C->B 
    ModuleUpdCost(taxonA, taxonB, RELATION_R2, RELATION_R2, reln_type, Reachability_Graph_Mat, Output_Text_File)    
    # check the non processed neighbors of taxonA that can be connected in future via edge C=A
    # check the possible connection C->B
    ModuleUpdCost(taxonA, taxonB, RELATION_R3, RELATION_R2, reln_type, Reachability_Graph_Mat, Output_Text_File)
    # add - sourya
    # check the non processed neighbors of taxonA that can be connected in future via edge C><A
    # check the possible connection C><B
    ModuleUpdCost(taxonA, taxonB, RELATION_R4, RELATION_R4, reln_type, Reachability_Graph_Mat, Output_Text_File)
    # end add - sourya
    ##-----------------------------------------
    ## check non processed neighborhood of taxonB
    ##-----------------------------------------
    # check the non processed neighbors of node B that can be connected in future via edge B->C
    # check the possible connection A->C
    ModuleUpdCost(taxonB, taxonA, RELATION_R1, RELATION_R1, reln_type, Reachability_Graph_Mat, Output_Text_File)    
    # check the non processed neighbors of node B that can be connected in future via edge B=C
    # check the possible connection A->C 
    ModuleUpdCost(taxonB, taxonA, RELATION_R3, RELATION_R1, reln_type, Reachability_Graph_Mat, Output_Text_File)    
      
  # case 2 - A<-B connection 
  if (reln_type == RELATION_R2):    
    ##-----------------------------------------
    ## check non processed neighborhood of taxonA 
    ##-----------------------------------------
    # check the non processed neighbors of node A that can be connected in future via edge A->C
    # check the possible connection B->C
    ModuleUpdCost(taxonA, taxonB, RELATION_R1, RELATION_R1, reln_type, Reachability_Graph_Mat, Output_Text_File)    
    # check the non processed neighbors of node A that can be connected in future via edge A=C
    # check the possible connection B->C 
    ModuleUpdCost(taxonA, taxonB, RELATION_R3, RELATION_R1, reln_type, Reachability_Graph_Mat, Output_Text_File)
    ##-----------------------------------------
    ## check non processed neighborhood of taxonB
    ##-----------------------------------------
    # check the non processed neighbors of node B that can be connected in future via edge C->B
    # establish the possible connection C->A
    ModuleUpdCost(taxonB, taxonA, RELATION_R2, RELATION_R2, reln_type, Reachability_Graph_Mat, Output_Text_File)
    # check the neighbors of node B that can be connected in future via edge C=B
    # establish the connection C->A
    ModuleUpdCost(taxonB, taxonA, RELATION_R3, RELATION_R2, reln_type, Reachability_Graph_Mat, Output_Text_File)
    # add - sourya
    # check the non processed neighbors of taxonB that can be connected in future via edge C><B
    # check the possible connection C><A 
    ModuleUpdCost(taxonB, taxonA, RELATION_R4, RELATION_R4, reln_type, Reachability_Graph_Mat, Output_Text_File)
    # end add - sourya 
    
  # case 3 - A=B connection 
  if (reln_type == RELATION_R3):    
    ##-----------------------------------------
    ## check non processed neighborhood of taxonA
    ##-----------------------------------------
    # check the non processed neighbors of node A that can be connected in future via edge A->C
    # check the possible connection B->C
    ModuleUpdCost(taxonA, taxonB, RELATION_R1, RELATION_R1, reln_type, Reachability_Graph_Mat, Output_Text_File)
    # check the non processed neighbors of node A that can be connected in future via edge C->A 
    # establish the connection C->B
    ModuleUpdCost(taxonA, taxonB, RELATION_R2, RELATION_R2, reln_type, Reachability_Graph_Mat, Output_Text_File)
    # check the non processed neighbors of node A that can be connected in future via edge A=C
    # check the possible connection B=C
    ModuleUpdCost(taxonA, taxonB, RELATION_R3, RELATION_R3, reln_type, Reachability_Graph_Mat, Output_Text_File)
    # add - sourya
    # check the non processed neighbors of node A that can be connected in future via edge A><C
    # check the possible connection B><C
    ModuleUpdCost(taxonA, taxonB, RELATION_R4, RELATION_R4, reln_type, Reachability_Graph_Mat, Output_Text_File)
    # end add - sourya
    ##-----------------------------------------
    ## check non processed neighborhood of taxonB 
    ##-----------------------------------------
    # check the non processed neighbors of node B that can be connected in future via edge B->C
    # check the connection A->C
    ModuleUpdCost(taxonB, taxonA, RELATION_R1, RELATION_R1, reln_type, Reachability_Graph_Mat, Output_Text_File)
    # check the non processed neighbors of node B that can be connected in future via edge C->B 
    # establish the possible connection C->A
    ModuleUpdCost(taxonB, taxonA, RELATION_R2, RELATION_R2, reln_type, Reachability_Graph_Mat, Output_Text_File)    
    # check the non processed neighbors of node B that can be connected in future via edge B=C
    # check the possible connection A=C
    ModuleUpdCost(taxonB, taxonA, RELATION_R3, RELATION_R3, reln_type, Reachability_Graph_Mat, Output_Text_File)    
    # add - sourya
    # check the non processed neighbors of taxonB that can be connected in future via edge B><C
    # check the possible connection A><C  
    ModuleUpdCost(taxonB, taxonA, RELATION_R4, RELATION_R4, reln_type, Reachability_Graph_Mat, Output_Text_File)
    # end add - sourya
		      
  # case 4 - A><B connection (no connection)
  if (reln_type == RELATION_R4):    
    ##-----------------------------------------
    ## check non processed neighborhood of taxonA 
    ##-----------------------------------------
    # check the non processed neighbors of node A that can be connected in future via edge A->C
    # check the possible connection B><C
    ModuleUpdCost(taxonA, taxonB, RELATION_R1, RELATION_R4, reln_type, Reachability_Graph_Mat, Output_Text_File)
    # check the non processed neighbors of node A that can be connected in future via edge A=C
    # check the possible connection B><C      
    ModuleUpdCost(taxonA, taxonB, RELATION_R3, RELATION_R4, reln_type, Reachability_Graph_Mat, Output_Text_File)
    ##-----------------------------------------
    ## check non processed neighborhood of taxonB 
    ##-----------------------------------------
    # check the non processed neighbors of node B that can be connected in future via edge B->C
    # check the possible connection A><C
    ModuleUpdCost(taxonB, taxonA, RELATION_R1, RELATION_R4, reln_type, Reachability_Graph_Mat, Output_Text_File)
    # check the non processed neighbors of node B that can be connected in future via edge B=C
    # check the possible connection A><C
    ModuleUpdCost(taxonB, taxonA, RELATION_R3, RELATION_R4, reln_type, Reachability_Graph_Mat, Output_Text_File)    
  
#------------------------------------------------------
""" this code section implements the max priority queue """
#------------------------------------------------------
# parent node of current node
def Parent(idx):
  return int((idx-1) / 2)

# left child node of current node  
def Left(idx):
  return (2*idx+1)

# right child node of current node
def Right(idx):
  return (2*idx+2)
  
#------------------------------------------------------------------------
# version 4
#-------------------------
"""
def Lower_Score_Value(inp_queue, i, j):
  #---------------------
  # two cases 
  #---------------------
  # case 1 - when both of the scores are non-positive
  if (inp_queue[i][3] <= 0) and (inp_queue[j][3] <= 0):
    if (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
      return 1      
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])):
      return 0      
    elif (inp_queue[i][3] < inp_queue[j][3]):
      return 1
    elif (inp_queue[i][3] > inp_queue[j][3]):
      return 0
    elif ((inp_queue[i][2] == RELATION_R4) and (inp_queue[j][2] != RELATION_R4)):
      # no edge has high priority
      return 0
    elif ((inp_queue[j][2] == RELATION_R4) and (inp_queue[i][2] != RELATION_R4)):
      return 1
    elif ((inp_queue[i][2] == RELATION_R3) and ((inp_queue[j][2] == RELATION_R1) \
      or (inp_queue[j][2] == RELATION_R2))):
      return 1      
    elif ((inp_queue[j][2] == RELATION_R3) and ((inp_queue[i][2] == RELATION_R1) \
      or (inp_queue[i][2] == RELATION_R2))):
      return 0
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
      return 1
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])):
      return 0    
  else:
    # case 2 - otherwise we follow previous algorithm 
    if (inp_queue[i][3] < inp_queue[j][3]):
      return 1
    elif (inp_queue[i][3] == inp_queue[j][3]):
      # for tie case of cost
      if ((inp_queue[i][2] == RELATION_R4) and (inp_queue[j][2] != RELATION_R4)):
	# change - sourya
	# no edge has high priority
	return 0
      elif ((inp_queue[j][2] == RELATION_R4) and (inp_queue[i][2] != RELATION_R4)):
	return 1
      # end change - sourya
      elif ((inp_queue[i][2] == RELATION_R3) and ((inp_queue[j][2] == RELATION_R1) or (inp_queue[j][2] == RELATION_R2))):
	# bi directed edge is set to low priority, compared to other directed edges - add - sourya
	return 1      
      elif ((inp_queue[j][2] == RELATION_R3) and ((inp_queue[i][2] == RELATION_R1) or (inp_queue[i][2] == RELATION_R2))):
	return 0
      elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	    < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
	# higher edge priority have greater priority
	return 1
      elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])\
	    < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])):
	return 0    
      elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	    < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
	# higher edge frequency have greater priority - add - sourya
	return 1      
      elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])\
	    < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])):
	return 0      
      
  return 0
"""
#------------------------------------------------------------------------

def Resolve_R4_and_R4(inp_queue, i, j):
  key_i = (inp_queue[i][0], inp_queue[i][1])
  reln_i = inp_queue[i][2]
  key_j = (inp_queue[j][0], inp_queue[j][1])
  reln_j = inp_queue[j][2]
  
  # first we prioritize the consensus relation
  if (reln_i in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j not in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
    return i
  elif (reln_i not in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
    return j
  # next we prioritize the relation with higher frequency
  elif (TaxaPair_Reln_Dict[key_i]._GetEdgeWeight(reln_i) < TaxaPair_Reln_Dict[key_j]._GetEdgeWeight(reln_j)):
    return j
  elif (TaxaPair_Reln_Dict[key_i]._GetEdgeWeight(reln_i) > TaxaPair_Reln_Dict[key_j]._GetEdgeWeight(reln_j)):
    return i
  # next we prioritize the relation with higher priority
  elif (TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i) < TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)):
    return j
  elif (TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i) > TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)):
    return i  
  # last, we select the relation with high sum of XL
  elif (TaxaPair_Reln_Dict[key_i]._GetAvgTreeXL() < TaxaPair_Reln_Dict[key_j]._GetAvgTreeXL()):
    return j
  else:
    return i

def Resolve_R4_and_R3(inp_queue, i, j):
  key_i = (inp_queue[i][0], inp_queue[i][1])
  reln_i = inp_queue[i][2]
  key_j = (inp_queue[j][0], inp_queue[j][1])
  reln_j = inp_queue[j][2]
  
  # first we prioritize the consensus relation
  if (reln_i in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j not in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
    return i
  elif (reln_i not in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
    return j
  # next we prioritize the relation with higher frequency
  elif (TaxaPair_Reln_Dict[key_i]._GetEdgeWeight(reln_i) < TaxaPair_Reln_Dict[key_j]._GetEdgeWeight(reln_j)):
    return j
  elif (TaxaPair_Reln_Dict[key_i]._GetEdgeWeight(reln_i) > TaxaPair_Reln_Dict[key_j]._GetEdgeWeight(reln_j)):
    return i 
  # next we prioritize the relation with higher priority
  elif (TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i) < TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)):
    return j
  elif (TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i) > TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)):
    return i  
  elif (inp_queue[i][0] == inp_queue[j][0]) and (inp_queue[i][1] == inp_queue[j][1]):
    # same couplets - prioritize R4
    return i
  else:
    # for different couplet, prioritize R3
    return j

def Resolve_R4_and_R1_R2(inp_queue, i, j):
  key_i = (inp_queue[i][0], inp_queue[i][1])
  reln_i = inp_queue[i][2]
  key_j = (inp_queue[j][0], inp_queue[j][1])
  reln_j = inp_queue[j][2]

  # first we prioritize the consensus relation
  if (reln_i in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j not in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
    return i
  elif (reln_i not in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
    return j
  # next we prioritize the relation with higher frequency  
  elif (TaxaPair_Reln_Dict[key_i]._GetEdgeWeight(reln_i) < TaxaPair_Reln_Dict[key_j]._GetEdgeWeight(reln_j)):
    return j
  elif (TaxaPair_Reln_Dict[key_i]._GetEdgeWeight(reln_i) > TaxaPair_Reln_Dict[key_j]._GetEdgeWeight(reln_j)):
    return i
  # next we prioritize the relation with higher priority
  elif (TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i) < TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)):
    return j
  elif (TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i) > TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)):
    return i  
  elif (inp_queue[i][0] == inp_queue[j][0]) and (inp_queue[i][1] == inp_queue[j][1]):
    # same couplets - prioritize R4
    return i
  else:
    ## for different couplets, prioritize R1 / R2 relation
    #return j
    # last, we select the relation with high sum of XL
    if (TaxaPair_Reln_Dict[key_i]._GetAvgTreeXL() < TaxaPair_Reln_Dict[key_j]._GetAvgTreeXL()):
      return j
    else:
      return i
  
def Resolve_R3_and_R1_R2(inp_queue, i, j):
  key_i = (inp_queue[i][0], inp_queue[i][1])
  reln_i = inp_queue[i][2]
  key_j = (inp_queue[j][0], inp_queue[j][1])
  reln_j = inp_queue[j][2]
  
  # first we prioritize the consensus relation
  if (reln_i in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j not in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
    return i
  elif (reln_i not in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
    return j
  # next we prioritize the relation with higher frequency    
  elif (TaxaPair_Reln_Dict[key_i]._GetEdgeWeight(reln_i) < TaxaPair_Reln_Dict[key_j]._GetEdgeWeight(reln_j)):
    return j
  elif (TaxaPair_Reln_Dict[key_i]._GetEdgeWeight(reln_i) > TaxaPair_Reln_Dict[key_j]._GetEdgeWeight(reln_j)):
    return i 
  # next we prioritize the relation with higher priority
  elif (TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i) < TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)):
    return j
  elif (TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i) > TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)):
    return i    
  # R3 is set to low priority
  return j  

def Resolve_R3_and_R3(inp_queue, i, j):
  key_i = (inp_queue[i][0], inp_queue[i][1])
  reln_i = inp_queue[i][2]
  key_j = (inp_queue[j][0], inp_queue[j][1])
  reln_j = inp_queue[j][2]
  
  # first we prioritize the consensus relation
  if (reln_i in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j not in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
    return i
  elif (reln_i not in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
    return j
  # next we prioritize the relation with higher frequency    
  elif (TaxaPair_Reln_Dict[key_i]._GetEdgeWeight(reln_i) < TaxaPair_Reln_Dict[key_j]._GetEdgeWeight(reln_j)):
    return j
  elif (TaxaPair_Reln_Dict[key_i]._GetEdgeWeight(reln_i) > TaxaPair_Reln_Dict[key_j]._GetEdgeWeight(reln_j)):
    return i 
  # next we prioritize the relation with higher priority
  elif (TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i) < TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)):
    return j
  elif (TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i) > TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)):
    return i    
  
  return i

def Resolve_R1_R2_and_R2_R1(inp_queue, i, j):
  key_i = (inp_queue[i][0], inp_queue[i][1])
  reln_i = inp_queue[i][2]
  key_j = (inp_queue[j][0], inp_queue[j][1])
  reln_j = inp_queue[j][2]
  
  # first we prioritize the consensus relation
  if (reln_i in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j not in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
    return i
  elif (reln_i not in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
    return j  
  # next we prioritize the relation with higher frequency    
  elif (TaxaPair_Reln_Dict[key_i]._GetEdgeWeight(reln_i) < TaxaPair_Reln_Dict[key_j]._GetEdgeWeight(reln_j)):
    return j
  elif (TaxaPair_Reln_Dict[key_i]._GetEdgeWeight(reln_i) > TaxaPair_Reln_Dict[key_j]._GetEdgeWeight(reln_j)):
    return i
  # next we prioritize the relation with higher priority
  elif (TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i) < TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)):
    return j
  elif (TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i) > TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)):
    return i    
  # last, we select the relation with high sum of XL
  elif (TaxaPair_Reln_Dict[key_i]._GetAvgTreeXL() < TaxaPair_Reln_Dict[key_j]._GetAvgTreeXL()):
    return j
  
  return i
  

##---------------------------------------------------------

#def Resolve_R4_and_R4(inp_queue, i, j):
  #if (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  #< TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
    #return j
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  #> TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
    #return i
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  #< TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
    #return j
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  #> TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
    #return i  
  ## we select the relation with high sum of XL
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetAvgTreeXL() < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetAvgTreeXL()):
    #return j
  #else:
    #return i

#def Resolve_R4_and_R3(inp_queue, i, j):
  #if (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  #< TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
    #return j
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  #> TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
    #return i  
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  #< TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
    #return j
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  #> TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
    #return i  
  #elif (inp_queue[i][0] == inp_queue[j][0]) and (inp_queue[i][1] == inp_queue[j][1]):
    ## same couplets - prioritize R4
    #return i
  #else:
    ## for different couplet, prioritize R3
    #return j

#def Resolve_R4_and_R1_R2(inp_queue, i, j):
  #if (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  #< TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
    #return j
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  #> TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
    #return i
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  #< TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
    #return j
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  #> TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
    #return i  
  #elif (inp_queue[i][0] == inp_queue[j][0]) and (inp_queue[i][1] == inp_queue[j][1]):
    ## same couplets - prioritize R4
    #return i
  #else:
    ### for different couplets, prioritize R1 / R2 relation
    ##return j
    #if (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetAvgTreeXL() < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetAvgTreeXL()):
      #return j
    #else:
      #return i
  
#def Resolve_R3_and_R1_R2(inp_queue, i, j):
  #if (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  #< TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
    #return j
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  #> TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
    #return i 
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  #< TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
    #return j
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  #> TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
    #return i    
  ## R3 is set to low priority
  #return j  

#def Resolve_R3_and_R3(inp_queue, i, j):
  ## selection between two R3 relations
  ## frequency metric
  #if (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  #< TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
    #return j
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  #< TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
    #return j  
  
  #return i

#def Resolve_R1_R2_and_R2_R1(inp_queue, i, j):
  #if (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  #< TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
    #return j
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  #> TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
    #return i
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  #< TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
    #return j
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  #> TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
    #return i    
  #elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetAvgTreeXL() < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetAvgTreeXL()):
    #return j
  
  #return i



#------------------------------------------------------------------------
# version 1B
# new version by Sourya - 8th July 2015
#-------------------------

def Greater_Score_Value(inp_queue, i, j):
  if (inp_queue[i][3] < inp_queue[j][3]):
    return j
  elif (inp_queue[i][3] > inp_queue[j][3]):
    return i
  else:	#if (inp_queue[i][3] == inp_queue[j][3]):
    # for tie case of cost
    # case A - if both relations are R4
    if ((inp_queue[i][2] == RELATION_R4) and (inp_queue[j][2] == RELATION_R4)):
      return Resolve_R4_and_R4(inp_queue, i, j)
    elif ((inp_queue[i][2] == RELATION_R4) and (inp_queue[j][2] == RELATION_R3)):
      return Resolve_R4_and_R3(inp_queue, i, j)
    elif ((inp_queue[i][2] == RELATION_R3) and (inp_queue[j][2] == RELATION_R4)):
      return Resolve_R4_and_R3(inp_queue, j, i)
    elif ((inp_queue[i][2] == RELATION_R4) and ((inp_queue[j][2] == RELATION_R1) or (inp_queue[j][2] == RELATION_R2))):
      return Resolve_R4_and_R1_R2(inp_queue, i, j)
    elif ((inp_queue[j][2] == RELATION_R4) and ((inp_queue[i][2] == RELATION_R1) or (inp_queue[i][2] == RELATION_R2))):
      return Resolve_R4_and_R1_R2(inp_queue, j, i)
    elif ((inp_queue[i][2] == RELATION_R3) and ((inp_queue[j][2] == RELATION_R1) or (inp_queue[j][2] == RELATION_R2))):
      return Resolve_R3_and_R1_R2(inp_queue, i, j)
    elif ((inp_queue[j][2] == RELATION_R3) and ((inp_queue[i][2] == RELATION_R1) or (inp_queue[i][2] == RELATION_R2))):
      return Resolve_R3_and_R1_R2(inp_queue, j, i)
    elif (inp_queue[i][2] == RELATION_R3) and (inp_queue[j][2] == RELATION_R3):
      return Resolve_R3_and_R3(inp_queue, i, j)
    elif ((inp_queue[i][2] == RELATION_R1) and (inp_queue[j][2] == RELATION_R2)) or ((inp_queue[i][2] == RELATION_R2) and (inp_queue[j][2] == RELATION_R1)): 
      return Resolve_R1_R2_and_R2_R1(inp_queue, i, j)
    
  return i

#------------------------------------------------------------------------
# version 1A
# new version by Sourya - 8th July 2015
#-------------------------
"""
def Lower_Score_Value(inp_queue, i, j):
  if (inp_queue[i][3] < inp_queue[j][3]):
    return 1
  elif (inp_queue[i][3] > inp_queue[j][3]):
    return 0
  else:	#if (inp_queue[i][3] == inp_queue[j][3]):
    # for tie case of cost
    # case A - if both relations are R4
    if ((inp_queue[i][2] == RELATION_R4) and (inp_queue[j][2] == RELATION_R4)):
      # we select the relation with high sum of XL
      if (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetSumXL() < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetSumXL()):
	return 1
      else:
	return 0
    elif ((inp_queue[i][2] == RELATION_R4) and (inp_queue[j][2] != RELATION_R4)):
      # check if these are same couplet
      if (inp_queue[i][0] == inp_queue[j][0]) and (inp_queue[i][1] == inp_queue[j][1]):
	# here no relation will have higher priority
	return 0
      else:
	# no edge has low priority
	return 1
    elif ((inp_queue[j][2] == RELATION_R4) and (inp_queue[i][2] != RELATION_R4)):
      # check if these are same couplet
      if (inp_queue[i][0] == inp_queue[j][0]) and (inp_queue[i][1] == inp_queue[j][1]):
	# here no relation will have higher priority
	return 1
      else:
	# no edge has low priority
	return 0
    elif ((inp_queue[i][2] == RELATION_R3) and ((inp_queue[j][2] == RELATION_R1) or (inp_queue[j][2] == RELATION_R2))):
      # bi directed edge is set to low priority, compared to other directed edges
      return 1      
    elif ((inp_queue[j][2] == RELATION_R3) and ((inp_queue[i][2] == RELATION_R1) or (inp_queue[i][2] == RELATION_R2))):
      return 0
    elif ((inp_queue[i][2] == RELATION_R1) and (inp_queue[j][2] == RELATION_R2)) or ((inp_queue[i][2] == RELATION_R2) and (inp_queue[j][2] == RELATION_R1)): 
      # we select the relation with high sum of XL
      if (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetSumXL() < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetSumXL()):
	return 1
      else:
	return 0
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
      # higher edge priority have greater priority
      return 1
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])):
      return 0    
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
      # higher edge frequency have greater priority - add - sourya
      return 1      
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])):
      return 0      
    
  return 0
"""
#------------------------------------------------------------------------
# version 1 
#-------------------------

# this function is used for cost based sorting of the inp_queue in the unweighted supertree
# this function compares two elements of the heap
# returns 1 if element in the i'th index is lower than element in the j'th index
"""
def Lower_Score_Value(inp_queue, i, j):
  if (inp_queue[i][3] < inp_queue[j][3]):
    return j
  elif (inp_queue[i][3] > inp_queue[j][3]):
    return i
  else:	#if (inp_queue[i][3] == inp_queue[j][3]):
    # for tie case of cost
    if ((inp_queue[i][2] == RELATION_R4) and (inp_queue[j][2] != RELATION_R4)):
      # no edge has low priority
      return j
    elif ((inp_queue[j][2] == RELATION_R4) and (inp_queue[i][2] != RELATION_R4)):	#added
      return i
    elif ((inp_queue[i][2] == RELATION_R3) and ((inp_queue[j][2] == RELATION_R1) or (inp_queue[j][2] == RELATION_R2))):
      # bi directed edge is set to low priority, compared to other directed edges - add - sourya
      return j      
    elif ((inp_queue[j][2] == RELATION_R3) and ((inp_queue[i][2] == RELATION_R1) or (inp_queue[i][2] == RELATION_R2))):	#added
      return i
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
      # higher edge priority have greater priority
      return j
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])):	#added
      return i    
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
      # higher edge frequency have greater priority - add - sourya
      return j      
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])):	#added
      return i      
    
  return i
"""
#------------------------------------------------------------------------
# old version
"""
# this function compares two elements of the heap
# returns 1 if element in the i'th index is lower than element in the j'th index
def Lower_Score_Value(inp_queue, i, j):
  if (inp_queue[i][3] < inp_queue[j][3]):
    return j
  elif (inp_queue[i][3] == inp_queue[j][3]):
    # for tie case of cost
    if ((inp_queue[i][2] == RELATION_R4) and (inp_queue[j][2] != RELATION_R4)):
      # no edge has low priority
      return j
    elif ((inp_queue[i][2] == RELATION_R3) and ((inp_queue[j][2] == RELATION_R1) or (inp_queue[j][2] == RELATION_R2))):	# add - sourya
      # bi directed edge is set to low priority, compared to other directed edges - add - sourya
      return j      
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
      # higher edge priority have greater priority
      return j
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
      # higher edge frequency have greater priority - add - sourya
      return j      
  return i
"""
#------------------------------------------------------------------------
  
# this function exchanges two elements in the heap
def Exchange_Elem(inp_queue, i, j):
  temp_key = inp_queue[i][0]
  inp_queue[i][0] = inp_queue[j][0]
  inp_queue[j][0] = temp_key
  temp_key = inp_queue[i][1]
  inp_queue[i][1] = inp_queue[j][1]
  inp_queue[j][1] = temp_key
  temp_edge = inp_queue[i][2]
  inp_queue[i][2] = inp_queue[j][2]
  inp_queue[j][2] = temp_edge
  temp_val = inp_queue[i][3]
  inp_queue[i][3] = inp_queue[j][3]
  inp_queue[j][3] = temp_val

# maintain max heap property
# note: heap_size may not be the actual length of the queue
# but the working length (on which the remaining sorting operation will take place)
def Max_Heapify(inp_queue, idx, heap_size):
  l = Left(idx)
  r = Right(idx)
  if (l < heap_size) and (Greater_Score_Value(inp_queue, idx, l) == l):
    largest_idx = l
  else:
    largest_idx = idx
  if (r < heap_size) and (Greater_Score_Value(inp_queue, largest_idx, r) == r):
    largest_idx = r
  if (largest_idx != idx):
    Exchange_Elem(inp_queue, idx, largest_idx)
    Max_Heapify(inp_queue, largest_idx, heap_size)

# extract the current maximum and also pop the element from the heap
def Heap_Extract_Max(inp_queue):
  if (len(inp_queue) < 1):
    print 'underflow of max priority queue'
  # 1st element is the maximum
  # comment - sourya
  # this only assigns the pointer which will lead incorrect result
  # when the heap is modified
  # if the 
  #max_elem = inp_queue[0]
  # add - sourya
  max_elem = list(inp_queue[0])
  # end add - sourya
  # replace the first element with the last element of the queue
  Exchange_Elem(inp_queue, 0, len(inp_queue) - 1)
  # delete the last element of the queue
  del inp_queue[len(inp_queue) - 1]
  heap_size = len(inp_queue)
  # call the max_heapify function to maintain the heap property
  # 0 is the starting index of the list storing the heap structure
  Max_Heapify(inp_queue, 0, heap_size)
  return max_elem
  
# this function builds the priority queue (max heap property)
def Build_Max_Heap(inp_queue):
  heap_size = len(inp_queue)
  for idx in range(int(len(inp_queue) / 2), -1, -1):
    Max_Heapify(inp_queue, idx, heap_size)
    
# this is the heap sort algorithm
def Sort_Priority_Queue(inp_queue):
  Build_Max_Heap(inp_queue)
  heap_size = len(inp_queue)
  ## should be commented - sourya
  #for idx in range((len(inp_queue) - 1), 0, -1):
    #Exchange_Elem(inp_queue, 0, idx)
    #heap_size = heap_size - 1
    #Max_Heapify(inp_queue, idx, heap_size)
  
def Queue_Increase_Score(inp_queue, idx, target_val):
  inp_queue[idx][3] = target_val
  while (idx > 0) and (Greater_Score_Value(inp_queue, Parent(idx), idx) == idx):
    Exchange_Elem(inp_queue, idx, Parent(idx))
    idx = Parent(idx)
    
def Queue_Decrease_Score(inp_queue, idx, target_val):
  inp_queue[idx][3] = target_val
  Max_Heapify(inp_queue, idx, len(inp_queue))
      
