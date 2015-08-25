#!/usr/bin/env python

import Header
from Header import *
import Cost_Update
from Cost_Update import *
import ReachGraph_Update
from ReachGraph_Update import *
import Conflict_Detect
from Conflict_Detect import *
import UtilFunc
from UtilFunc import *
  
#-------------------------------------------------------
""" this function processes individual queues (score lists) to construct the final supertree 
    Single_Reln_no_Conflict_Queue_process: if 1, then points to the single relation queue of score metrics
	      otherwise, points to the multi relation queue of score metrics
    it can be collection of taxa pairs exhibiting single relation instance
    or can be taxa pairs exhibiting multi relation instance """
def Proc_Queue(Reachability_Graph_Mat, Single_Reln_no_Conflict_Queue_process, Output_Text_File, DYNAMIC_SCORE_UPDATE):
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
      if (Single_Reln_no_Conflict_Queue_process == 1):
	fp.write('\n ===>> NON CONFLICTING QUEUE -- ')
      else:
	fp.write('\n ===>> CONFLICTING QUEUE -- ')      
      fp.write(' current extracted max element: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
	' edge type: ' + str(src_to_dest_edge_type) + ' conn score: ' + str(conn_score))
      fp.close()
    
    """ if the current score metric based relation does not induce a cycle to the existing configuration of the final supertree
    then include the connection in it """
    conflict_detection = Possible_Conflict_Curr_Reln(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, src_to_dest_edge_type, Output_Text_File)
    
    if (conflict_detection == 0):      
      """ current element does not create a cycle / conflict
      not that it is already present in the supertree 
      valid connection is found - append this connection to the final formed tree """
      if (DEBUG_LEVEL > 0):
	fp = open(Output_Text_File, 'a')    
	if (Single_Reln_no_Conflict_Queue_process == 1):
	  fp.write('\n ==>>>>>>>>> NEW CONN --- NON CONFLICTING QUEUE -- ')
	else:
	  fp.write('\n ==>>>>>>>>> NEW CONN --- CONFLICTING QUEUE -- ')
	  
	fp.write(' nodes to be connected: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
	      ' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) + ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
	      ' edge type: ' + str(src_to_dest_edge_type) + ' conn score: ' + str(conn_score))
	fp.close()
      
      # also update the reachability graph information
      Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, src_to_dest_edge_type, DYNAMIC_SCORE_UPDATE)
      
      # this function is called to update the edge score between the two nodes, and also on their neighborhood
      # only the non processed neighborhood is considered 
      if (DYNAMIC_SCORE_UPDATE == True):
	for list_idx in range(len(COUPLET_CONNECTED_LIST)):
	  taxonA = COUPLET_CONNECTED_LIST[list_idx][0]
	  taxonB = COUPLET_CONNECTED_LIST[list_idx][1]
	  reln_type = COUPLET_CONNECTED_LIST[list_idx][2]
	  UpdateEdgeCost_Conn_Reln(taxonA, taxonB, reln_type, Reachability_Graph_Mat, Output_Text_File)
	
	# now reset the COUPLET_CONNECTED_LIST structure
	del COUPLET_CONNECTED_LIST[:]
      
    else:
      # conflict is detected
      # remove relation information from individual taxon structure
      Taxa_Info_Dict[src_taxa_label]._RemSelectedReln(src_to_dest_edge_type, dest_taxa_label)
      Taxa_Info_Dict[dest_taxa_label]._RemSelectedReln(FindComplementaryReln(src_to_dest_edge_type), src_taxa_label)
                  
  return Reachability_Graph_Mat
