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

##-------------------------------------------------------
#"""
#this function processes couplets which have single R4 relation
#"""
#def Process_Include_R4_Reln(Reachability_Graph_Mat, Output_Text_File):
	#for l in TaxaPair_Reln_Dict:
		#if (TaxaPair_Reln_Dict[l]._Check_Single_Reln_OnlyAllowed(RELATION_R4) == True):
			##or (TaxaPair_Reln_Dict[l]._CheckTargetRelnMajorityConsensus(RELATION_R4) == True):
			#src_taxa_idx = l[0]
			#src_taxa_label = COMPLETE_INPUT_TAXA_LIST[src_taxa_idx]
			#dest_taxa_idx = l[1]
			#dest_taxa_label = COMPLETE_INPUT_TAXA_LIST[dest_taxa_idx]
			#conn_score = TaxaPair_Reln_Dict[(src_taxa_idx, dest_taxa_idx)]._GetEdgeCost_ConnReln(RELATION_R4)
			
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n ===>> NEW CONN --- CONFLICTING QUEUE -- ' + str(src_taxa_idx) + '(' + str(src_taxa_label) + \
					#') and ' + str(dest_taxa_idx) + '(' + str(dest_taxa_label) + ')' + \
						#' relation type: ' + str(RELATION_R4) + ' conn score: ' + str(conn_score))
				#fp.close()
			
			## also update the reachability graph information
			#Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R4)
	
	#return Reachability_Graph_Mat

#-------------------------------------------------------
""" 
this function processes individual queues (score lists) to construct the final supertree 
    bool_nonconflict_queue: if 1, then points to the single relation queue of score metrics
	      otherwise, points to the multi relation queue of score metrics
    it can be collection of taxa pairs exhibiting single relation instance
    or can be taxa pairs exhibiting multi relation instance 
"""
def Proc_Queue(Reachability_Graph_Mat, bool_nonconflict_queue, Output_Text_File):
	if (bool_nonconflict_queue == 1):
		Inp_Queue = Cost_List_Taxa_Pair_Single_Reln
	else:
		Inp_Queue = Cost_List_Taxa_Pair_Multi_Reln

	while (0 < len(Inp_Queue)):
		""" 
		extract the 1st element of "Inp_Queue" 
		since it is sorted to have max cost at the beginning 
		"""
		outlist = Heap_Extract_Max(Inp_Queue)
		
		src_taxa_idx = outlist[0]
		src_taxa_label = COMPLETE_INPUT_TAXA_LIST[src_taxa_idx]
		dest_taxa_idx = outlist[1]
		dest_taxa_label = COMPLETE_INPUT_TAXA_LIST[dest_taxa_idx]
		reln_type = outlist[2]
		conn_score = TaxaPair_Reln_Dict[(src_taxa_idx, dest_taxa_idx)]._GetEdgeCost_ConnReln(reln_type)

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			if (bool_nonconflict_queue == 1):
				fp.write('\n ===>> NON CONFLICTING QUEUE -- ')
			else:
				fp.write('\n ===>> CONFLICTING QUEUE -- ')      
			fp.write(' current extracted max element: ' + str(src_taxa_idx) + '(' + str(src_taxa_label) + \
				') and ' + str(dest_taxa_idx) + '(' + str(dest_taxa_label) + ')' + \
								' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
			fp.close()

		""" 
		if the current score metric based relation does not induce a cycle to the existing configuration of the final supertree
		then include the connection in it 
		"""
		conflict_detection = Possible_Conflict_Curr_Reln(src_taxa_idx, dest_taxa_idx, Reachability_Graph_Mat, reln_type, Output_Text_File)
		
		if (conflict_detection == 0):   
			""" 
			current element does not create a cycle / conflict
			not that it is already present in the supertree 
			valid connection is found - append this connection to the final formed tree 
			"""
			if (DEBUG_LEVEL > 0):
				fp = open(Output_Text_File, 'a')    
				if (bool_nonconflict_queue == 1):
					fp.write('\n ==>>>>>>>>> NEW CONN --- NON CONFLICTING QUEUE -- ')
				else:
					fp.write('\n ==>>>>>>>>> NEW CONN --- CONFLICTING QUEUE -- ')
				fp.write(' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
				fp.close()
			
			# also update the reachability graph information
			Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, reln_type)
			
	return Reachability_Graph_Mat
