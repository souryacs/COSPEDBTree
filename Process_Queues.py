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
""" 
this function processes individual queues (score lists) to construct the final supertree 
    bool_nonconflict_queue: if 1, then points to the single relation queue of score metrics
	      otherwise, points to the multi relation queue of score metrics
    it can be collection of taxa pairs exhibiting single relation instance
    or can be taxa pairs exhibiting multi relation instance 
"""
def Proc_Queue(Reachability_Graph_Mat, bool_nonconflict_queue, Output_Text_File):
	if (bool_nonconflict_queue == 2):
		Inp_Queue = Queue_Score_Conflict_Couplet
	elif (bool_nonconflict_queue == 1):
		Inp_Queue = Queue_Score_2Reln_2Tree
	else:
		Inp_Queue = Queue_Score_1Reln_2Tree

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
		ckey = (src_taxa_idx, dest_taxa_idx)
		conn_score = TaxaPair_Reln_Dict[ckey]._GetEdgeCost_ConnReln(reln_type)

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			if (bool_nonconflict_queue == 0):
				fp.write('\n ===>> QUEUE 1 RELN > 1 TREE -- ')
			elif (bool_nonconflict_queue == 1):
				fp.write('\n ===>> QUEUE 2 RELN > 1 TREE -- ')
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
				if (bool_nonconflict_queue == 0):
					fp.write('\n ==>>>>>>>>> NEW CONN --- QUEUE 1 RELN > 1 TREE -- ')
				elif (bool_nonconflict_queue == 1):
					fp.write('\n ==>>>>>>>>> NEW CONN --- QUEUE 2 RELN > 1 TREE -- ')
				else:
					fp.write('\n ==>>>>>>>>> NEW CONN --- CONFLICTING QUEUE -- ')
				fp.write(' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
				#fp.write('  Allowed reln: ' + str(TaxaPair_Reln_Dict[ckey]._GetAllowedRelnList()) + '   Freq: ')
				#for r in TaxaPair_Reln_Dict[ckey]._GetAllowedRelnList():
					#fp.write('  ' + str(TaxaPair_Reln_Dict[ckey]._GetEdgeWeight(r)))
				#fp.write('   Pseudo reln status : ' + str(TaxaPair_Reln_Dict[ckey]._GetFreqPseudoR1(0)) + '/' + \
					#str(TaxaPair_Reln_Dict[ckey]._GetFreqPseudoR1(1)) + '/' + \
						#str(TaxaPair_Reln_Dict[ckey]._GetFreqPseudoR1(2)))
				fp.close()
			
			# also update the reachability graph information
			Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, reln_type)
			
	return Reachability_Graph_Mat
