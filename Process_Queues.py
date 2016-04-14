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

##---------------------------------------------------------------------
#"""
#this function processes the set of taxa belonging to the allowed taxa set of the "Inp_Taxa" 
#with respect to the "inp_reln"
#"""
#def Process_allowed_reln_set(Inp_Taxa, ReachMat, Outfile, inp_reln):
	#"""
	#first reset the support score queue
	#"""
	#Cost_List_Taxa_Pair_Multi_Reln[:] = []
	#"""
	#extract the allowed set of taxa corresponding to the "inp_reln"
	#"""
	#inp_reln_allowed_set = Taxa_Info_Dict[Inp_Taxa]._GetAllowedRelnSet(inp_reln)
	#if (len(inp_reln_allowed_set) > 0):
		#"""
		#form the support score queue with the allowed taxa set of R4 relations
		#"""
		#for t in inp_reln_allowed_set:
			#if (Inp_Taxa < t):
				#target_key = (Inp_Taxa, t)
				#target_reln = inp_reln
			#else:
				#target_key = (t, Inp_Taxa)
				#target_reln = Complementary_Reln(inp_reln)
				
			#if target_key in TaxaPair_Reln_Dict:
				#sublist = [target_key[0], target_key[1], target_reln, TaxaPair_Reln_Dict[target_key]._GetEdgeCost_ConnReln(target_reln)]
				#Cost_List_Taxa_Pair_Multi_Reln.append(sublist)
				
		#if (len(Cost_List_Taxa_Pair_Multi_Reln) > 0):
			#"""
			#sort the support score queue
			#"""
			#Sort_Priority_Queue(Cost_List_Taxa_Pair_Multi_Reln)
			
			#"""
			#now within a loop, process this support score queue, to establish the couplet based connections
			#"""
			#while (0 < len(Cost_List_Taxa_Pair_Multi_Reln)):
				#""" 
				#extract the 1st element of "Inp_Queue" 
				#since it is sorted to have max cost at the beginning 
				#"""
				#outlist = Heap_Extract_Max(Cost_List_Taxa_Pair_Multi_Reln)
				#src_taxa_idx = outlist[0]
				#src_taxa_label = COMPLETE_INPUT_TAXA_LIST[src_taxa_idx]
				#dest_taxa_idx = outlist[1]
				#dest_taxa_label = COMPLETE_INPUT_TAXA_LIST[dest_taxa_idx]
				#reln_type = outlist[2]
				#conn_score = TaxaPair_Reln_Dict[(src_taxa_idx, dest_taxa_idx)]._GetEdgeCost_ConnReln(reln_type)

				#if (DEBUG_LEVEL >= 2):
					#fp = open(Outfile, 'a')
					#fp.write('\n ===>> CONFLICTING QUEUE -- ')      
					#fp.write(' current extracted max element: ' + str(src_taxa_idx) + '(' + str(src_taxa_label) + \
						#') and ' + str(dest_taxa_idx) + '(' + str(dest_taxa_label) + ')' + \
							#' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
					#fp.close()

				#""" 
				#if the current score metric based relation does not induce a cycle to 
				#the existing configuration of the final supertree
				#then include the connection in it 
				#"""
				#conflict_detection = Possible_Conflict_Curr_Reln(src_taxa_idx, dest_taxa_idx, ReachMat, reln_type, Outfile)
				
				#if (conflict_detection == 0):
					#if (Possible_Future_Conflict(src_taxa_idx, dest_taxa_idx, ReachMat, reln_type, Outfile) == 0):
						#""" 
						#current element does not create a cycle / conflict
						#not that it is already present in the supertree 
						#valid connection is found - append this connection to the final formed tree 
						#"""
						#if (DEBUG_LEVEL > 0):
							#fp = open(Outfile, 'a')    
							#fp.write('\n ==>>>>>>>>> NEW CONN --- CONFLICTING QUEUE -- ')
							#fp.write(' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
							#fp.close()
						
						## also update the reachability graph information
						#ReachMat = AdjustReachGraph(ReachMat, src_taxa_idx, dest_taxa_idx, reln_type)
				
	#return ReachMat

##---------------------------------------------------------------------
#"""
#this function processes one taxa at a time, mainly its connections with respect to the input gene trees, 
#to establish the couplet based connections in the output supertree
#"""
#def Process_Single_Taxa(Inp_Taxa, ReachMat, Outfile):
	#ReachMat = Process_allowed_reln_set(Inp_Taxa, ReachMat, Outfile, RELATION_R4)
	#ReachMat = Process_allowed_reln_set(Inp_Taxa, ReachMat, Outfile, RELATION_R1)
	##ReachMat = Process_allowed_reln_set(Inp_Taxa, ReachMat, Outfile, RELATION_R2)
	#ReachMat = Process_allowed_reln_set(Inp_Taxa, ReachMat, Outfile, RELATION_R3)
	
	#return ReachMat

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
			if 1:	#(Possible_Future_Conflict(src_taxa_idx, dest_taxa_idx, Reachability_Graph_Mat, reln_type, Output_Text_File) == 0):
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
