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

##----------------------------------------------------------------
#"""
#this function checks whether any alternative relation (other than R3) can be established 
#between the clusters containing taxa1 and taxa2
#"""
#def Find_Alternate_Reln_ClustPair(taxa1, taxa2, Reachability_Graph_Mat, outfile):
	#clust1 = Taxa_Info_Dict[taxa1]._Get_Taxa_Part_Clust_Idx()
	#clust2 = Taxa_Info_Dict[taxa2]._Get_Taxa_Part_Clust_Idx()
	#taxa_list1 = Cluster_Info_Dict[clust1]._GetSpeciesList()
	#taxa_list2 = Cluster_Info_Dict[clust2]._GetSpeciesList()
	#R1_score = R2_score = R4_score = 0
	#for x1 in taxa_list1:
		#x1_idx = COMPLETE_INPUT_TAXA_LIST.index(x1)
		#for x2 in taxa_list2:  
			#x2_idx = COMPLETE_INPUT_TAXA_LIST.index(x2)
			#if (x1_idx < x2_idx):
				#target_key = (x1_idx, x2_idx)
				#compl_reln = False
			#else:
				#target_key = (x2_idx, x1_idx)
				#compl_reln = True
			#if target_key in TaxaPair_Reln_Dict:
				#if (compl_reln == False):
					#R1_score = R1_score + TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R1)
					#R2_score = R2_score + TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R2)
					#R4_score = R4_score + TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R4)
				#else:
					#R1_score = R1_score + TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R2)
					#R2_score = R2_score + TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R1)
					#R4_score = R4_score + TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R4)
	
	#"""
	#assemble the scores in the score list and sort the score list in descending order
	#"""
	#score_list = [[RELATION_R1, R1_score], [RELATION_R2, R2_score], [RELATION_R4, R4_score]]
	#score_list.sort(key=lambda x: x[1], reverse=True)
	#if (DEBUG_LEVEL >= 2):
		#fp = open(outfile, 'a')
		#fp.write('\n --- Selecting alternative relation between the cluster pair ')
		#fp.write('\n after sorting the scoring list : ')
		#for i in range(len(score_list)):
			#fp.write('\n relation: ' + str(score_list[i][0]) + ' score: ' + str(score_list[i][1]))
		#fp.close()
	
	#for i in range(len(score_list)):
		#target_reln = score_list[i][0]
		#if (Possible_Conflict_Curr_Reln(taxa1, taxa2, Reachability_Graph_Mat, target_reln, outfile) == 0):
			#return target_reln
		#else:
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n --- Conflict found while checking the alternative relation: ' + str(target_reln))
				#fp.close()
	
	#return UNDEFINED_RELATION

#-------------------------------------------------------
"""
this function checks whether a given pair of taxa clusters can be merged (establish relation R3 between them)
"""
def Check_Merge_Clust_Possible(clust1, clust2):
	no_of_couplets = 0
	taxa_list1 = Cluster_Info_Dict[clust1]._GetSpeciesList()
	taxa_list2 = Cluster_Info_Dict[clust2]._GetSpeciesList()
	for x1 in taxa_list1:
		x1_idx = COMPLETE_INPUT_TAXA_LIST.index(x1)
		for x2 in taxa_list2:  
			x2_idx = COMPLETE_INPUT_TAXA_LIST.index(x2)
			if (x1_idx < x2_idx):
				target_key = (x1_idx, x2_idx)
			else:
				target_key = (x2_idx, x1_idx)
			if target_key in TaxaPair_Reln_Dict:
				no_of_couplets = no_of_couplets + 1
				if (RELATION_R3 not in TaxaPair_Reln_Dict[target_key]._GetAllowedRelnList()):
					return False
				
	if (no_of_couplets == 0):
		return False
	
	return True


#-------------------------------------------------------
""" 
this function processes support score queue containing only R3 relation information
the objective is to create cluster of taxa
"""
def Proc_Queue_RelnR3(Reachability_Graph_Mat, Output_Text_File, inp_no):
	if (inp_no == 1):
		Inp_Queue = Queue_Score_R3_SingleReln
	else:
		Inp_Queue = Queue_Score_R3_MajCons
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
		reln_freq = outlist[3]
		conn_score = outlist[4]
	
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			if (inp_no == 1):
				fp.write('\n ===>> SUPPORT SCORE QUEUE -- NonConflict R3 RELATION -- ')      
			else:
				fp.write('\n ===>> SUPPORT SCORE QUEUE -- Majority Consensus R3 RELATION -- ')      
			fp.write(' current extracted max element: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
					' relation type: ' + str(reln_type) + ' frequency: ' + str(reln_freq) + ' conn score: ' + str(conn_score))
			fp.close()

		"""
		there is no concept of conflict in this case
		we just check whether R3 relation is predominant among all taxa pairs within this pair of cluster
		"""
		clust1 = Taxa_Info_Dict[src_taxa_idx]._Get_Taxa_Part_Clust_Idx()
		clust2 = Taxa_Info_Dict[dest_taxa_idx]._Get_Taxa_Part_Clust_Idx()
		
		if (clust1 == clust2):
			if (DEBUG_LEVEL > 0):
				fp = open(Output_Text_File, 'a')    
				fp.write('\n ---- Already in the same cluster index --- ')
				fp.close()
		else:
			r3_possible = Check_Merge_Clust_Possible(clust1, clust2)
			if (r3_possible == True):
				"""
				the cluster pair can be merged (relation R3)
				"""
				if (DEBUG_LEVEL > 0):
					fp = open(Output_Text_File, 'a')    
					fp.write('\n ==>>>>>>>>> NEW CONN --- CONFLICTING QUEUE -- relation type: ' + str(reln_type) \
						+ ' frequency: ' + str(reln_freq) + ' conn score: ' + str(conn_score))
					fp.close()
				"""
				also update the reachability graph information
				"""
				Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, clust1, clust2, reln_type, Output_Text_File)
	
	return Reachability_Graph_Mat

#-------------------------------------------------------
""" 
this function processes the support score queue designed to contain the frequencies and 
support scores for individual relations between a pair of cluster
"""
def Proc_Queue_Clust(Reachability_Graph_Mat, Output_Text_File):
	Inp_Queue = Queue_Score_Cluster_Pair

	while (0 < len(Inp_Queue)):
		""" 
		extract the 1st element of "Inp_Queue" 
		since it is sorted to have max cost at the beginning 
		"""
		outlist = Heap_Extract_Max(Inp_Queue)
		clust1 = outlist[0]
		clust2 = outlist[1]
		reln_type = outlist[2]
		reln_freq = outlist[3]
		conn_score = outlist[4]

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n ===>> CLUSTER BASED SUPPORT SCORE QUEUE -- ')      
			fp.write(' current extracted max element (cluster pair): ' + str(clust1) + ' and ' + str(clust2) + \
					' relation type: ' + str(reln_type) + '  relation freq: ' + str(reln_freq) + '  conn score: ' + str(conn_score))
			fp.close()

		"""
		if the current extracted relation does not induce a conflict to the existing configuration of the DAG, 
		include the connection in it 
		"""
		conflict_detection = Possible_Conflict_Curr_Reln(clust1, clust2, Reachability_Graph_Mat, reln_type, Output_Text_File)
		
		if (conflict_detection == 0):
			""" 
			current element does not create a cycle / conflict
			not that it is already present in the supertree 
			valid connection is found - append this connection to the final formed tree 
			"""
			if (DEBUG_LEVEL > 0):
				fp = open(Output_Text_File, 'a')    
				fp.write('\n ==>>>>>>>>> NEW CONN --- CONFLICTING QUEUE -- relation type: ' + str(reln_type) \
					+ ' relation freq: ' + str(reln_freq) + ' conn score: ' + str(conn_score))
				fp.close()

			"""
			also update the reachability graph information
			"""
			Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, clust1, clust2, reln_type, Output_Text_File)
			
	return Reachability_Graph_Mat
