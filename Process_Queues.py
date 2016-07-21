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
this function checks whether a given pair of taxa clusters can be merged (establish relation R3 between them)
"""
def Check_Merge_Couplet_Possible(clust1, clust2, outfile):
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
					if (DEBUG_LEVEL >= 2):
						fp = open(outfile, 'a')
						fp.write('\n Here the clusters have taxon ' + str(x1) + ' and ' + str(x2) + \
							'  which do not have R3 as their allowed relation -- no merge cluster possible ')
						fp.close()
					return False
				
	if (no_of_couplets == 0):
		return False
	
	return True

#-------------------------------------------------------
""" 
this function processes support score queue containing couplet based relations
when the input relation is R3, the objective is to create cluster of taxa
otherwise, for other cases, corresponding couplets are merged
"""
def Proc_Queue_Couplet_Reln_R3(Output_Text_File, inp_no):
	"""
	select input support score queue according to the number given in input parameter
	"""
	if (inp_no == 1):
		Inp_Queue = Queue_Score_R3_SingleReln
	else:	#if (inp_no == 2):
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
		
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			if (inp_no == 1):
				fp.write('\n ===>> SUPPORT SCORE QUEUE -- NonConflict R3 RELATION -- ')      
			else:	#if (inp_no == 2):
				fp.write('\n ===>> SUPPORT SCORE QUEUE -- Majority Consensus R3 RELATION -- ')      
			fp.write(' current extracted max element: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
					' relation type: ' + str(reln_type) + ' frequency: ' + str(reln_freq))	# + ' conn score: ' + str(conn_score))
			fp.close()

		"""
		there is no concept of conflict in this case
		we just check whether R3 relation is predominant among all taxa pairs within this pair of cluster
		"""
		clust1 = Taxa_Info_Dict[src_taxa_idx]._Get_Taxa_Part_Clust_Idx()
		clust2 = Taxa_Info_Dict[dest_taxa_idx]._Get_Taxa_Part_Clust_Idx()
			
		"""
		for R3 input relations - case of cluster merging
		"""
		if (clust1 == clust2):
			if (DEBUG_LEVEL > 0):
				fp = open(Output_Text_File, 'a')    
				fp.write('\n ---- Already in the same cluster index --- ')
				fp.close()
		else:
			r3_possible = Check_Merge_Couplet_Possible(clust1, clust2, Output_Text_File)
			if (r3_possible == True):
				"""
				the cluster pair can be merged (relation R3)
				"""
				if (DEBUG_LEVEL > 0):
					fp = open(Output_Text_File, 'a')    
					fp.write('\n ==>>>>>>>>> NEW CONN --- CONFLICTING QUEUE -- relation type: ' + str(reln_type) \
						+ ' frequency: ' + str(reln_freq))	# + ' conn score: ' + str(conn_score))
					fp.close()
				"""
				update the cluster connectivity information
				"""
				if (clust1 < clust2):
					Copy_Cluster_Content(clust1, clust2, Output_Text_File)
				else:
					Copy_Cluster_Content(clust2, clust1, Output_Text_File)

	return 

#-------------------------------------------------------
""" 
this function processes the support score queue designed to contain the frequencies and 
support scores for individual relations between a pair of cluster
"""
def Proc_Queue_Clust(Reachability_Graph_Mat, Output_Text_File, inp_no):
	if (inp_no == 1):
		Inp_Queue = Queue_Score_Cluster_Pair_NonConflict
	else:
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
			if (inp_no == 1):
				fp.write('\n\n ===>> CLUSTER BASED (*** NONCONFLICTING ***) SUPPORT SCORE QUEUE -- ')      
			else:
				fp.write('\n\n ===>> CLUSTER BASED (*** CONFLICTING ***) SUPPORT SCORE QUEUE -- ')      
			fp.write(' current extracted max element (cluster pair): ' + str(clust1) + ' and ' + str(clust2) + \
					' relation type: ' + str(reln_type) + '  relation freq: ' + str(reln_freq) + '  conn score: ' + str(conn_score))
			fp.close()

		"""
		Note: we have re-written the conflict detection routine
		if the current extracted relation does not induce a conflict to the existing configuration of the DAG, 
		include the connection in it 
		"""
		conflict_detection = Possible_Conflict_Curr_Reln(clust1, clust2, Reachability_Graph_Mat, reln_type, Output_Text_File)
		
		if (conflict_detection == 0):
			""" 
			current element does not create a cycle / conflict
			not that it is already present in the supertree 
			"""
			if (DEBUG_LEVEL > 0):
				fp = open(Output_Text_File, 'a')    
				if (inp_no == 1):
					fp.write('\n ==>>>>>>>>> NEW CONN --- NONCONFLICTING QUEUE')
				else:
					fp.write('\n ==>>>>>>>>> NEW CONN --- CONFLICTING QUEUE')
				fp.write('-- relation type: ' + str(reln_type) \
					+ '  from clust: ' + str(clust1) + ' to clust: ' + str(clust2) \
						+ ' relation freq: ' + str(reln_freq) + ' conn score: ' + str(conn_score))
				fp.close()
			"""
			also update the reachability graph information
			"""
			Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, clust1, clust2, reln_type, Output_Text_File)
		
		else:

			if (conflict_detection == 1):
				clust_pair_key = (clust1, clust2)
				Cluster_Pair_Info_Dict[clust_pair_key]._RemovePossibleReln(reln_type)
				if (DEBUG_LEVEL > 0):
					fp = open(Output_Text_File, 'a')    
					fp.write('\n Conflict detection output: ' + str(conflict_detection) + '  removed possible relation: ' + str(reln_type))
					fp.close()
				
				"""
				this condition is enforced only if the relation is not at all applicable to the cluster pair
				that's why we have checked whether the value of "conflict_detection" is 1, not 2
				
				then check whether there exists only one possible relation among this cluster pair and 
				that relation is R4
				if the R4 relation is non-conflicting then apply the relation between this pair of cluster
				"""
				allowed_reln_list = Cluster_Pair_Info_Dict[clust_pair_key]._GetPossibleRelnList()
				
				if (reln_type != RELATION_R4):
					if ((len(allowed_reln_list) == 1) and (RELATION_R4 in allowed_reln_list)) or (len(allowed_reln_list) == 0):
						if (Possible_Conflict_Curr_Reln(clust1, clust2, Reachability_Graph_Mat, RELATION_R4, Output_Text_File) == 0):
							if (DEBUG_LEVEL > 0):
								fp = open(Output_Text_File, 'a') 
								fp.write('\n Only one allowed relation is remaining  and that is the relation R4 - it is non-conflicting as well')
								if (inp_no == 1):
									fp.write('\n ==>>>>>>>>> NEW CONN --- NONCONFLICTING QUEUE')
								else:
									fp.write('\n ==>>>>>>>>> NEW CONN --- CONFLICTING QUEUE')
								fp.write('-- relation type: ' + str(RELATION_R4) \
									+ '  from clust: ' + str(clust1) + ' to clust: ' + str(clust2) \
										+ ' relation freq: ' + str(Cluster_Pair_Info_Dict[clust_pair_key]._GetFreq(RELATION_R4)) \
											+ ' conn score: ' + str(Cluster_Pair_Info_Dict[clust_pair_key]._GetSupportScore(RELATION_R4)))
								fp.close()
							"""
							also update the reachability graph information
							"""
							Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, clust1, clust2, RELATION_R4, Output_Text_File)
						
	return Reachability_Graph_Mat
