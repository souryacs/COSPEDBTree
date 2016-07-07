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
#This function is a module of the function Check_Merge_Clust_Possible

#@input: x1_idx (index of one taxa)
				#x2_idx (index of other taxa)

#the following condition is checked:

#suppose there exists other taxon z (not equal to x1 or x2)
#suppose, one of the following cases happen:

#1) a) R1(x1,z) is a Consensus relation, b) R2(x2,z) is consensus c) frequency of R1(x2,z) = 0
#2) a) R1(x1,z) is a Consensus relation, b) R2(x2,z) is consensus c) frequency of R2(x1,z) = 0
#3) a) R2(x1,z) is a Consensus relation, b) R1(x2,z) is consensus c) frequency of R2(x2,z) = 0
#4) a) R2(x1,z) is a Consensus relation, b) R1(x2,z) is consensus c) frequency of R1(x1,z) = 0

#In such a case, merging of the clusters is not performed

#@ return type:
#True: if merging of clusters is not allowed
#False: if merging of clusters is allowed

#"""
#def ConflictingR1R2Reln(x1_idx, x2_idx, outfile):
	#for z_idx in range(len(COMPLETE_INPUT_TAXA_LIST)):
		#if (z_idx == x1_idx) or (z_idx == x2_idx):
			#continue
		
		#"""
		#determine the frequency of R1 relation from x1 to z
		#also determine if the relation is consensus or not
		#"""
		#freq_r1_x1_idx_z_idx = GetRelnFreq(x1_idx, z_idx, RELATION_R1)
		#consenus_r1_x1_idx_z_idx = CheckConsensusReln(x1_idx, z_idx, RELATION_R1)
		
		#"""
		#determine the frequency of R2 relation from x1 to z
		#also determine if the relation is consensus or not
		#"""
		#freq_r2_x1_idx_z_idx = GetRelnFreq(x1_idx, z_idx, RELATION_R2)
		#consenus_r2_x1_idx_z_idx = CheckConsensusReln(x1_idx, z_idx, RELATION_R2)
		
		#"""
		#determine the frequency of R1 relation from x2 to z
		#also determine if the relation is consensus or not
		#"""
		#freq_r1_x2_idx_z_idx = GetRelnFreq(x2_idx, z_idx, RELATION_R1)
		#consenus_r1_x2_idx_z_idx = CheckConsensusReln(x2_idx, z_idx, RELATION_R1)
		
		#"""
		#determine the frequency of R2 relation from x2 to z
		#also determine if the relation is consensus or not
		#"""
		#freq_r2_x2_idx_z_idx = GetRelnFreq(x2_idx, z_idx, RELATION_R2)
		#consenus_r2_x2_idx_z_idx = CheckConsensusReln(x2_idx, z_idx, RELATION_R2)
		
		#"""
		#cases a to d are mentioned in these conditions
		#"""
		#if (consenus_r1_x1_idx_z_idx == 1) and (consenus_r2_x2_idx_z_idx == 1) and (freq_r1_x2_idx_z_idx == 0):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n x1: ' + str(COMPLETE_INPUT_TAXA_LIST[x1_idx]) + ' x2: ' + str(COMPLETE_INPUT_TAXA_LIST[x2_idx]) + \
					#' z: ' + str(COMPLETE_INPUT_TAXA_LIST[z_idx]) + ' case a -- no cluster merging possible')
				#fp.close()
			#return True
		
		#if (consenus_r1_x1_idx_z_idx== 1) and (consenus_r2_x2_idx_z_idx == 1) and (freq_r2_x1_idx_z_idx == 0):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n x1: ' + str(COMPLETE_INPUT_TAXA_LIST[x1_idx]) + ' x2: ' + str(COMPLETE_INPUT_TAXA_LIST[x2_idx]) + \
					#' z: ' + str(COMPLETE_INPUT_TAXA_LIST[z_idx]) + ' case b -- no cluster merging possible')
				#fp.close()
			#return True

		#if (consenus_r2_x1_idx_z_idx == 1) and (consenus_r1_x2_idx_z_idx == 1) and (freq_r2_x2_idx_z_idx == 0):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n x1: ' + str(COMPLETE_INPUT_TAXA_LIST[x1_idx]) + ' x2: ' + str(COMPLETE_INPUT_TAXA_LIST[x2_idx]) + \
					#' z: ' + str(COMPLETE_INPUT_TAXA_LIST[z_idx]) + ' case c -- no cluster merging possible')
				#fp.close()
			#return True

		#if (consenus_r2_x1_idx_z_idx == 1) and (consenus_r1_x2_idx_z_idx == 1) and (freq_r1_x1_idx_z_idx == 0):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n x1: ' + str(COMPLETE_INPUT_TAXA_LIST[x1_idx]) + ' x2: ' + str(COMPLETE_INPUT_TAXA_LIST[x2_idx]) + \
					#' z: ' + str(COMPLETE_INPUT_TAXA_LIST[z_idx]) + ' case d -- no cluster merging possible')
				#fp.close()
			#return True

		##if (freq_r1_x1_idx_z_idx > 0) and (freq_r2_x1_idx_z_idx == 0) and (freq_r1_x2_idx_z_idx == 0) and (freq_r2_x2_idx_z_idx > 0):
			##return True
		##if (freq_r1_x1_idx_z_idx == 0) and (freq_r2_x1_idx_z_idx > 0) and (freq_r1_x2_idx_z_idx > 0) and (freq_r2_x2_idx_z_idx == 0):
			##return True
		
	#return False


##-------------------------------------------------------
#"""
#This function is a module of the function Check_Merge_Clust_Possible

#@input: x1_idx (index of one taxa x1)
				#x2_idx (index of other taxa x2)

#target relation: R3(x1, x2)

#the following condition is checked:

#suppose there exists other taxon z (not equal to x1 or x2)
#suppose, one of the following cases happen:

#1) a) R1(x1,z) is consensus, b) R3(x2,z) is allowed 
#2) a) R1(x2,z) is consensus, b) R3(x1,z) is allowed 
#3) a) R1(x1,z) is single allowed, b) R2(x2, z) is consensus
#4) a) R2(x1,z) is single allowed, b) R1(x2, z) is consensus

#In such a case, merging of the clusters is not performed

#@ return type:
#True: if merging of clusters is not allowed
#False: if merging of clusters is allowed

#"""
#def ConflictingR1R2Reln(x1_idx, x2_idx, outfile):
	#for z_idx in range(len(COMPLETE_INPUT_TAXA_LIST)):
		#if (z_idx == x1_idx) or (z_idx == x2_idx):
			#continue
		
		#if (CheckConsensusReln(x1_idx, z_idx, RELATION_R1) == True) and (CheckAllowedReln(x2_idx, z_idx, RELATION_R3) == True):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n x1: ' + str(COMPLETE_INPUT_TAXA_LIST[x1_idx]) + ' x2: ' + str(COMPLETE_INPUT_TAXA_LIST[x2_idx]) + \
					#' z: ' + str(COMPLETE_INPUT_TAXA_LIST[z_idx]) + ' case a -- no couplet merging possible')
				#fp.close()
			#return True
		
		#if (CheckConsensusReln(x2_idx, z_idx, RELATION_R1) == True) and (CheckAllowedReln(x1_idx, z_idx, RELATION_R3) == True):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n x1: ' + str(COMPLETE_INPUT_TAXA_LIST[x1_idx]) + ' x2: ' + str(COMPLETE_INPUT_TAXA_LIST[x2_idx]) + \
					#' z: ' + str(COMPLETE_INPUT_TAXA_LIST[z_idx]) + ' case b -- no couplet merging possible')
				#fp.close()
			#return True

		#if (CheckSingleAllowedReln(x1_idx, z_idx, RELATION_R1) == True) and (CheckConsensusReln(x2_idx, z_idx, RELATION_R2) == True):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n x1: ' + str(COMPLETE_INPUT_TAXA_LIST[x1_idx]) + ' x2: ' + str(COMPLETE_INPUT_TAXA_LIST[x2_idx]) + \
					#' z: ' + str(COMPLETE_INPUT_TAXA_LIST[z_idx]) + ' case c -- no couplet merging possible')
				#fp.close()
			#return True

		#if (CheckSingleAllowedReln(x1_idx, z_idx, RELATION_R2) == True) and (CheckConsensusReln(x2_idx, z_idx, RELATION_R1) == True):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n x1: ' + str(COMPLETE_INPUT_TAXA_LIST[x1_idx]) + ' x2: ' + str(COMPLETE_INPUT_TAXA_LIST[x2_idx]) + \
					#' z: ' + str(COMPLETE_INPUT_TAXA_LIST[z_idx]) + ' case d -- no couplet merging possible')
				#fp.close()
			#return True
		
	#return False

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
				## add - sourya
				#if (ConflictingR1R2Reln(x1_idx, x2_idx, outfile) == True):
					#return False
				## end add - sourya
				
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
		#conn_score = outlist[4]
		
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
#"""
#this function checks whether the "inp_reln" from "clust1" to "clust2" 
#is allowed and Consensus
#"""
#def Check_allowed_consensus_Reln(clust1, clust2, inp_reln):
	#if (clust1 < clust2):
		#clust_pair_key = (clust1, clust2)
		#if clust_pair_key in Cluster_Pair_Info_Dict:
			#return Cluster_Pair_Info_Dict[clust_pair_key]._CheckAllowedConsensus(inp_reln)
	#else:
		#clust_pair_key = (clust2, clust1)
		#if clust_pair_key in Cluster_Pair_Info_Dict:
			#return Cluster_Pair_Info_Dict[clust_pair_key]._CheckAllowedConsensus(Complementary_Reln(inp_reln))

#-------------------------------------------------------
#"""
#this function  checks whether a R1 relation can be established from the clust1 to the clust2
#"""
#def Check_R1_Reln_Possible(clust1, clust2, ReachMat, outfile, reln_list):
	#clust2_R2_list = Cluster_Info_Dict[clust2]._GetClustRelnList(RELATION_R2)

	#if (DEBUG_LEVEL >= 2):
		#fp = open(outfile, 'a')
		#fp.write('\n Target relation: ' + str(clust1) + '->' + str(clust2))
		#fp.close()

	#if (len(clust2_R2_list) == 0):
		#"""
		#there exists no cluster x as of now  such that x->clust2 is True
		#so, we can establish clust1->clust2
		#"""
		#if (DEBUG_LEVEL >= 2):
			#fp = open(outfile, 'a')
			#fp.write('\n Here the cluster ' + str(clust2) + '  has no R2 relation -- return True')
			#fp.close()
		#return reln_list, True
	#else:
		#"""
		#explore all the clusters x such that x->clust2 is True
		#"""
		##Target relation: 7->6
		
		#flag = True
		#for i in range(len(clust2_R2_list)):
			#x = clust2_R2_list[i]
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n Examining the in edge (R2) cluster: ' + str(x))
				#fp.close()
			#"""
			#given: x->clust2
			#objective: clust1->clust2
			#case 1: clust1->x is already present. So, clust1->clust2 would already be True. As the relation is not present yet, 
			#so clust1->x is not True yet
			#case 2: clust1->x is not present, but is a Consensus and allowed relation. Then, only 
			#clust1->x need to be established, since it will automatically imply the proposed relation.
			#case 3: x->clust1 is True. Then, the proposed relation is allowed.
			#case 4: x->clust1 is a Consensus and allowed relation. Then, x->clust1 and the proposed relation need to be set.
			#"""
			#if (ReachMat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(clust1)] == 1) and \
				#(ReachMat[CURRENT_CLUST_IDX_LIST.index(clust1)][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
				#"""
				#case 3 is satisfied
				#"""
				#if (DEBUG_LEVEL >= 2):
					#fp = open(outfile, 'a')
					#fp.write('\n Here the cluster ' + str(clust2) + '  has R2 relation cluster : ' + str(x) + \
						#' and the relation ' + str(x) + '->' + str(clust1) + '  is already true - so proceed' )
					#fp.close()
					
			#elif (Check_allowed_consensus_Reln(clust1, x, RELATION_R1) == True):
				#"""
				#case 2 is satisfied
				#"""
				#subl = [clust1, x, RELATION_R1]
				#reln_list.append(subl)
				#if (DEBUG_LEVEL >= 2):
					#fp = open(outfile, 'a')
					#fp.write('\n Here the cluster ' + str(clust2) + '  has R2 relation cluster : ' + str(x) + \
						#' and the relation ' + str(clust1) + '->' + str(x) + ' is allowed and Consensus - update reln_list and proceed' )
					#fp.close()

			#elif (Check_allowed_consensus_Reln(x, clust1, RELATION_R1) == True):
				#"""
				#case 4 is satisfied
				#"""
				#subl = [x, clust1, RELATION_R1]
				#reln_list.append(subl)
				#if (DEBUG_LEVEL >= 2):
					#fp = open(outfile, 'a')
					#fp.write('\n Here the cluster ' + str(clust2) + '  has R2 relation cluster : ' + str(x) + \
						#' and the relation ' + str(x) + '->' + str(clust1) + ' is allowed and Consensus - update reln_list and proceed' )
					#fp.close()

			#else:
				#flag = False
				#break
	
	#if (flag == False):
		#if (DEBUG_LEVEL >= 2):
			#fp = open(outfile, 'a')
			#fp.write('\n Target relation: ' + str(clust1) + '->' + str(clust2) + ' could not be established' )
			#fp.close()
	
	#return reln_list, flag

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
		outlist = Heap_Extract_Max(Inp_Queue)	#, inp_no)
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
		# comment - sourya
		conflict_detection = Possible_Conflict_Curr_Reln(clust1, clust2, Reachability_Graph_Mat, reln_type, Output_Text_File)
		# add - sourya
		#conflict_detection = Check_Conflict(clust1, clust2, Reachability_Graph_Mat, reln_type, Output_Text_File)
		
		if (conflict_detection == 0):
			""" 
			current element does not create a cycle / conflict
			not that it is already present in the supertree 
			"""
			#"""
			#we apply one more checking, provided the relation is either R1 or R2
			#"""
			#Conn_Allowed = True
			#Conn_Allowed1 = True
			
			#if (reln_type == RELATION_R1): 
				##Conn_Allowed = CheckR1RelnConflict(Reachability_Graph_Mat, clust1, clust2, Output_Text_File)
				#Conn_Allowed = CheckR1RelnNoPossibleCycle(Reachability_Graph_Mat, clust1, clust2, Output_Text_File)
				#Conn_Allowed1 = CheckR1RelnMPP(Reachability_Graph_Mat, clust1, clust2, Output_Text_File)
			#elif (reln_type == RELATION_R2):
				##Conn_Allowed = CheckR1RelnConflict(Reachability_Graph_Mat, clust2, clust1, Output_Text_File)
				#Conn_Allowed = CheckR1RelnNoPossibleCycle(Reachability_Graph_Mat, clust2, clust1, Output_Text_File)
				#Conn_Allowed1 = CheckR1RelnMPP(Reachability_Graph_Mat, clust2, clust1, Output_Text_File)
				
			#if (Conn_Allowed == True) and (Conn_Allowed1 == True):
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
			
			#else:
				
				## add - sourya - current relation is not allowed
				## remove the relation from the allowed relation list between this cluster pair
				#clust_pair_key = (clust1, clust2)
				#Cluster_Pair_Info_Dict[clust_pair_key]._RemovePossibleReln(reln_type)
				#if (DEBUG_LEVEL > 0):
					#fp = open(Output_Text_File, 'a')    
					#fp.write('\n Current relation could not be established --- removed possible relation: ' + str(reln_type))
					#fp.close()
		
		else:

			#"""
			#current relation either produces a cycle / conflict 
			#or there exists already some relation between this cluster pair
			#first remove the relation from the set of allowed relations between this cluster pair
			#"""
			#clust_pair_key = (clust1, clust2)
			#Cluster_Pair_Info_Dict[clust_pair_key]._RemovePossibleReln(reln_type)

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
