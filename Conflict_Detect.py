#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

##-------------------------------------------------------
#"""
#this function checks whether there exists a taxa D such that:
#1) src_reln(A,D) and dest_reln(B,D) are the only allowed relations between them, and 
#2) both have not been formed, or one have been formed and other have not been.
#3) these relations conflict with the proposed relation between A and B
#"""
#def PossibleConflict(ReachMat, A, B, A_reln, B_reln, ReachMat_A_idx, ReachMat_B_idx):
	#D_list = list((set(Taxa_Info_Dict[A]._GetAllowedRelnSet(A_reln))).intersection(set(Taxa_Info_Dict[B]._GetAllowedRelnSet(B_reln))))
	#if (len(D_list) > 0):
		#for D in D_list:
			#if (Taxa_Info_Dict[A]._CheckSingleAllowedReln(D, A_reln)) and (Taxa_Info_Dict[B]._CheckSingleAllowedReln(D, B_reln)):
				#D_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(Taxa_Info_Dict[D]._Get_Taxa_Part_Clust_Idx())
				#if (D_reach_mat_idx != ReachMat_A_idx) and (D_reach_mat_idx != ReachMat_B_idx) and \
					#(((ReachMat[D_reach_mat_idx][ReachMat_A_idx] == 0) and (ReachMat[ReachMat_A_idx][D_reach_mat_idx] == 0)) or \
							#((ReachMat[D_reach_mat_idx][ReachMat_B_idx] == 0) and (ReachMat[ReachMat_B_idx][D_reach_mat_idx] == 0))):
					#return 1
	
	#return 0

##-------------------------------------------------------
#"""
#this function checks whether establishment of the "target_reln_type"
#from the "src_taxa_idx" to the "dest_taxa_idx" may lead to future conflict option
#in such a case, a conflict case is returned
#"""
#def Possible_Future_Conflict(src_taxa_idx, dest_taxa_idx, Reachability_Graph_Mat, target_reln_type, Output_Text_File):
	#if (target_reln_type == RELATION_R2):
		#src_taxa_clust_idx = Taxa_Info_Dict[dest_taxa_idx]._Get_Taxa_Part_Clust_Idx()
		#dest_taxa_clust_idx = Taxa_Info_Dict[src_taxa_idx]._Get_Taxa_Part_Clust_Idx()
	#else:
		#src_taxa_clust_idx = Taxa_Info_Dict[src_taxa_idx]._Get_Taxa_Part_Clust_Idx()
		#dest_taxa_clust_idx = Taxa_Info_Dict[dest_taxa_idx]._Get_Taxa_Part_Clust_Idx()
	
	#"""
	#we find the indices of the reachability matrix corresponding to individual cluster indices 
	#"""
	#reach_mat_src_taxa_clust_idx = CURRENT_CLUST_IDX_LIST.index(src_taxa_clust_idx)
	#reach_mat_dest_taxa_clust_idx = CURRENT_CLUST_IDX_LIST.index(dest_taxa_clust_idx)
	
	#if (target_reln_type == RELATION_R2) or (target_reln_type == RELATION_R1):
		#"""
		#if A->B is to be established
		#and suppose there exists (only single allowed relations) of the following combination: 
		#(without being connected)
		#1) D->A and B->D / B><D
		#2) D><A and D->B / B->D
		#3) D><B and D->A / A->D
		#4) B->D and D->A / A><D
		#then return a future conflict case
		#"""
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R2, RELATION_R1, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 1')
				#fp.close()
			#return 1
		
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R2, RELATION_R4, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 2')
				#fp.close()
			#return 1

		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R4, RELATION_R2, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 3')
				#fp.close()
			#return 1

		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R4, RELATION_R1, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 4')
				#fp.close()
			#return 1

		#if (PossibleConflict(Reachability_Graph_Mat, dest_taxa_idx, src_taxa_idx, RELATION_R4, RELATION_R2, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 5')
				#fp.close()
			#return 1

		#if (PossibleConflict(Reachability_Graph_Mat, dest_taxa_idx, src_taxa_idx, RELATION_R4, RELATION_R1, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 6')
				#fp.close()
			#return 1

		#if (PossibleConflict(Reachability_Graph_Mat, dest_taxa_idx, src_taxa_idx, RELATION_R1, RELATION_R2, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 7')
				#fp.close()
			#return 1

		#if (PossibleConflict(Reachability_Graph_Mat, dest_taxa_idx, src_taxa_idx, RELATION_R1, RELATION_R4, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 8')
				#fp.close()
			#return 1

	#if (target_reln_type == RELATION_R3):
		#"""
		#if A==B is to be established
		#and suppose there exists (only single allowed relations) of the following combination: 
		#(without being connected)
		#1) D->A and B->D / B><D
		#2) A->D and D->B / B><D
		#3) D><A and D->B / B->D
		#4) D><B and D->A / A->D
		#5) B->D and D->A / A><D
		#6) D->B and A->D / A><D
		#then return a future conflict case
		#"""
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R2, RELATION_R1, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 1')
				#fp.close()
			#return 1
		
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R2, RELATION_R4, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 2')
				#fp.close()
			#return 1
		
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R1, RELATION_R2, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 3')
				#fp.close()
			#return 1
		
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R1, RELATION_R4, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 4')
				#fp.close()
			#return 1
		
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R4, RELATION_R2, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 5')
				#fp.close()
			#return 1
		
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R4, RELATION_R1, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 6')
				#fp.close()
			#return 1
		
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R2, RELATION_R4, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 7')
				#fp.close()
			#return 1
		
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R1, RELATION_R4, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 8')
				#fp.close()
			#return 1
		
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R2, RELATION_R1, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 9')
				#fp.close()
			#return 1
		
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R4, RELATION_R1, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 10')
				#fp.close()
			#return 1
		
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R1, RELATION_R2, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 11')
				#fp.close()
			#return 1
		
		#if (PossibleConflict(Reachability_Graph_Mat, src_taxa_idx, dest_taxa_idx, RELATION_R4, RELATION_R2, \
			#reach_mat_src_taxa_clust_idx, reach_mat_dest_taxa_clust_idx) == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n target_reln_type : ' + str(target_reln_type) + '  possible future conflict - case 12')
				#fp.close()
			#return 1
		
	##if (target_reln_type == RELATION_R4):

	#return 0

#-------------------------------------------------------
""" 
this function checks whether the current relation (given as an input) is feasible to the input supertree configuration
if the current relation is already existing in the supertree, then this function returns 2
if the current relation produces a cycle to the existing configuration of the final supertree, then it returns 1
otherwise the function returns 0 
"""
def Possible_Conflict_Curr_Reln(src_taxa_idx, dest_taxa_idx, Reachability_Graph_Mat, target_reln_type, Output_Text_File):
  
	if (target_reln_type == RELATION_R2):
		src_taxa_clust_idx = Taxa_Info_Dict[dest_taxa_idx]._Get_Taxa_Part_Clust_Idx()
		dest_taxa_clust_idx = Taxa_Info_Dict[src_taxa_idx]._Get_Taxa_Part_Clust_Idx()
	else:
		src_taxa_clust_idx = Taxa_Info_Dict[src_taxa_idx]._Get_Taxa_Part_Clust_Idx()
		dest_taxa_clust_idx = Taxa_Info_Dict[dest_taxa_idx]._Get_Taxa_Part_Clust_Idx()
	#--------------------------------------------------------------------
	"""
	check whether the couplet is already resolved
	"""
	"""
	if the taxa pair is already part of the same cluster
	that is, they are already related via equivalence partition, then return
	"""
	if (src_taxa_clust_idx == dest_taxa_clust_idx):
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n already in same cluster index')
			fp.close()
		return 1
	
	"""
	we find the indices of the reachability matrix corresponding to individual cluster indices 
	"""
	reach_mat_src_taxa_clust_idx = CURRENT_CLUST_IDX_LIST.index(src_taxa_clust_idx)
	reach_mat_dest_taxa_clust_idx = CURRENT_CLUST_IDX_LIST.index(dest_taxa_clust_idx)

	"""
	case 1 - if the clusters are already related (depicted in the reachability matrix) then return 1
	"""
	if (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][reach_mat_src_taxa_clust_idx] > 0) or \
			(Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][reach_mat_dest_taxa_clust_idx] > 0):
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			if (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][reach_mat_dest_taxa_clust_idx] == 1):
				fp.write('\n target_reln_type: ' + str(target_reln_type) + ' already related via relation r1')
			elif (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][reach_mat_dest_taxa_clust_idx] == 2):
				fp.write('\n target_reln_type: ' + str(target_reln_type) + ' already related via relation r4')
			elif (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][reach_mat_src_taxa_clust_idx] == 1):
				fp.write('\n target_reln_type: ' + str(target_reln_type) + ' already related via relation r2')
			fp.close()
		return 1
	
	#--------------------------------------------------------------------
	# check for the transitive conflict
	if (target_reln_type == RELATION_R2) or (target_reln_type == RELATION_R1):
		
		"""
		if A->B is to be established
		and there exists D->A 
		then if B->D or B><D then return a conflict
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R2):
			if (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1):	#> 0):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n possible conflict - case 1')
					fp.close()
				return 1
		
		"""
		if A->B is to be established
		and there exists B->E 
		then if E->A or E><A then return a conflict  
		"""
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
			if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] == 1):	# > 0):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n possible conflict - case 2')
					fp.close()
				return 1
		
		"""
		if A->B is to be established
		and there exists D->A and B->E 
		then if E->D or E=D or E><D then return a conflict  
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R2):
			for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
				if (x == y) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] == 1):	#> 0):
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n possible conflict - case 3')
						fp.close()
					return 1
				
		"""
		if A->B is to be established
		and there exists D><A
		then it would be D><B afterwards 
		so if there is D->B, B->D or D=B then return a conflict  
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R4):
			#if (dest_taxa_clust_idx == x) \
				#or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] == 1) \
				#or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1):
			if (dest_taxa_clust_idx == x):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n possible conflict - case 4')
					fp.close()
				return 1
		
		"""
		if A->B is to be established
		and there exists D><A
		then for all B->E
		it would be D><E afterwards 
		so if D->E or E->D or D=E then return a conflict  
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R4):
			for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
				#if (x == y) \
					#or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
					#or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 1):
				if (x == y):
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n possible conflict - case 5')
						fp.close()
					return 1
		
	elif (target_reln_type == RELATION_R3):
		
		"""
		if A=B is to be established
		and there exists B->D 
		then if D->A or D><A or D=A then return a conflict    
		"""
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
			if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] == 1):	#> 0):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n possible conflict - case 6')
					fp.close()
				return 1
		
		"""
		if A=B is to be established
		and there exists D->B 
		then if A->D or A><D or A=D then return a conflict    	
		"""
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R2):
			if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1):	#> 0):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n possible conflict - case 7')
					fp.close()
				return 1
		
		"""
		if A=B is to be established
		and there exists B><D 
		then if A->D or D->A or A=D then return a conflict    	
		"""
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R4):
			#if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
					#or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] == 1):
			if (x == src_taxa_clust_idx):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n possible conflict - case 8')
					fp.close()
				return 1

		"""
		if A=B is to be established
		and there exists A->D 
		then if D->B or D><B or D=B then return a conflict    	
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
			if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] == 1):	#> 0):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n possible conflict - case 9')
					fp.close()
				return 1
		
		"""
		if A=B is to be established
		and there exists D->A 
		then if B->D or B><D or B=D then return a conflict    	
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R2):
			if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1):	#> 0):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n possible conflict - case 10')
					fp.close()
				return 1
		
		"""
		if A=B is to be established
		and there exists A><D 
		then if D->B or B->D or B=D then return a conflict    	
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R4):
			#if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
					#or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] == 1):
			if (x == dest_taxa_clust_idx):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n possible conflict - case 11')
					fp.close()
				return 1

	else:	#if (target_reln_type == RELATION_R4):
		"""
		construct the out neighborhood of src_cluster
		it will contain the cluster itself and all the other clusters connected via out edges from this cluster
		"""
		src_clust_out_neighb = []
		src_clust_out_neighb.append(src_taxa_clust_idx)
		src_clust_out_neighb.extend(Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R1))
		
		"""
		construct the out neighborhood of dest cluster
		it will contain the cluster itself and all the other clusters connected via out edges from this cluster    
		"""
		dest_clust_out_neighb = []
		dest_clust_out_neighb.append(dest_taxa_clust_idx)
		dest_clust_out_neighb.extend(Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1))
		
		"""
		if mutually any pair of taxa belonging to respective out edge clusters
		are themselves related via any relationships other than NO EDGE then this connection is not possible
		"""
		for x in src_clust_out_neighb:
			for y in dest_clust_out_neighb:
				#if (x == y) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
					#or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 1):
				if (x == y):
					return 1
		
	return 0
