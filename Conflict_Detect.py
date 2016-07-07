#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

##-------------------------------------------------------
#"""
#this function checks about the possible connection clust1->clust2
#and checks whether it will yield to the MPP (2) problem
#that is, x->y and z->y and (x,z) are related via the relation R4
#Returns:
#True if the connection is possible
#False if no such connection is possible
#"""
#def CheckR1RelnMPP(ReachMat, clust1, clust2, outfile):
	#"""
	#check the already existing connections z->clust2, if available
	#"""
	#reachmat_idx_clust1 = CURRENT_CLUST_IDX_LIST.index(clust1)
	#clust2_R2_list = Cluster_Info_Dict[clust2]._GetClustRelnList(RELATION_R2)
	#for z in clust2_R2_list:
		#"""
		#if there exists z >< clust1 aready, 
		#or only R4(z, clust1) is allowed
		#then we discard the connection 
		#"""
		#reachmat_idx_z = CURRENT_CLUST_IDX_LIST.index(z)
		#if (ReachMat[reachmat_idx_z][reachmat_idx_clust1] == 2):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n *** Within function CheckR1RelnMPP --- Possible conflict--- target relation: ' + str(clust1) + '->' + str(clust2) + \
						#' But there is a cluster ' + str(z) + ' such that ' + str(z) + '->' + str(clust2) + \
							#' And ' + str(z) + '><' + str(clust1) + ' is already established')
				#fp.close()
			#return False
		
		#allowed_reln_list = GetAllowedRelnClusterPair(clust1, clust2)
		#if (len(allowed_reln_list) == 1) and (RELATION_R4 in allowed_reln_list):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n *** Within function CheckR1RelnMPP --- Possible conflict--- target relation: ' + str(clust1) + '->' + str(clust2) + \
						#' But there is a cluster ' + str(z) + ' such that ' + str(z) + '->' + str(clust2) + \
							#' And ' + str(z) + '><' + str(clust1) + ' is only allowed')
				#fp.close()
			#return False
	
	#return True

##------------------------------------------------------------
#"""
#this function checks whether a R1 relation from "clust1" to "clust2" is feasible
#Returns:
#True if the connection is possible
#False if no such connection is possible
#"""
#def CheckR1RelnNoPossibleCycle(ReachMat, clust1, clust2, outfile):
	#"""
	#target relation: clust1->clust2
	#we check whether there exists a connection clust2->clust3 already
	#such that clust3->clust1 is allowed (may not be included)
	#"""
	#reachmat_idx_clust1 = CURRENT_CLUST_IDX_LIST.index(clust1)
	#reachmat_idx_clust2 = CURRENT_CLUST_IDX_LIST.index(clust2)
	#for i in range(len(CURRENT_CLUST_IDX_LIST)):
		#if (ReachMat[reachmat_idx_clust2][i] == 1):
			#reachmat_idx_clust3 = i
			#clust3 = CURRENT_CLUST_IDX_LIST[reachmat_idx_clust3]
			#if (CheckAllowedRelnClusterPair(clust3, clust1, RELATION_R1) == True):
				#if (DEBUG_LEVEL >= 2):
					#fp = open(outfile, 'a')
					#fp.write('\n *** Within function CheckR1RelnNoPossibleCycle --- Possible conflict --- target relation: ' + str(clust1) + '->' + str(clust2) + \
						#' But there is a cluster ' + str(clust3) + ' such that ' + str(clust2) + '->' + str(clust3) + \
							#' And ' + str(clust3) + '->' + str(clust1) + ' is allowed')
					#fp.close()
				#return False

	#return True

##--------------------------------------------------
#"""
#this function checks whether a R1 relation from "clust1" to "clust2" is feasible
#"""
#def CheckR1RelnConflict(ReachMat, clust1, clust2, outfile):
	#"""
	#explore the already existing R2 (in edge) list of the cluster "clust2"
	#since the cluster "clust1" is a possible candidate to be inserted in that list
	#"""
	#clust2_R2_list = Cluster_Info_Dict[clust2]._GetClustRelnList(RELATION_R2)
	#for cx in clust2_R2_list:
		#"""
		#check if there exists a cluster cx in clust2_R2_list
		#such that there exists no possibility of any directed edge between cx and clust1
		#i.e. cx->clust1 or cx<-clust1 is not possible by any means
		#"""
		#if (cx < clust1):
			#clust_pair_key = (cx, clust1)
		#else:
			#clust_pair_key = (clust1, cx)
		#if clust_pair_key in Cluster_Pair_Info_Dict:
			#cx_clust1_allowed_list = Cluster_Pair_Info_Dict[clust_pair_key]._GetPossibleRelnList()
		#else:
			#cx_clust1_allowed_list = []
		
		#"""
		#we check whether cx->clust1 or cx<-clust1 is already not established
		#and also the relations are not allowed by their respective cluster pair configuration
		#"""
		#if (ReachMat[CURRENT_CLUST_IDX_LIST.index(cx)][CURRENT_CLUST_IDX_LIST.index(clust1)] != 1) \
			#and (ReachMat[CURRENT_CLUST_IDX_LIST.index(clust1)][CURRENT_CLUST_IDX_LIST.index(cx)] != 1):
				#if (RELATION_R1 not in cx_clust1_allowed_list) and (RELATION_R2 not in cx_clust1_allowed_list):
					#if (DEBUG_LEVEL >= 2):
						#fp = open(outfile, 'a')
						#fp.write('\n Within function -- CheckR1RelnConflict :   target relation: ' + str(clust1) + '->' + str(clust2) + \
							#'  but there exists a cluster ' + str(cx) + '  such that there exists no directed edge possibility between them - so return False ')
						#fp.close()
					#return False
	
	#return True

##-------------------------------------------------------
#""" 
#this function checks whether the current relation (given as an input) between the input pair of cluster 
#is feasible to the existing DAG
#first we check whether the cluster pair is already connected
#otherwise, we check whether the connection transitively implies a connection (between some cluster pair) which is not allowed 
#"""
#def Check_Conflict(clust1, clust2, Reachability_Graph_Mat, target_reln_type, Output_Text_File):
	#"""
	#depending on the "target_reln_type" determine the src cluster from which a relation R1 to the dest cluster is sought
	#"""
	#if (target_reln_type == RELATION_R2):
		#src_taxa_clust_idx = clust2
		#dest_taxa_clust_idx = clust1
	#else:
		#src_taxa_clust_idx = clust1
		#dest_taxa_clust_idx = clust2
	
	##------------------------------------------------------------
	#"""
	#check whether the cluster pair is already resolved
	#"""
	#"""
	#first check whether clust1 = clust2 - same cluster
	#"""
	#if (src_taxa_clust_idx == dest_taxa_clust_idx):
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n already in same cluster index')
			#fp.close()
		#return 1
	
	#"""
	#we find the indices of the reachability matrix corresponding to individual cluster indices 
	#"""
	#reach_mat_src_taxa_clust_idx = CURRENT_CLUST_IDX_LIST.index(src_taxa_clust_idx)
	#reach_mat_dest_taxa_clust_idx = CURRENT_CLUST_IDX_LIST.index(dest_taxa_clust_idx)

	#"""
	#check whether R1 relation is already established
	#"""
	#if (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][reach_mat_dest_taxa_clust_idx] == 1):
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n target_reln_type: ' + str(target_reln_type) + ' already related via relation r1')
			#fp.close()
		#return 1

	#"""
	#check whether R2 relation is already established
	#"""
	#if (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][reach_mat_src_taxa_clust_idx] == 1):
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n target_reln_type: ' + str(target_reln_type) + ' already related via relation r2')
			#fp.close()
		#return 1
	
	##----------------------------------------------------------
	#"""
	#now check whether the proposed connection transitively implies certain connection which is not allowed 
	#between the corresponding cluster pair
	#"""
	
	#"""
	#if A->B is to be established
	#and there exists D->A 
	#then if B->D then return a conflict
	#"""
	#for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R2):
		#if (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n possible conflict - case 1 - target: ' + str(src_taxa_clust_idx) + '->' + str(dest_taxa_clust_idx) + \
					#' given ' + str(x) + '->' + str(src_taxa_clust_idx) + ' and ' + str(dest_taxa_clust_idx) + '->' + str(x))
				#fp.close()
			#return 1
	
	#"""
	#if A->B is to be established
	#and there exists B->E 
	#then if E->A then return a conflict  
	#"""
	#for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
		#if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] == 1):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n possible conflict - case 2 - target: ' + str(src_taxa_clust_idx) + '->' + str(dest_taxa_clust_idx) + \
					#' given ' + str(dest_taxa_clust_idx) + '->' + str(x) + ' and ' + str(x) + '->' + str(src_taxa_clust_idx))
				#fp.close()
			#return 1
	
	#"""
	#if A->B is to be established
	#and there exists D->A and B->E 
	#then if E->D or E=D then return a conflict  
	#"""
	#for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R2):
		#for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
			#if (x == y) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] == 1):
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n possible conflict - case 3 - target: ' + str(src_taxa_clust_idx) + '->' + str(dest_taxa_clust_idx) + \
					#' given ' + str(x) + '->' + str(src_taxa_clust_idx) + ' and ' \
						#+ str(dest_taxa_clust_idx) + '->' + str(y) + ' and ' + str(y) + '-> / =' + str(x))
					#fp.close()
				#return 1
	
	#"""
	#case 4 - for A->B connection
	#suppose D->A is already established
	#so, D->B will be implied
	#but if R1 is not allowed between D and B, reject A->B
	#"""
	#for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R2):
		#if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] == 0):
			#if (CheckAllowedRelnClusterPair(x, dest_taxa_clust_idx, RELATION_R1) == False):
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n possible conflict - case 4 - target: ' + str(src_taxa_clust_idx) + '->' + str(dest_taxa_clust_idx) + \
					#' given ' + str(x) + '->' + str(src_taxa_clust_idx) + ' but ' \
						#+ str(x) + '->' + str(dest_taxa_clust_idx) + '  is not allowed ')
					#fp.close()
				#return 1
	
	#"""
	#case 5 - for A->B connection
	#suppose B->E is already established
	#so, A->E will be implied
	#but if R1 is not allowed between A and E, reject A->B
	#"""
	#for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
		#if (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
			#if (CheckAllowedRelnClusterPair(src_taxa_clust_idx, x, RELATION_R1) == False):
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n possible conflict - case 4 - target: ' + str(src_taxa_clust_idx) + '->' + str(dest_taxa_clust_idx) + \
					#' given ' + str(dest_taxa_clust_idx) + '->' + str(x) + ' but ' \
						#+ str(src_taxa_clust_idx) + '->' + str(x) + '  is not allowed ')
					#fp.close()
				#return 1
	
	#"""
	#case 6 - for A->B connection
	#suppose D->A and B->E are already established
	#so, D->E will be implied
	#but if R1 is not allowed between D and E, reject A->B
	#"""
	#for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R2):
		#for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
			#if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 0):
				#if (CheckAllowedRelnClusterPair(x, y, RELATION_R1) == False):
					#if (DEBUG_LEVEL >= 2):
						#fp = open(Output_Text_File, 'a')
						#fp.write('\n possible conflict - case 4 - target: ' + str(src_taxa_clust_idx) + '->' + str(dest_taxa_clust_idx) + \
						#' given ' + str(x) + '->' + str(src_taxa_clust_idx) + ' and ' + \
							#str(dest_taxa_clust_idx) + '->' + str(y) + ' but ' \
							#+ str(x) + '->' + str(y) + '  is not allowed ')
						#fp.close()
					#return 1
	
	#return 0

#-------------------------------------------------------
""" 
this function checks whether the current relation (given as an input) is feasible to the input supertree configuration
if the current relation is already existing in the supertree, then this function returns 2
if the current relation produces a cycle to the existing configuration of the final supertree, then it returns 1
otherwise the function returns 0 
"""
def Possible_Conflict_Curr_Reln(clust1, clust2, Reachability_Graph_Mat, target_reln_type, Output_Text_File):
  
	if (target_reln_type == RELATION_R2):
		src_taxa_clust_idx = clust2
		dest_taxa_clust_idx = clust1
	else:
		src_taxa_clust_idx = clust1
		dest_taxa_clust_idx = clust2
	#--------------------------------------------------------------------
	"""
	check whether the cluster pair is already resolved
	"""
	"""
	first check whether clust1 = clust2 - same cluster
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
		return 2
	
	#--------------------------------------------------------------------
	# check for the transitive conflict
	if (target_reln_type == RELATION_R2) or (target_reln_type == RELATION_R1):
		
		"""
		if A->B is to be established
		and there exists D->A 
		then if B->D or B><D then return a conflict
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R2):
			if (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] > 0):	#== 1):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n possible conflict - case 1 - target: ' + str(src_taxa_clust_idx) + '->' + str(dest_taxa_clust_idx) + \
						' given ' + str(x) + '->' + str(src_taxa_clust_idx) + ' and ' + str(dest_taxa_clust_idx) + '-> / ><' + str(x))
					fp.close()
				return 1
		
		"""
		if A->B is to be established
		and there exists B->E 
		then if E->A or E><A then return a conflict  
		"""
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
			if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] > 0):	#== 1):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n possible conflict - case 2 - target: ' + str(src_taxa_clust_idx) + '->' + str(dest_taxa_clust_idx) + \
						' given ' + str(dest_taxa_clust_idx) + '->' + str(x) + ' and ' + str(x) + '-> / ><' + str(src_taxa_clust_idx))
					fp.close()
				return 1
		
		"""
		if A->B is to be established
		and there exists D->A and B->E 
		then if E->D or E=D or E><D then return a conflict  
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R2):
			for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
				if (x == y) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] > 0):	#== 1):
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n possible conflict - case 3 - target: ' + str(src_taxa_clust_idx) + '->' + str(dest_taxa_clust_idx) + \
						' given ' + str(x) + '->' + str(src_taxa_clust_idx) + ' and ' \
							+ str(dest_taxa_clust_idx) + '->' + str(y) + ' and ' + str(y) + '-> / = / ><' + str(x))
						fp.close()
					return 1
				
		"""
		if A->B is to be established
		and there exists D><A
		then it would be D><B afterwards 
		so if there is D->B, B->D or D=B then return a conflict  
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R4):
			if (dest_taxa_clust_idx == x) \
				or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] == 1) \
				or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1):
			#if (dest_taxa_clust_idx == x):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n possible conflict - case 4 - target: ' + str(src_taxa_clust_idx) + '->' + str(dest_taxa_clust_idx) + \
						' given ' + str(src_taxa_clust_idx) + '><' + str(x) + ' and ' + str(x) + '-> / = / <-' + str(dest_taxa_clust_idx))
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
				if (x == y) \
					or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
					or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 1):
				#if (x == y):
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n possible conflict - case 5 - target: ' + str(src_taxa_clust_idx) + '->' + str(dest_taxa_clust_idx) + \
						' given ' + str(src_taxa_clust_idx) + '><' + str(x) + \
							' and ' + str(dest_taxa_clust_idx) + '->' + str(y) + ' and ' + str(x) + '-> / = / <-' + str(y))
						fp.close()
					return 1
		
	elif (target_reln_type == RELATION_R3):
		
		"""
		if A=B is to be established
		and there exists B->D 
		then if D->A or D><A or D=A then return a conflict    
		"""
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
			if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] > 0):	#== 1):
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
			if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] > 0):	#== 1):
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
			if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
					or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] == 1):
			#if (x == src_taxa_clust_idx):
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
			if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] > 0):	#== 1):
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
			if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] > 0):	#== 1):
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
			if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
					or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] == 1):
			#if (x == dest_taxa_clust_idx):
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
				if (x == y) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
					or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 1):
				#if (x == y):
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n possible conflict - case 12')
						fp.close()
					return 1
		
	return 0
