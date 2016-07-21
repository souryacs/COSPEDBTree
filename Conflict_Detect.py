#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

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
