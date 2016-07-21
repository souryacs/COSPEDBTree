#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *
import Conflict_Detect
from Conflict_Detect import *

#------------------------------------------------
"""
this function removes reln_type connection from clust1 to clust2
"""
def Remove_ClusterPairConn(clust1, clust2, reln_type):
	Cluster_Info_Dict[clust1]._RemoveRelnInstance(reln_type, clust2)
	Cluster_Info_Dict[clust2]._RemoveRelnInstance(Complementary_Reln(reln_type), clust1)
	return

#-------------------------------------------------
"""
checks the "score_list_cx" and finds the maximum / minimum, depending on the measure employed
"""
def Find_Clust_Idx_Max_Criterion(score_list_cx, measure):
	target_clust_idx = 0
	for i in range(1, len(score_list_cx)):
		if ((measure == 1) or (measure == 2)):
			if (score_list_cx[i][1] < score_list_cx[target_clust_idx][1]):
				target_clust_idx = i
		else:
			if (score_list_cx[i][1] > score_list_cx[target_clust_idx][1]):
				target_clust_idx = i
	
	return target_clust_idx

#-----------------------------------------------------        
"""
this function solves multiple parent problem (C2)
by uniquely selecting one particular parent
"""
def SelectUniqueParent_Directed(Output_Text_File, measure):
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')

	for cx in Cluster_Info_Dict:
		"""
		list of clusters cz such that cz->cx is satisfied
		"""
		cx_R2_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2)
		if (len(cx_R2_list) > 1):
			if (DEBUG_LEVEL >= 2):
				fp.write('\n\n ***** In function SelectUniqueParent_Directed ---- Examining cluster -- ' + str(cx))
			
			"""
			initialize one scoring list corresponding to the cluster cx
			the scoring list contains scores corresponding to the relations cz->cx for all possible clusters cz
			"""
			score_list_cx = []
			for cz in cx_R2_list:
				cz_score = FindAvgDistanceMeasure(Cluster_Info_Dict[cz]._GetSpeciesList(), \
					Cluster_Info_Dict[cx]._GetSpeciesList(), measure, 2, 1)	#0
				
				"""
				create a list containing the cluster cz information
				and also the support score and internode count measure
				"""
				temp_subl = [cz, cz_score]
				score_list_cx.append(temp_subl)
				if (DEBUG_LEVEL >= 2):
					fp.write('\n --- cluster (R2 reln): ' + str(cz) + '  Support score: ' + str(cz_score))
			
			"""
			determine the target cluster (according to the measure employed) which is the best
			"""
			target_score_clust_idx = Find_Clust_Idx_Max_Criterion(score_list_cx, measure)
			
			"""
			target cluster (best measure)
			"""
			target_score_clust = score_list_cx[target_score_clust_idx][0]
			
			if (DEBUG_LEVEL >= 2):
				fp.write('\n Initially selected target_score_clust_idx: ' + str(target_score_clust_idx) + \
					'  target_score_clust: ' + str(target_score_clust))
			
			#--------------------------------------------------------
			# this portion can be added - sourya
			#--------------------------------------------------------
			"""
			first we check whether there is a cluster "clust2"
			such that the R1 and R2 relation lists of clust2 are subsets of corresponding in target_score_clust
			"""
			subset_clust_list = []
			for i in range(len(score_list_cx)):
				if (i == target_score_clust_idx):
					continue
				if (CheckSubsetClust(target_score_clust, score_list_cx[i][0]) == True):
					templ = [score_list_cx[i][0], score_list_cx[i][1]]
					subset_clust_list.append(templ)
					if (DEBUG_LEVEL >= 2):
						fp.write('\n Here the cluster ' + str(score_list_cx[i][0]) + \
							' has R1 & R2 lists which are subsets of the corresponding in ' + str(target_score_clust) + \
								'  appending in the list of subset_clust_list')
			
			if (DEBUG_LEVEL >= 2):
				fp.write('\n Formed subset_clust_list: ' + str(subset_clust_list))
			
			"""
			if the subset_clust_list is not empty:
			sort the list according to the support score measure (descending order)
			and select the cluster with the highest support score as the parent of the current cluster
			"""
			if (len(subset_clust_list) > 0):
				if (measure == 1) or (measure == 2):
					subset_clust_list.sort(key=lambda x: x[1])
				else:
					subset_clust_list.sort(key=lambda x: x[1], reverse=True)
				target_clust_subset = subset_clust_list[0][0]
				if (DEBUG_LEVEL >= 2):
					fp.write('\n target_clust_subset: ' + str(target_clust_subset))

				"""
				first establish a directed edge connection from the cluster "target_score_clust"
				to the cluster "target_clust_subset"
				"""
				Cluster_Info_Dict[target_score_clust]._AddRelnInstance(RELATION_R1, target_clust_subset)
				Cluster_Info_Dict[target_clust_subset]._AddRelnInstance(RELATION_R2, target_score_clust)
				
				"""
				then remove all the clusters except the "target_clust_subset"
				from the parent list of cx
				"""
				for i in range(len(score_list_cx)):
					target_delete_clust = score_list_cx[i][0]
					if (target_delete_clust != target_clust_subset):
						Remove_ClusterPairConn(cx, target_delete_clust, RELATION_R2)
						if (DEBUG_LEVEL >= 2):
							fp.write('\n Removed connection ' + str(target_delete_clust) + '->' + str(cx))
			
				#--------------------------------------------------------
				# end
				#--------------------------------------------------------
			else:
				"""
				otherwise, we have already selected target_score_clust_idx and target_score_clust
				remove other clusters from the parent list
				"""
				for i in range(len(score_list_cx)):
					target_delete_clust = score_list_cx[i][0]
					if (target_delete_clust != target_score_clust):
						Remove_ClusterPairConn(cx, target_delete_clust, RELATION_R2)
						if (DEBUG_LEVEL >= 2):
							fp.write('\n Removed connection ' + str(target_delete_clust) + '->' + str(cx))

	# close the text file
	if (DEBUG_LEVEL >= 2):
		fp.close()
		
	return 

#-----------------------------------------------------        
"""
this function returns the root node for the final supertree 
for a depth first forest, multiple root nodes can be possible - 
so it returns the node with lowest indegree (or lowest R2 + possible R2 list)
"""
def Extract_Node_Min_Indeg(outfile, clust_dict):
	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')
		fp.write('\n\n In function --- Extract_Node_Min_Indeg \n\n ')
	"""
	stores the index of cluster having the minimum indegree
	"""
	min_indeg_node_idx = -1
	"""
	indicator that at least one cluster with minimum indegree (so far) has been found
	"""
	valid_node_found = 0
	
	for i in Cluster_Info_Dict:
		if (Cluster_Info_Dict[i]._GetExploredStatus() == 0):
			"""
			we check the clusters which have not been explored yet
			"""
			if (valid_node_found == 0):
				min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
				min_indeg_node_idx = i
				valid_node_found = 1
				if (DEBUG_LEVEL >= 2):
					fp.write('\n Minimum indegree cluster so far: ' + str(i))
			elif (valid_node_found == 1) and (Cluster_Info_Dict[i]._Get_Indegree() < min_indeg):
				min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
				min_indeg_node_idx = i
				if (DEBUG_LEVEL >= 2):
					fp.write('\n Minimum indegree cluster: ' + str(i))
			elif (valid_node_found == 1) and (Cluster_Info_Dict[i]._Get_Indegree() == min_indeg):
				if (Cluster_Info_Dict[i]._Get_Outdegree() > Cluster_Info_Dict[min_indeg_node_idx]._Get_Outdegree()):
					min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
					min_indeg_node_idx = i
					if (DEBUG_LEVEL >= 2):
						fp.write('\n Minimum indegree cluster: ' + str(i))
	
	return min_indeg_node_idx

#-----------------------------------------------------  
""" 
this function performs transitive reduction of a graph (transitive closure) 
and subsequently modifies the cluster of nodes
in terms of the edge connectivity, to make it free of redunant edges 
"""
def CompressDirectedGraph(Reachability_Graph_Mat, Outfile):
	no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
	
	# open the output text file
	if (DEBUG_LEVEL >= 2):
		fp = open(Outfile, 'a')
		fp.write('\n\n *** Inside the function  --- CompressDirectedGraph  **** \n\n')
	
	"""
	transitive reduction for the following case: 
	A->B, B->C, A->C ---------- remove A->C
	"""
	for j in range(no_of_clusters):
		clust_j = CURRENT_CLUST_IDX_LIST[j]
		for i in range(no_of_clusters):
			if (i == j):
				continue
			clust_i = CURRENT_CLUST_IDX_LIST[i]
			# A->B case
			if (Reachability_Graph_Mat[i][j] == 1):
				for k in range(no_of_clusters):
					if (i == k) or (j == k):
						continue
					clust_k = CURRENT_CLUST_IDX_LIST[k]
					# A->C and B->C case
					if (Reachability_Graph_Mat[j][k] == 1) and (Reachability_Graph_Mat[i][k] == 1):
						# we do not remove the matrix based connectivity information
						# Reachability_Graph_Mat[i][k] = 0
						
						# remove the edge from the cluster node directory
						Remove_ClusterPairConn(clust_i, clust_k, RELATION_R1)
						
						if (DEBUG_LEVEL >= 2):
							fp.write('\n ---- ' + str(clust_i) + ' -> ' + str(clust_j) + \
								' , ' + str(clust_j) + ' -> ' + str(clust_k) + ', and ' \
								+ str(clust_i) + ' -> ' + str(clust_k) + ' --- so removing ' + str(clust_i) + ' -> ' + str(clust_k))

	# close the output text file
	if (DEBUG_LEVEL >= 2):
		fp.write('\n\n *** Outside the function  --- CompressDirectedGraph  **** \n\n')
		fp.close()

	return

#-----------------------------------------------------
""" 
this function creates one new cluster with the given index value
also, it inserts one specified taxa in that cluster 
"""
def Create_Cluster_Taxa_Label(target_clust_idx, target_taxa_label):
	# create the cluster
	Cluster_Info_Dict.setdefault(target_clust_idx, Cluster_node(target_taxa_label))
	# include the cluster idx in the global list CURRENT_CLUST_IDX_LIST
	CURRENT_CLUST_IDX_LIST.append(target_clust_idx)
	# mention the cluster index in the taxa information
	taxa_key = COMPLETE_INPUT_TAXA_LIST.index(target_taxa_label)
	Taxa_Info_Dict[taxa_key]._Set_Clust_Idx_taxa_Part(target_clust_idx)

#-----------------------------------------------------
""" 
this function appends one specified taxon on a given cluster 
"""
def Append_Cluster_Taxa_Label(target_clust_idx, target_taxa_label):
	target_clust_spec_list = Cluster_Info_Dict[target_clust_idx]._GetSpeciesList()
	if target_taxa_label not in target_clust_spec_list:
		Cluster_Info_Dict[target_clust_idx]._Append_taxa(target_taxa_label)
		# mention the cluster index in the taxa information
		taxa_key = COMPLETE_INPUT_TAXA_LIST.index(target_taxa_label)
		Taxa_Info_Dict[taxa_key]._Set_Clust_Idx_taxa_Part(target_clust_idx)  

