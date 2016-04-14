#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#-----------------------------------------------------
""" 
this function computes the score of the relation R1 from clust1 to clust2
@parameters:
	clust1 and clust2 are taxa clusters containing one or more taxa
	MPP_SOLVE_METRIC: if 1, priority of relation R1 is used as the score measure
										if 2, excess gene leaf count (normalized) is used as the score measure
	DIST_MAT_TYPE: used when XL based score is used
									variation of XL measure is indicated by this parameter
	Both score measures are normalized by the number of support trees
"""
def ComputeScore(clust1, clust2, Output_Text_File, MPP_SOLVE_METRIC, DIST_MAT_TYPE):
  
	fp = open(Output_Text_File, 'a')
  
	if (DEBUG_LEVEL > 2):
		fp.write('\n ===>>> Score --- cluster 1 - taxa set: ' \
			+ str(Cluster_Info_Dict[clust1]._GetSpeciesList())\
			+ ' cluster 2 taxa set ' + str(Cluster_Info_Dict[clust2]._GetSpeciesList()))
		
	score_val = 0
	couplet_count = 0
	
	for t1 in Cluster_Info_Dict[clust1]._GetSpeciesList():
		t1_idx = COMPLETE_INPUT_TAXA_LIST.index(t1)
		for t2 in Cluster_Info_Dict[clust2]._GetSpeciesList():
			t2_idx = COMPLETE_INPUT_TAXA_LIST.index(t2)
			
			"""
			formation of the couplet key
			"""
			if (t1_idx < t2_idx):
				target_key = (t1_idx, t2_idx)
				target_reln = RELATION_R1
			else:
				target_key = (t2_idx, t1_idx)
				target_reln = RELATION_R2
			
			"""
			accumulate the score
			"""
			if target_key in TaxaPair_Reln_Dict:
				couplet_count = couplet_count + 1
				if (MPP_SOLVE_METRIC == 1):
					score_val = score_val + \
						(TaxaPair_Reln_Dict[target_key]._GetConnPrVal(target_reln) * 1.0) / TaxaPair_Reln_Dict[target_key]._GetNoSupportTrees()
				else:
					score_val = score_val + TaxaPair_Reln_Dict[target_key]._GetNormalizedXLSumGeneTrees(DIST_MAT_TYPE)
			else:
				if (DEBUG_LEVEL >= 2):
					fp.write('\n score compute -- key pair ' + str(t1) + ',' + str(t2) + ' does not exist ')

	if (DEBUG_LEVEL >= 2):
		fp.write('\n couplet_count: ' + str(couplet_count)) 
		if (couplet_count > 0):
			fp.write(' pairwise score of this cluster pair is : ' + str((score_val * 1.0) / couplet_count))
		
	fp.close()
	
	if (couplet_count > 0):
		return (score_val * 1.0) / couplet_count
	else:
		return 0

#-----------------------------------------------------    
""" 
this function solves multiple parent problem (C2)
by uniquely selecting one particular parent
the selection is carried out using a scoring mechanism 
@parameters:
	MPP_SOLVE_METRIC: if 1, priority of relation R1 is used as the score measure
										if 2, excess gene leaf count (normalized) is used as the score measure
	DIST_MAT_TYPE: used when XL based score is used
									variation of XL measure is indicated by this parameter
"""
def SolveMultipleParentC2Problem(Output_Text_File, MPP_SOLVE_METRIC, DIST_MAT_TYPE):
	# open the Output_Text_File
	fp = open(Output_Text_File, 'a')
	
	for cx in Cluster_Info_Dict:
		if (DEBUG_LEVEL > 2):
			fp.write('\n ***** Examining cluster -- ')
			Cluster_Info_Dict[cx]._PrintClusterInfo(cx, Output_Text_File)    
		
		if (Cluster_Info_Dict[cx]._Get_Indegree() > 1):
			"""
			for the current cluster cx
			take note of all its parent clusters (indexed by cz in the iterations)
			"""
			if (DEBUG_LEVEL > 2):
				fp.write('\n ***** Examining cluster with more than one indegree -- before in edge list fixing: ')
				Cluster_Info_Dict[cx]._PrintClusterInfo(cx, Output_Text_File)
			
			"""
			initialize one dictionary keyed by cluster indices cz
			cz signifies one parent node of the current cluster cx
			"""
			scoring_dict = dict()
			for cz in Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2):
				scoring_dict.setdefault(cz, 0)
			
			"""
			now for each of the parent clusters cz of the current cluster cx, 
			compute the R1 score from cz to cx
			"""
			for cz in Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2):
				scoring_dict[cz] = ComputeScore(cz, cx, Output_Text_File, MPP_SOLVE_METRIC, DIST_MAT_TYPE)
			
			"""
			all such R1 scores from cz to cx (where cz iteratively points to one ancestor of cx)
			are saved in a list named "Scoring_List"
			items of this list is a key-value pair
			key: the cluster cz
			value: the score measure
			"""
			Scoring_List = []
			for cz in Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2):
				if (DEBUG_LEVEL > 2):
					fp.write('\n scoring dict elem: ' + str(cz) + ' score: ' + str(Scoring_Dict[cz]))
				temp_subl = [cz, scoring_dict[cz]]
				Scoring_List.append(temp_subl)
			
			if (DEBUG_LEVEL > 2):
				fp.write('\n --- before sorting the scoring list --- ')
				for i in range(len(Scoring_List)):
					fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(Scoring_List[i][0]) + ' score: ' + str(Scoring_List[i][1]))

			"""
			sort the scoring list in ascending order
			"""
			Scoring_List.sort(key=lambda x: x[1])
			if (DEBUG_LEVEL > 2):
				fp.write('\n --- after sorting the scoring list --- ')
				for i in range(len(Scoring_List)):
					fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(Scoring_List[i][0]) + ' score: ' + str(Scoring_List[i][1]))
		
			if (MPP_SOLVE_METRIC == 1):
				"""
				for the priority measure, remove all except the last element from the 
				scoring list of the current cluster cx
				the last element contains the highest priority measure
				"""
				for i in range(len(Scoring_List) - 1):
					target_delete_clust_idx = Scoring_List[i][0]
					Cluster_Info_Dict[cx]._RemoveRelnInstance(RELATION_R2, target_delete_clust_idx)
					Cluster_Info_Dict[target_delete_clust_idx]._RemoveRelnInstance(RELATION_R1, cx)
			else:
				"""
				for the XL based measure, remove all except the first element from the 
				scoring list of the current cluster cx
				the first element contains the lowest XL measure
				"""
				for i in range(1, len(Scoring_List)):
					target_delete_clust_idx = Scoring_List[i][0]
					Cluster_Info_Dict[cx]._RemoveRelnInstance(RELATION_R2, target_delete_clust_idx)
					Cluster_Info_Dict[target_delete_clust_idx]._RemoveRelnInstance(RELATION_R1, cx)
			
			if (DEBUG_LEVEL > 2):
				fp.write('\n ***** Examining cluster with more than one indegree -- after in edge list fixing: ')
				Cluster_Info_Dict[cx]._PrintClusterInfo(cx, Output_Text_File)

	# close the Output_Text_File
	fp.close()
	return

#-----------------------------------------------------        
"""
this function returns the root node for the final supertree 
for a depth first forest, multiple root nodes can be possible - 
so it returns the node with 0 indegree 
"""
def Extract_Node_Min_Indeg(no_of_clusters):
	min_indeg_node_idx = -1
	valid_node_found = 0
	for i in Cluster_Info_Dict:
		if (Cluster_Info_Dict[i]._GetExploredStatus() == 0):
			if (valid_node_found == 0):
				min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
				min_indeg_node_idx = i
				valid_node_found = 1
			elif (valid_node_found == 1) and (Cluster_Info_Dict[i]._Get_Indegree() < min_indeg):
				min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
				min_indeg_node_idx = i
			elif (valid_node_found == 1) and (Cluster_Info_Dict[i]._Get_Indegree() == min_indeg)\
				and (Cluster_Info_Dict[i]._Get_Outdegree() > Cluster_Info_Dict[min_indeg_node_idx]._Get_Outdegree()):    
				min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
				min_indeg_node_idx = i
		
	return min_indeg_node_idx

#-----------------------------------------------------  
""" 
this function performs transitive reduction of a graph (transitive closure) 
and subsequently modifies the cluster of nodes
in terms of the edge connectivity, to make it free of redunant edges 
"""
def CompressDirectedGraph(Reachability_Graph_Mat):
	no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
	
	# transitive reduction
	for j in range(no_of_clusters):
		for i in range(no_of_clusters):
			# A->B case
			if (Reachability_Graph_Mat[i][j] == 1):
				for k in range(no_of_clusters):
					# A->C and B->C case
					if (Reachability_Graph_Mat[j][k] == 1) and (Reachability_Graph_Mat[i][k] == 1):
						# comment - sourya - check
						#Reachability_Graph_Mat[i][k] = 0
						
						# remove the edge from the cluster node directory
						clust_i = CURRENT_CLUST_IDX_LIST[i]
						clust_k = CURRENT_CLUST_IDX_LIST[k]
						Cluster_Info_Dict[clust_i]._RemoveRelnInstance(RELATION_R1, clust_k)
						Cluster_Info_Dict[clust_k]._RemoveRelnInstance(RELATION_R2, clust_i)

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
	if target_taxa_label not in Cluster_Info_Dict[target_clust_idx]._GetSpeciesList():
		Cluster_Info_Dict[target_clust_idx]._Append_taxa(target_taxa_label)
		# mention the cluster index in the taxa information
		taxa_key = COMPLETE_INPUT_TAXA_LIST.index(target_taxa_label)
		Taxa_Info_Dict[taxa_key]._Set_Clust_Idx_taxa_Part(target_clust_idx)  

