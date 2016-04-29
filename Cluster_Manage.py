#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#------------------------------------------------
"""
this function removes reln_type connection from clust1 to clust2
"""
def Remove_ClusterPairConn(clust1, clust2, reln_type):
	Cluster_Info_Dict[clust1]._RemoveRelnInstance(reln_type, clust2)
	Cluster_Info_Dict[clust2]._RemoveRelnInstance(Complementary_Reln(reln_type), clust1)
	return

#----------------------------------------------------
"""
this function establishes a possible R1 connection from clust1 to clust2
"""
def Connect_PossibleR1Reln_ClusterPair(clust1, clust2):
	Cluster_Info_Dict[clust1]._AddPossibleR1(clust2)
	Cluster_Info_Dict[clust2]._AddPossibleR2(clust1)
	return

#-----------------------------------------------------
# add - sourya
"""
this function checks and connects two cluster pairs with pseudo R1 / R2 connectivity
"""
def CheckPseudoR1RelnClust(clust1, clust2):
	r1_reln = False
	r2_reln = False
	
	for t1 in Cluster_Info_Dict[clust1]._GetSpeciesList():
		t1_idx = COMPLETE_INPUT_TAXA_LIST.index(t1)
		for t2 in Cluster_Info_Dict[clust2]._GetSpeciesList():
			t2_idx = COMPLETE_INPUT_TAXA_LIST.index(t2)
			"""
			formation of the couplet key
			"""
			if (t1_idx < t2_idx):
				target_key = (t1_idx, t2_idx)
				complement_operation = False
			else:
				target_key = (t2_idx, t1_idx)
				complement_operation = True
				
			if target_key in TaxaPair_Reln_Dict:
				supp_tree = TaxaPair_Reln_Dict[target_key]._GetNoSupportTrees()
				fr4 = TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R4)
				pr3 = TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(2)
				if (complement_operation == False):
					pr1 = TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(0)
					pr2 = TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(1)
				else:
					pr1 = TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(1)
					pr2 = TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(0)
					
				"""
				sourya - previously the number of support trees were restricted as > 1
				now we have included the condition >= 1 as well
				"""
				if (supp_tree >= 1) and (fr4 == (pr1 + pr3)):
					r1_reln = True
					
				if (supp_tree >= 1) and (fr4 == (pr2 + pr3)):
					r2_reln = True
					
				if (r1_reln == True) and (r2_reln == True):
					break
	
	"""
	we add pseudo R1 / R2 edges, depending on whether it satisfies the given relations
	"""
	if (r1_reln == True):
		Connect_PossibleR1Reln_ClusterPair(clust1, clust2)

	if (r2_reln == True):
		Connect_PossibleR1Reln_ClusterPair(clust2, clust1)
		
	return

# end add - sourya
#-----------------------------------------------------
""" 
this function adds an edge between a pair of clusters (of taxa) 
it also updates the entries of reachability matrix 
"""
def Connect_ClusterPair(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, reln_type, nodeA_clust_idx, nodeB_clust_idx):
	if (reln_type == RELATION_R1):
		"""
		adjust the clusters
		"""
		Cluster_Info_Dict[nodeA_clust_idx]._AddRelnInstance(RELATION_R1, nodeB_clust_idx)
		Cluster_Info_Dict[nodeB_clust_idx]._AddRelnInstance(RELATION_R2, nodeA_clust_idx)
		"""
		update the reachability matrix
		"""
		Reachability_Graph_Mat[nodeA_reach_mat_idx][nodeB_reach_mat_idx] = 1
		#Reachability_Graph_Mat[nodeB_reach_mat_idx][nodeA_reach_mat_idx] = -1	# add - sourya
	elif (reln_type == RELATION_R4):
		"""
		adjust the clusters
		"""
		Cluster_Info_Dict[nodeA_clust_idx]._AddRelnInstance(RELATION_R4, nodeB_clust_idx)
		Cluster_Info_Dict[nodeB_clust_idx]._AddRelnInstance(RELATION_R4, nodeA_clust_idx)    
		"""
		update the reachability matrix
		"""
		Reachability_Graph_Mat[nodeA_reach_mat_idx][nodeB_reach_mat_idx] = 2
		Reachability_Graph_Mat[nodeB_reach_mat_idx][nodeA_reach_mat_idx] = 2
		# add - sourya
		"""
		we check whether pseudo R1 relation can be established between these clusters
		"""
		CheckPseudoR1RelnClust(nodeA_clust_idx, nodeB_clust_idx)
		# end add - sourya

#-----------------------------------------------------        
"""
solves multiple parent problem from the possible R1 / R2 lists 
the selection is carried out using one of the following three measures:
	1) Priority of R1 relation
	2) XL based 
	3) internode count based scoring mechanism
This function is applicable for directed in edge (Relation R2) only
"""
def Solve_MPP_PossibleR1R2(Reachability_Graph_Mat, Output_Text_File):
	#--------------------------------------------------------------
	"""
	now we check the clusters cx having following properties:
	1) No cluster cy exists such that cy -> cx holds
	2) There exists at least two clusters c1 and c2 such that c1--->cx and c2--->cx
	in such a case, we select among the multiple dashed edges (possible R2 lists), the suitable one or more possible R2 candidate
	"""
	for cx in Cluster_Info_Dict:
		"""
		explore only the clusters having no in degree (no cluster is connected with R2 relation)
		"""
		if 1:	#(Cluster_Info_Dict[cx]._Get_Indegree() == 0):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n ***** Examining cluster with zero indegree -- ' + str(cx))
				#fp.close()
			#if (DEBUG_LEVEL > 2):
				#Cluster_Info_Dict[cx]._PrintClusterInfo(cx, Output_Text_File)    
			
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n ***** In function - Solve_MPP_PossibleR1R2 - Examining cluster -- ' + str(cx))
				fp.close()
			
			# open the text file
			fp = open(Output_Text_File, 'a')

			"""
			explore all clusters belonging to the possible R2 list of the cluster cx
			that is, for each cluster cz in the possible R2 list, cz--->cx relation is established with respect to input gene trees
			"""
			cx_possible_R2_list = Cluster_Info_Dict[cx]._GetPossibleR2List()
			"""
			first we check only those clusters (cx) such that the length of cx_possible_R2_list is > 0, 
			in such a case, we have at least one cluster cz such that cz--->cx holds
			"""
			if (len(cx_possible_R2_list) > 1):
				"""
				initialize one scoring list corresponding to the cluster cx
				"""
				score_list_cx = []
				"""
				explore all candidate clusters cz belonging to the possible R2 list of the cluster cx
				"""
				for cz in cx_possible_R2_list:
					#cz_score1 = AvgR1RelnPriority(Cluster_Info_Dict[cz]._GetSpeciesList(), Cluster_Info_Dict[cx]._GetSpeciesList())
					#cz_score2 = FindAvgXL(Cluster_Info_Dict[cz]._GetSpeciesList(), Cluster_Info_Dict[cx]._GetSpeciesList(), DIST_MAT_TYPE, 2, 1)
					cz_score3 = FindAvgInternodeCount(Cluster_Info_Dict[cz]._GetSpeciesList(), Cluster_Info_Dict[cx]._GetSpeciesList(), 2, 1)
					if (DEBUG_LEVEL > 2):
						#fp.write('\n --- element (R2 reln): ' + str(cz) + ' R1 reln priority score: ' + str(cz_score1) + \
							#' XL score: ' + str(cz_score2) + ' Internode count score: ' + str(cz_score3))
						fp.write('\n --- element (R2 reln): ' + str(cz) + ' Internode count score: ' + str(cz_score3))
						
					#temp_subl = [cz, cz_score1, cz_score2, cz_score3]
					temp_subl = [cz, cz_score3]
					score_list_cx.append(temp_subl)
			
				"""
				sort the scoring list in ascending order, with respect to the internode count
				"""
				#score_list_cx.sort(key=lambda x: x[3])
				score_list_cx.sort(key=lambda x: x[1])
				if (DEBUG_LEVEL >= 2):
					fp.write('\n --- after sorting the scoring list corresponding to the cluster : ' + str(cx))
					for i in range(len(score_list_cx)):
						#fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(score_list_cx[i][0]) + \
							#' R1 reln priority score: ' + str(score_list_cx[i][1]) + \
								#' XL score: ' + str(score_list_cx[i][2]) + ' Internode count score: ' + str(score_list_cx[i][3]))
						fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(score_list_cx[i][0]) + \
							' Internode count score: ' + str(score_list_cx[i][1]))
			
				"""
				we note the min internode count measure according to the sorted list
				we only delete the clusters (from the possible R2 list) if they have internode count greater than this minimum value
				"""
				#min_internode_count = score_list_cx[0][3]
				min_internode_count = score_list_cx[0][1]
				for i in range(1, len(score_list_cx)):
					#if (score_list_cx[i][3] > min_internode_count):
					if (score_list_cx[i][1] > min_internode_count):
						# delete this entry from the possible R2 list
						target_delete_clust_idx = score_list_cx[i][0]
						Cluster_Info_Dict[cx]._RemovePossibleR2(target_delete_clust_idx)
						Cluster_Info_Dict[target_delete_clust_idx]._RemovePossibleR1(cx)
						if (DEBUG_LEVEL >= 2):
							fp.write('\n Removed possible R1 edge from the cluster ' + str(target_delete_clust_idx) + '  to the cluster: ' + str(cx))
			
			# close the text file
			fp.close()

	return

#-----------------------------------------------------        
"""
solves no parent problem - checks the clusters having no parent at all 
(empty R2 relation list) and empty possible R2 relation list as well
in such a case, it finds a candidate cluster which will be placed as a R2 or 
possible R2 relation candidate
"""
def Solve_NPP_NoPossibleR2(Reachability_Graph_Mat, Output_Text_File):
	
	#--------------------------------------------------------------
	"""
	first we check the clusters cx having following properties:
	1) No cluster cy exists such that cy -> cx holds
	2) There exists no clusters c1 such that c1--->cx 
	"""
	for cx in Cluster_Info_Dict:
		"""
		clusters having no indegree (no cluster is connected with R2 relation)
		"""
		if (Cluster_Info_Dict[cx]._Get_Indegree() == 0):
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n ***** Examining cluster with zero indegree -- ' + str(cx))
				fp.close()
			if (DEBUG_LEVEL > 2):
				Cluster_Info_Dict[cx]._PrintClusterInfo(cx, Output_Text_File)    
			
			"""
			possible R2 list of the cluster cx --- that is, collection of cz such that cz--->cx holds
			"""
			cx_possible_R2_list = Cluster_Info_Dict[cx]._GetPossibleR2List()
			"""
			here we check only those cases such that no cz exists
			"""
			if (len(cx_possible_R2_list) == 0):

				# open the text file
				fp = open(Output_Text_File, 'a')
				
				if (DEBUG_LEVEL >= 2):
					fp.write('\n ***** The cluster has no element in possible R2 list as well -- ')
				
				"""
				it is a boolean flag, indicating that cx has no in edge or possible R2 cluster
				"""
				flag = False
				
				"""
				first find clusters x having the following property:
				let cx->y (Relation R1, if holds)
				and x->y (parent of y, excluding cx)
				for such x, if x and cx are not related by any means, we include x->cx
				"""
				cx_R1_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R1)
				if (len(cx_R1_list) > 0):
					for y in cx_R1_list:
						y_R2_list = Cluster_Info_Dict[y]._GetClustRelnList(RELATION_R2)
						if (len(y_R2_list) > 0):
							for x in y_R2_list:
								if (x != cx):
									if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(cx)] == 0) and \
										(Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(cx)][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
										if (DEBUG_LEVEL >= 2):
											fp.write('\n We have ' + str(cx) + '->' + str(y) + '  and  ' + str(x) + '->' + str(y) + \
												'  ===>>> So inserting ' + str(x) + '-> ' + str(cx) + '  since they are not related ')
										# establish x->cx
										Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
											CURRENT_CLUST_IDX_LIST.index(cx), RELATION_R1, x, cx)
										# set the flag as well
										flag = True
				
				"""
				if the flag is True at this point, then some cluster is assigned as the parent of cx
				otherwise we explore the clusters x having the following property:
				let cx--->y (Relation R1, if holds)
				and x->y (parent of y, excluding cx)
				for such x, if x and cx are not related by any means, we include x->cx
				"""
				if (flag == False):
					cx_possible_R1_list = Cluster_Info_Dict[cx]._GetPossibleR1List()
					if (len(cx_possible_R1_list) > 0):
						for y in cx_possible_R1_list:
							y_R2_list = Cluster_Info_Dict[y]._GetClustRelnList(RELATION_R2)
							if (len(y_R2_list) > 0):
								for x in y_R2_list:
									if (x != cx):
										if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(cx)] == 0) and \
											(Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(cx)][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
											# establish x->cx
											if (DEBUG_LEVEL >= 2):
												fp.write('\n We have ' + str(cx) + '-->' + str(y) + '  and  ' + str(x) + '->' + str(y) + \
													'  ===>>> So inserting ' + str(x) + '-> ' + str(cx) + '  since they are not related ')
											# establish x->cx
											Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
												CURRENT_CLUST_IDX_LIST.index(cx), RELATION_R1, x, cx)
											flag = True
				
				"""
				if the flag is True at this point, then some cluster is assigned as the parent of cx
				otherwise, again we explore all the clusters such that cx--->y
				this time we explore the clusters such that x--->y (and x is not cx)
				we establish the relation x--->cx
				"""
				if (flag == False):
					cx_possible_R1_list = Cluster_Info_Dict[cx]._GetPossibleR1List()
					for y in cx_possible_R1_list:
						y_possible_R2_list = Cluster_Info_Dict[y]._GetPossibleR2List()
						if (len(y_possible_R2_list) > 0):
							for x in y_possible_R2_list:
								if (x != cx):
									if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(cx)] == 0) and \
										(Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(cx)][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
										if (DEBUG_LEVEL >= 2):
											fp.write('\n We have ' + str(cx) + '-->' + str(y) + '  and  ' + str(x) + '--->' + str(y) + \
												'  ===>>> So inserting ' + str(x) + '----> ' + str(cx) + '  since they are not related ')
										# establish x--->cx
										Connect_PossibleR1Reln_ClusterPair(x, cx)
										flag = True
	
	
				# close the text file
				fp.close()

	return

#-----------------------------------------------------        
"""
this function solves multiple parent problem (C2)
by uniquely selecting one particular parent
the selection is carried out using one of the following three measures:
	1) Priority of R1 relation
	2) XL based 
	3) internode count based scoring mechanism
This function is applicable for directed in edge (Relation R2) only
"""
def SelectUniqueParent_Directed(Output_Text_File):

	for cx in Cluster_Info_Dict:
		if (DEBUG_LEVEL > 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n ***** Examining cluster -- ' + str(cx))
			fp.close()
			Cluster_Info_Dict[cx]._PrintClusterInfo(cx, Output_Text_File)    

		# open the text file
		fp = open(Output_Text_File, 'a')
		
		"""
		list of clusters belonging to the R2 list of cx
		"""
		cx_R2_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2)
		if (len(cx_R2_list) > 1):
			
			"""
			initialize one scoring list corresponding to the cluster cx
			"""
			score_list_cx = []
			
			"""
			explore all the in edges (relation R2)
			"""
			for cz in cx_R2_list:
				cz_score = FindAvgInternodeCount(Cluster_Info_Dict[cz]._GetSpeciesList(), Cluster_Info_Dict[cx]._GetSpeciesList(), 2, 1)
				temp_subl = [cz, cz_score]
				score_list_cx.append(temp_subl)
				if (DEBUG_LEVEL >= 2):
					fp.write('\n --- element (R2 reln): ' + str(cz) + ' Internode score: ' + str(cz_score))
					
			"""
			sort the scoring list in descending order
			for internode count based measure, select the cluster having the highest average internode count
			"""
			score_list_cx.sort(key=lambda x: x[1], reverse=True)
			if (DEBUG_LEVEL >= 2):
				fp.write('\n --- after sorting the scoring list corresponding to the cluster : ' + str(cx))
				for i in range(len(score_list_cx)):
					fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(score_list_cx[i][0]) + ' score: ' + str(score_list_cx[i][1]))

			"""
			for the internode count based measure, remove all except the first element from the 
			scoring list of the current cluster cx
			since the first element contains the highest internode count measure
			"""
			max_internode_measure = score_list_cx[0][1]
			for i in range(1, len(score_list_cx)):
				if (score_list_cx[i][1] < max_internode_measure):
					target_delete_clust_idx = score_list_cx[i][0]
					Cluster_Info_Dict[cx]._RemoveRelnInstance(RELATION_R2, target_delete_clust_idx)
					Cluster_Info_Dict[target_delete_clust_idx]._RemoveRelnInstance(RELATION_R1, cx)
					if (DEBUG_LEVEL >= 2):
						fp.write('\n Removed out edge from the cluster ' + str(target_delete_clust_idx) + '  to the cluster: ' + str(cx))
		
		"""
		check whether the R2 list is now of length 1
		otherwise, analyze the clusters again
		"""
		cx_R2_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2)
		if (len(cx_R2_list) > 1):
			"""
			initialize one scoring list corresponding to the cluster cx
			"""
			score_list_cx = []
			"""
			explore all the in edges (relation R2)
			"""
			for cz in cx_R2_list:
				cz_score = len(Cluster_Info_Dict[cz]._GetClustRelnList(RELATION_R1))
				temp_subl = [cz, cz_score]
				score_list_cx.append(temp_subl)
				if (DEBUG_LEVEL >= 2):
					fp.write('\n --- element (R2 reln): ' + str(cz) + ' No of R1 reln clusters: ' + str(cz_score))

			"""
			sort the scoring list in ascending order
			for internode count based measure, select the cluster having the highest average internode count
			"""
			score_list_cx.sort(key=lambda x: x[1])
			if (DEBUG_LEVEL >= 2):
				fp.write('\n --- after sorting the scoring list corresponding to the cluster : ' + str(cx))
				for i in range(len(score_list_cx)):
					fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(score_list_cx[i][0]) + ' score: ' + str(score_list_cx[i][1]))

			"""
			remove all except the first element from the 
			scoring list of the current cluster cx
			"""
			for i in range(1, len(score_list_cx)):
				target_delete_clust_idx = score_list_cx[i][0]
				Cluster_Info_Dict[cx]._RemoveRelnInstance(RELATION_R2, target_delete_clust_idx)
				Cluster_Info_Dict[target_delete_clust_idx]._RemoveRelnInstance(RELATION_R1, cx)
				if (DEBUG_LEVEL >= 2):
					fp.write('\n Removed out edge from the cluster ' + str(target_delete_clust_idx) + '  to the cluster: ' + str(cx))

		# close the text file
		fp.close()
		
	return

#-----------------------------------------------------        
"""
this function returns the root node for the final supertree 
for a depth first forest, multiple root nodes can be possible - 
so it returns the node with lowest indegree (or lowest R2 + possible R2 list)
"""
def Extract_Node_Min_Indeg(outfile):
	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')
		fp.write('\n\n In function --- Extract_Node_Min_Indeg \n\n ')
	
	""" 
	total number of clusters
	"""
	no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
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
			# we check the clusters which have not been explored yet
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
				if (len(Cluster_Info_Dict[i]._GetDistinctPossibleR2List()) < len(Cluster_Info_Dict[min_indeg_node_idx]._GetDistinctPossibleR2List())):
					min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
					min_indeg_node_idx = i
					if (DEBUG_LEVEL >= 2):
						fp.write('\n Minimum indegree cluster: ' + str(i))
				elif (len(Cluster_Info_Dict[i]._GetDistinctPossibleR2List()) == len(Cluster_Info_Dict[min_indeg_node_idx]._GetDistinctPossibleR2List())) \
					and ((Cluster_Info_Dict[i]._Get_Outdegree() + len(Cluster_Info_Dict[i]._GetPossibleR1List())) > \
						(Cluster_Info_Dict[min_indeg_node_idx]._Get_Outdegree() + len(Cluster_Info_Dict[min_indeg_node_idx]._GetPossibleR1List()))):    
					min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
					min_indeg_node_idx = i
					if (DEBUG_LEVEL >= 2):
						fp.write('\n Minimum indegree cluster: ' + str(i))
	
	if (DEBUG_LEVEL >= 2):
		fp.write('\n\n Out of the function --- Extract_Node_Min_Indeg \n\n ')
		fp.close()
	
	return min_indeg_node_idx

#-----------------------------------------------------  
""" 
this function performs transitive reduction of a graph for the dashed edges 
"""
def CompressDAG_Dashed(Reachability_Graph_Mat, Clust_Possible_R1_Mat, Outfile):
	no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
	
	# open the output text file
	if (DEBUG_LEVEL >= 2):
		fp = open(Outfile, 'a')

	"""
	transitive reduction for the following case: 
	A--->B, B---->C, A---->C ---------- remove A---->C
	Note: we should not consider the case when A--->B and B--->A simultaneously
	"""
	for j in range(no_of_clusters):
		clust_j = CURRENT_CLUST_IDX_LIST[j]
		for i in range(no_of_clusters):
			if (i == j):
				continue
			clust_i = CURRENT_CLUST_IDX_LIST[i]
			# A--->B case
			if (Clust_Possible_R1_Mat[i][j] == 1) and (Clust_Possible_R1_Mat[j][i] == 0):
				for k in range(no_of_clusters):
					if (i == k) or (j == k):
						continue
					clust_k = CURRENT_CLUST_IDX_LIST[k]
					# A--->C and B--->C case
					if (Clust_Possible_R1_Mat[j][k] == 1) and (Clust_Possible_R1_Mat[k][j] == 0) \
						and (Clust_Possible_R1_Mat[i][k] == 1) and (Clust_Possible_R1_Mat[k][i] == 0):
						# comment - sourya - check
						#Clust_Possible_R1_Mat[i][k] = 0
						
						Cluster_Info_Dict[clust_i]._RemovePossibleR1(clust_k)
						Cluster_Info_Dict[clust_k]._RemovePossibleR2(clust_i)
						if (DEBUG_LEVEL >= 2):
							fp.write('\n ---- ' + str(clust_i) + ' ----> ' + str(clust_j) + ', ' + str(clust_j) + ' ----> ' + str(clust_k) + ', and ' \
								+ str(clust_i) + ' -----> ' + str(clust_k) + ' --- so removing ' + str(clust_i) + ' ----> ' + str(clust_k))


	"""
	transitive reduction for the following case: 
	A->B, B---->C, A->C ---------- remove A->C
	"""
	for j in range(no_of_clusters):
		clust_j = CURRENT_CLUST_IDX_LIST[j]
		for i in range(no_of_clusters):
			if (i == j):
				continue
			clust_i = CURRENT_CLUST_IDX_LIST[i]
			#if (Reachability_Graph_Mat[i][j] == 1):
			if clust_j in Cluster_Info_Dict[clust_i]._GetClustRelnList(RELATION_R1):
				for k in range(no_of_clusters):
					if (i == k) or (j == k):
						continue
					clust_k = CURRENT_CLUST_IDX_LIST[k]
					#if (Clust_Possible_R1_Mat[j][k] == 1) and (Clust_Possible_R1_Mat[k][j] == 0) and (Reachability_Graph_Mat[i][k] == 1):
					if clust_k in Cluster_Info_Dict[clust_j]._GetPossibleR1List():
						if clust_j not in Cluster_Info_Dict[clust_k]._GetPossibleR1List():
							if clust_k in Cluster_Info_Dict[clust_i]._GetClustRelnList(RELATION_R1):
						
								Remove_ClusterPairConn(clust_i, clust_k, RELATION_R1)
								
								if (DEBUG_LEVEL >= 2):
									fp.write('\n ---- ' + str(clust_i) + ' -> ' + str(clust_j) + ', ' + str(clust_j) + ' ----> ' + str(clust_k) + ', and ' \
										+ str(clust_i) + ' -> ' + str(clust_k) + ' --- so removing ' + str(clust_i) + ' -> ' + str(clust_k))


	"""
	transitive reduction for the following case: 
	A---->B, B->C, A->C ---------- remove A->C
	"""
	for j in range(no_of_clusters):
		clust_j = CURRENT_CLUST_IDX_LIST[j]
		for i in range(no_of_clusters):
			if (i == j):
				continue
			clust_i = CURRENT_CLUST_IDX_LIST[i]
			if clust_j in Cluster_Info_Dict[clust_i]._GetPossibleR1List():
				if clust_i not in Cluster_Info_Dict[clust_j]._GetPossibleR1List():
					for k in range(no_of_clusters):
						if (i == k) or (j == k):
							continue
						clust_k = CURRENT_CLUST_IDX_LIST[k]
						#if (Reachability_Graph_Mat[j][k] == 1) and (Reachability_Graph_Mat[i][k] == 1):
						if clust_k in Cluster_Info_Dict[clust_j]._GetClustRelnList(RELATION_R1):
							if clust_k in Cluster_Info_Dict[clust_i]._GetClustRelnList(RELATION_R1):
								
								Remove_ClusterPairConn(clust_i, clust_k, RELATION_R1)
								
								if (DEBUG_LEVEL >= 2):
									fp.write('\n ---- ' + str(clust_i) + ' ----> ' + str(clust_j) + ', ' + str(clust_j) + ' -> ' + str(clust_k) + ', and ' \
										+ str(clust_i) + ' -> ' + str(clust_k) + ' --- so removing ' + str(clust_i) + ' -> ' + str(clust_k))


	# close the output text file
	if (DEBUG_LEVEL >= 2):
		fp.close()

	return

#-----------------------------------------------------  
""" 
this function performs transitive reduction of a graph (transitive closure) 
and subsequently modifies the cluster of nodes
in terms of the edge connectivity, to make it free of redunant edges 
"""
def CompressDirectedGraph(Reachability_Graph_Mat, Clust_Possible_R1_Mat, Outfile):
	no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
	
	# open the output text file
	if (DEBUG_LEVEL >= 2):
		fp = open(Outfile, 'a')
	
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
						# comment - sourya - check
						#Reachability_Graph_Mat[i][k] = 0
						
						# remove the edge from the cluster node directory
						Remove_ClusterPairConn(clust_i, clust_k, RELATION_R1)
						
						if (DEBUG_LEVEL >= 2):
							fp.write('\n ---- ' + str(clust_i) + ' -> ' + str(clust_j) + \
								' , ' + str(clust_j) + ' -> ' + str(clust_k) + ', and ' \
								+ str(clust_i) + ' -> ' + str(clust_k) + ' --- so removing ' + str(clust_i) + ' -> ' + str(clust_k))


	#"""
	#transitive reduction for the following case: 
	#A---->B, B---->C, A->C ---------- remove A->C
	#"""
	#for j in range(no_of_clusters):
		#clust_j = CURRENT_CLUST_IDX_LIST[j]
		#for i in range(no_of_clusters):
			#if (i == j):
				#continue
			#clust_i = CURRENT_CLUST_IDX_LIST[i]
			#if clust_j in Cluster_Info_Dict[clust_i]._GetPossibleR1List():
				#for k in range(no_of_clusters):
					#if (i == k) or (j == k):
						#continue
					#clust_k = CURRENT_CLUST_IDX_LIST[k]
					#if (clust_k in Cluster_Info_Dict[clust_j]._GetPossibleR1List()) and (Reachability_Graph_Mat[i][k] == 1):
						#Cluster_Info_Dict[clust_i]._RemoveRelnInstance(RELATION_R1, clust_k)
						#Cluster_Info_Dict[clust_k]._RemoveRelnInstance(RELATION_R2, clust_i)
						#if (DEBUG_LEVEL >= 2):
							#fp.write('\n ---- ' + str(clust_i) + ' ----> ' + str(clust_j) + ', ' + str(clust_j) + ' ----> ' + str(clust_k) + ', and ' \
								#+ str(clust_i) + ' -> ' + str(clust_k) + ' --- so removing ' + str(clust_i) + ' -> ' + str(clust_k))

	# close the output text file
	if (DEBUG_LEVEL >= 2):
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
	if target_taxa_label not in Cluster_Info_Dict[target_clust_idx]._GetSpeciesList():
		Cluster_Info_Dict[target_clust_idx]._Append_taxa(target_taxa_label)
		# mention the cluster index in the taxa information
		taxa_key = COMPLETE_INPUT_TAXA_LIST.index(target_taxa_label)
		Taxa_Info_Dict[taxa_key]._Set_Clust_Idx_taxa_Part(target_clust_idx)  

