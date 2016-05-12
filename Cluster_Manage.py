#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#-----------------------------------------------
"""
this function checks whether any given pair of clusters has no relation at all
"""
def Check_Cluster_Not_Connected(Reachability_Graph_Mat, outfile):
	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')
		fp.write('\n\n\n **** Within function --- Check_Cluster_Not_Connected **** \n\n')
	
	no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
	for i in range(no_of_clusters - 1):
		for j in range(i+1, no_of_clusters):
			if (Reachability_Graph_Mat[i][j] == 0) and (Reachability_Graph_Mat[j][i] == 0):
				clust1 = CURRENT_CLUST_IDX_LIST[i]
				clust2 = CURRENT_CLUST_IDX_LIST[j]
				taxa_list1 = Cluster_Info_Dict[clust1]._GetSpeciesList()
				taxa_list2 = Cluster_Info_Dict[clust2]._GetSpeciesList()
				if (DEBUG_LEVEL >= 2):
					fp.write('\n\n\n Cluster pair ' + str(clust1) + ' and ' + str(clust2) + '  with reach mat entries '+ str(i) \
						+ ' and ' + str(j) + '  and species lists : ( ' + str(taxa_list1) + ')  and  (' + str(taxa_list2) + ') is not connected ------ ')
				
				"""
				obtain the frequency and support score measures for this pair of taxa clusters
				"""
				R1_freq, R2_freq, R4_freq, R1_score, R2_score, R4_score, couplet_count, \
					r1_reln_allowed, r2_reln_allowed = GetFreqScore_ClusterPair(taxa_list1, taxa_list2)
				fp.write('\n R1_freq: ' + str(R1_freq) + ' R2_freq: ' + str(R2_freq) + ' R4_freq: ' + str(R4_freq) + \
					' couplet_count: ' + str(couplet_count) + ' r1_reln_allowed: ' + str(r1_reln_allowed) + ' r2_reln_allowed: ' + str(r2_reln_allowed)) 
				
				"""
				if R1 relation or R2 relation is strictly consensus then apply the relation between this cluster pair
				"""
				if (r1_reln_allowed == False):
					R1_freq = 0
				if (r2_reln_allowed == False):
					R2_freq = 0
				
				if (R1_freq > R2_freq) and (R1_freq > R4_freq):
					fp.write('\n Currently creating R1 relation from the cluster ' + str(clust1) + ' to the cluster ' + str(clust2))
					Reachability_Graph_Mat = Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(clust1), \
						CURRENT_CLUST_IDX_LIST.index(clust2), RELATION_R1, clust1, clust2)
				
				if (R2_freq > R1_freq) and (R2_freq > R4_freq):
					fp.write('\n Currently creating R1 relation from the cluster ' + str(clust2) + ' to the cluster ' + str(clust1))
					Reachability_Graph_Mat = Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(clust2), \
						CURRENT_CLUST_IDX_LIST.index(clust1), RELATION_R1, clust2, clust1)
					
				if (R4_freq > R1_freq) and (R4_freq > R2_freq):
					fp.write('\n Currently creating R4 relation from the cluster ' + str(clust1) + ' to the cluster ' + str(clust2))
					Reachability_Graph_Mat = Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(clust1), \
						CURRENT_CLUST_IDX_LIST.index(clust2), RELATION_R4, clust1, clust2)

	if (DEBUG_LEVEL >= 2):
		fp.write('\n\n\n **** Out from the function --- Check_Cluster_Not_Connected **** \n\n')
		fp.close()
	
	return Reachability_Graph_Mat

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
"""
this function checks and connects two cluster pairs with pseudo R1 / R2 connectivity
"""
def CheckPseudoR1RelnClust(clust1, clust2, outfile):
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
				if (supp_tree >= 1) and (fr4 > 0) and (fr4 == (pr1 + pr3)):
					r1_reln = True
				if (supp_tree >= 1) and (fr4 > 0) and (fr4 == (pr2 + pr3)):
					r2_reln = True
				if (r1_reln == True) and (r2_reln == True):
					break
	
	if (DEBUG_LEVEL > 2):
		fp = open(outfile, 'a')
		fp.write('\n In function CheckPseudoR1RelnClust -- clusters : ' + str(clust1) + '  and  ' + str(clust2) + \
			' r1_reln: ' + str(r1_reln) + ' r2_reln: ' + str(r2_reln))
		fp.close()
		
	return r1_reln, r2_reln

#-----------------------------------------------------
"""
this function checks the cluster pairs having R4 relation between them
for such clusters, we check if there exists possible R1 / R2 relations between them
and update the cluster information accordingly
"""
def Check_PossibleR1R2Reln(Reachability_Graph_Mat, outfile):

	# total number of clusters
	no_of_clusters = len(CURRENT_CLUST_IDX_LIST)

	for i in range(no_of_clusters - 1):
		clust1 = CURRENT_CLUST_IDX_LIST[i]
		for j in range(i+1, no_of_clusters):
			clust2 = CURRENT_CLUST_IDX_LIST[j]
			if (Reachability_Graph_Mat[i][j] == 2) and (Reachability_Graph_Mat[j][i] == 2):
				if (DEBUG_LEVEL > 2):
					fp = open(outfile, 'a')
					fp.write('\n\n\n Checking possible R4 (R1/R2) relations between the clusters ' + str(clust1) + ' and ' + str(clust2))
					fp.close()
				"""
				this cluster pair is related by R4 relation
				so check whether there exists possible R1 / R2 relations between this pair of clusters
				"""
				r1_reln, r2_reln = CheckPseudoR1RelnClust(clust1, clust2, outfile)
				"""
				we add pseudo R1 / R2 edges, depending on whether it satisfies the given relations
				"""
				if (r1_reln == True):
					Connect_PossibleR1Reln_ClusterPair(clust1, clust2)
					if (DEBUG_LEVEL > 2):
						fp = open(outfile, 'a')
						fp.write('\n Connected possible R4 (R1) relations from the cluster ' + str(clust1) + ' to the cluster ' + str(clust2))
						fp.close()
				if (r2_reln == True):
					Connect_PossibleR1Reln_ClusterPair(clust2, clust1)
					if (DEBUG_LEVEL > 2):
						fp = open(outfile, 'a')
						fp.write('\n Connected possible R4 (R1) relations from the cluster ' + str(clust2) + ' to the cluster ' + str(clust1))
						fp.close()
	
	return

#-----------------------------------------------------
""" 
this function adds an edge between a pair of clusters (of taxa) 
it also updates the entries of reachability matrix 
"""
def Connect_ClusterPair(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, \
	reln_type, nodeA_clust_idx, nodeB_clust_idx):
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
		
	return Reachability_Graph_Mat

#-----------------------------------------------------        
"""
Associated function of Solve_MPP_PossibleR1R2 / RemovePossibleR2List_DiffParent
resolves the inp_score_list if its length is > 1
by discarding possible clusters

@parameters:
	curr_clust_idx: the cluster whose possible R2 relation clusters are checked
	inp_score_list: contains the scores and cluster indices
	ascending_order: if True, the scores in the inp_score_list are sorted in the ascending order
										if False, scores are sorted in the descending order
"""
def ResolveScoreList(curr_clust_idx, inp_score_list, outfile, ascending_order=True):
	
	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')
	
	if (ascending_order == True):
		"""
		here the first element of the score list has the minimum value
		we note this value
		we compare all the other score list elements according to this minimum value
		and delete the clusters having values greater than this minimum value
		"""
		min_val = inp_score_list[0][1]
		for i in range(1, len(inp_score_list)):
			if (inp_score_list[i][1] > min_val):
				"""
				delete this entry from the possible R2 list
				"""
				target_delete_clust_idx = inp_score_list[i][0]
				Cluster_Info_Dict[curr_clust_idx]._RemovePossibleR2(target_delete_clust_idx)
				Cluster_Info_Dict[target_delete_clust_idx]._RemovePossibleR1(curr_clust_idx)
				if (DEBUG_LEVEL >= 2):
					fp.write('\n Removed possible R2 edge from the cluster ' + str(curr_clust_idx) \
						+ '  to the cluster: ' + str(target_delete_clust_idx))
	else:
		"""
		here the first element of the score list has the maximum value
		we note this value
		we compare all the other score list elements according to this maximum value
		and delete the clusters having values less than this maximum value
		"""
		max_val = inp_score_list[0][1]
		for i in range(1, len(inp_score_list)):
			if (inp_score_list[i][1] < max_val):
				"""
				delete this entry from the possible R2 list
				"""
				target_delete_clust_idx = inp_score_list[i][0]
				Cluster_Info_Dict[curr_clust_idx]._RemovePossibleR2(target_delete_clust_idx)
				Cluster_Info_Dict[target_delete_clust_idx]._RemovePossibleR1(curr_clust_idx)
				if (DEBUG_LEVEL >= 2):
					fp.write('\n Removed possible R2 edge from the cluster ' + str(curr_clust_idx) \
						+ '  to the cluster: ' + str(target_delete_clust_idx))
	
	if (DEBUG_LEVEL >= 2):
		fp.close()
	return

#-----------------------------------------------------        
"""
this function 
"""
def ComputeClusterParentDistance(src_clust, dest_clust):
	dist = 0
	while (1):
		if (src_clust == dest_clust):
			break
		dist_clust_parent_list = Cluster_Info_Dict[dest_clust]._GetClustRelnList(RELATION_R2)
		if (len(dist_clust_parent_list) == 0):
			return -1
		else:
			dist = dist + 1
			dest_clust = dist_clust_parent_list[0]
	
	return dist

#-----------------------------------------------------        
"""
Resolve the possible R2 lists for individual clusters
Let A--->B and C---->B
we have to choose between A and C
Let X->B,   X->A,   Y->C,  and  X->Y
In such a case, we choose C--->B
"""
def RemovePossibleR2List_DiffParent(Reachability_Graph_Mat, outfile):
	
	for cx in Cluster_Info_Dict:
		if (Cluster_Info_Dict[cx]._Get_Indegree() > 0):
			if (DEBUG_LEVEL >= 2):
				fp = open(outfile, 'a')
				fp.write('\n\n\n ***** In function - RemovePossibleR2List_DiffParent - Examining cluster -- ' + str(cx))
				fp.close()
			cx_parent = (Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2))[0]
			cx_possible_R2_list = Cluster_Info_Dict[cx]._GetPossibleR2List()
			cx_possible_R2_list_copy = list(cx_possible_R2_list)
			if (DEBUG_LEVEL >= 2):
				fp = open(outfile, 'a')
				fp.write('\n cx_possible_R2_list_copy: ' + str(cx_possible_R2_list_copy))
				fp.close()
			if (len(cx_possible_R2_list_copy) > 0):
				"""
				define an empty list which will contain the difference of levels between cx_parent 
				and the parent cluster of individual cluster x within cx_possible_R2_list_copy
				"""
				score_cx = []
				for x in cx_possible_R2_list_copy:
					if (Cluster_Info_Dict[x]._Get_Indegree() > 0):
						x_parent = (Cluster_Info_Dict[x]._GetClustRelnList(RELATION_R2))[0]
						if (DEBUG_LEVEL >= 2):
							fp = open(outfile, 'a')
							fp.write('\n cx_parent: ' + str(cx_parent) + ' x: ' + str(x) + ' x_parent: ' + str(x_parent))
							fp.close()
						if (cx_parent != x_parent) and \
							(Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(cx_parent)][CURRENT_CLUST_IDX_LIST.index(x_parent)] != 1):
							Cluster_Info_Dict[cx]._RemovePossibleR2(x)
							Cluster_Info_Dict[x]._RemovePossibleR1(cx)
							if (DEBUG_LEVEL >= 2):
								fp = open(outfile, 'a')
								fp.write('\n Removed the cluster : ' + str(x) + '  from the possible R2 list of ' + str(cx))
								fp.close()
						else:
							"""
							either cx_parent = x_parent
							or, cx_parent -> x_parent
							"""
							dist = ComputeClusterParentDistance(cx_parent, x_parent)
							if (DEBUG_LEVEL >= 2):
								fp = open(outfile, 'a')
								fp.write('\n Computed distance between cx_parent: ' + str(cx_parent) \
									+ ' and x_parent: ' + str(x_parent) + '  is: ' + str(dist))
								fp.close()
							"""
							insert the distance and the cluster x information in the scoring list
							"""
							subl = [x, dist]
							score_cx.append(subl)
							
				"""
				sort the scoring list containing only the possible R2 relations in descending order, 
				with respect to the distance between x_parent and cx_parent
				"""
				if (len(score_cx) > 1):
					score_cx.sort(key=lambda x: x[1], reverse=True)
					if (DEBUG_LEVEL >= 2):
						fp = open(outfile, 'a')
						fp.write('\n --- Before ResolveScoreList -- Sorted scoring list (cluster parent distance) : ' + str(score_cx))
						fp.close()
					ResolveScoreList(cx, score_cx, outfile, False)
				
				#"""
				#if the highest distance between x_parent and cx_parent (where x --> cx is established) is > 0
				#then the following situation arises
				
				#cx_parent->cx
				#x_parent->x
				#x--> cx (and cx is not ---> x)
				#x_parent is a descendant of cx_parent
				
				#in such a case, we remove cx_parent from the parent list of cx
				#instead, we rely on x to traverse and find cx
				#or we establish x_parent as a parent of cx
				#"""
				#if (len(score_cx) > 0):
					#if (score_cx[0][1] > 0):
						#x = score_cx[0][0]
						#x_parent = (Cluster_Info_Dict[x]._GetClustRelnList(RELATION_R2))[0]
						#if (DEBUG_LEVEL >= 2):
							#fp = open(outfile, 'a')
							#fp.write('\n --- Max distance : ' + str(score_cx[0][1]) + '  Possible R2 cluster: ' + str(x))
							#fp.close()
						#"""
						#remove cx_parent from the parent list of cx
						#"""
						#Cluster_Info_Dict[cx]._RemoveRelnInstance(RELATION_R2, cx_parent)
						#if (DEBUG_LEVEL >= 2):
							#fp = open(outfile, 'a')
							#fp.write(' ---  removed the original parent: ' + str(cx_parent))
							#fp.close()
						#"""
						#if x -->cx and also x<---cx, then copy x_parent as a parent of cx
						#"""
						#if x in Cluster_Info_Dict[cx]._GetPossibleR1List():
							#Cluster_Info_Dict[cx]._AddRelnInstance(RELATION_R2, x_parent)
							#if (DEBUG_LEVEL >= 2):
								#fp = open(outfile, 'a')
								#fp.write(' ---  Possible R1 cluster: ' + str(x) + ' addition of parent node : ' + str(x_parent))
								#fp.close()
	
	return
	

#-----------------------------------------------------        
"""
solves multiple parent problem from the possible R1 / R2 lists 
the selection is carried out using one of the following measure:
	internode count based scoring mechanism
This function is applicable for directed in edge (Relation R2) only
"""
def Solve_MPP_PossibleR1R2(Output_Text_File, Cluster_Dict_Backup):
	#--------------------------------------------------------------
	"""
	now we check the clusters cx having following properties:
	1) No cluster cy exists such that cy -> cx holds
	2) There exists at least two clusters c1 and c2 such that c1--->cx and c2--->cx
	in such a case, we select among the multiple dashed edges (possible R2 lists), the suitable one or more possible R2 candidate
	"""
	for cx in Cluster_Info_Dict:
		
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n\n\n ***** In function - Solve_MPP_PossibleR1R2 - Examining cluster -- ' + str(cx))
			fp.close()
		
		"""
		explore all clusters belonging to the possible R2 list of the cluster cx
		that is, for each cluster cz in the possible R2 list, 
		cz--->cx relation is established with respect to input gene trees
		"""
		cx_possible_R2_list = Cluster_Info_Dict[cx]._GetPossibleR2List()
		cx_possible_R1_list = Cluster_Info_Dict[cx]._GetPossibleR1List()
		cx_R2_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2)

		#-----------------------------------------
		"""
		explore the clusters belonging to the R2 list of the current cluster
		"""
		if (len(cx_R2_list) > 0):
			score_list_cx_R2 = []
			# open the text file
			fp = open(Output_Text_File, 'a')
			for cz in cx_R2_list:
				cz_score3 = FindAvgInternodeCount(Cluster_Info_Dict[cz]._GetSpeciesList(), Cluster_Info_Dict[cx]._GetSpeciesList(), 2, 1)
				temp_subl = [cz, cz_score3]
				score_list_cx_R2.append(temp_subl)
				if (DEBUG_LEVEL >= 2):
					fp.write('\n --- element (true R2 reln): ' + str(cz) + ' Internode count score: ' + str(cz_score3))
			# close the text file
			fp.close()
		#-----------------------------------------
		
		"""
		first we check only those clusters (cx) such that the length of cx_possible_R2_list is > 0, 
		in such a case, we have at least one cluster cz such that cz--->cx holds
		"""
		if (len(cx_possible_R2_list) > 0):	#1):
			
			# open the text file
			fp = open(Output_Text_File, 'a')
			
			"""
			initialize two scoring lists corresponding to the cluster cx
			one list contains the clusters connected to cx via both possible R1 and possible R2 list
			another list contains the clusters connected to cx via only possible R2 list
			"""
			score_list_cx_both = []
			score_list_cx_single = []
			"""
			explore all candidate clusters cz belonging to the possible R2 list of the cluster cx
			"""
			for cz in cx_possible_R2_list:
				cz_score3 = FindAvgInternodeCount(Cluster_Info_Dict[cz]._GetSpeciesList(), Cluster_Info_Dict[cx]._GetSpeciesList(), 2, 1)
				if (DEBUG_LEVEL > 2):
					fp.write('\n --- element (possible R4 (R2) reln): ' + str(cz) + ' Internode count score: ' + str(cz_score3))
				temp_subl = [cz, cz_score3]
				if (cz in cx_possible_R1_list):
					score_list_cx_both.append(temp_subl)
				else:
					score_list_cx_single.append(temp_subl)
		
			"""
			sort the scoring list containing both possible R1 and R2 relations in ascending order, 
			with respect to the internode count
			"""
			if (len(score_list_cx_both) > 0):
				score_list_cx_both.sort(key=lambda x: x[1])
				if (DEBUG_LEVEL >= 2):
					fp.write('\n --- after sorting the scoring list (both possible R1 and R2 relations) corresponding to the cluster : ' + str(cx))
					for i in range(len(score_list_cx_both)):
						fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(score_list_cx_both[i][0]) + \
							' Internode count score: ' + str(score_list_cx_both[i][1]))
			
			"""
			sort the scoring list containing only the possible R2 relations in ascending order, 
			with respect to the internode count
			"""
			if (len(score_list_cx_single) > 0):
				score_list_cx_single.sort(key=lambda x: x[1])
				if (DEBUG_LEVEL >= 2):
					fp.write('\n --- after sorting the scoring list (only possible R2 relations) corresponding to the cluster : ' + str(cx))
					for i in range(len(score_list_cx_single)):
						fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(score_list_cx_single[i][0]) + \
							' Internode count score: ' + str(score_list_cx_single[i][1]))
			
			# close the text file
			fp.close()
		
			if (len(score_list_cx_both) > 1):
				ResolveScoreList(cx, score_list_cx_both, Output_Text_File)
			
			if (len(score_list_cx_single) > 1):
				ResolveScoreList(cx, score_list_cx_single, Output_Text_File)

			#-------------------------------
			# add - sourya
			
			"""
			if X-->Y and Z->Y and X, Y are not related by any means
			then we choose one of the parents (or pseudo parent) for Y
			"""
			
			if (len(cx_R2_list) > 0) and (len(cx_possible_R2_list) > 0):
				flag = False
				cx_R2_clust = score_list_cx_R2[0][0]
				if (len(score_list_cx_both) > 0):
					cx_R2_both_clust = score_list_cx_both[0][0]
					
					if ((cx_R2_both_clust in Cluster_Dict_Backup[cx_R2_clust]._GetClustRelnList(RELATION_R1)) or \
						((cx_R2_both_clust in Cluster_Dict_Backup[cx_R2_clust]._GetClustRelnList(RELATION_R2)) or \
							(cx_R2_both_clust in Cluster_Dict_Backup[cx_R2_clust]._GetPossibleR1List())) or \
								(cx_R2_both_clust in Cluster_Dict_Backup[cx_R2_clust]._GetPossibleR2List())):
						flag = True 

				if (len(score_list_cx_single) > 0):
					cx_R2_single_clust = score_list_cx_single[0][0]

					if ((cx_R2_single_clust in Cluster_Dict_Backup[cx_R2_clust]._GetClustRelnList(RELATION_R1)) or \
						((cx_R2_single_clust in Cluster_Dict_Backup[cx_R2_clust]._GetClustRelnList(RELATION_R2)) or \
							(cx_R2_single_clust in Cluster_Dict_Backup[cx_R2_clust]._GetPossibleR1List())) or \
								(cx_R2_single_clust in Cluster_Dict_Backup[cx_R2_clust]._GetPossibleR2List())):
						flag = True 
				
				if (flag == False):
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n\n ** The cluster has one parent cluster and at least one possible R2 cluster ---- without any connection between them')
						fp.close()
						
					if (len(score_list_cx_both) > 0):
						if (score_list_cx_R2[0][1] > score_list_cx_both[0][1]):
							# remove cx_R2_clust from the parent list of cx
							Remove_ClusterPairConn(cx, cx_R2_clust, RELATION_R2)
							if (DEBUG_LEVEL >= 2):
								fp = open(Output_Text_File, 'a')
								fp.write('\n\n ** Removed parent cluster ' + str(cx_R2_clust) + '  from the current cluster')
								fp.close()
			
					if (len(score_list_cx_single) > 0):
						if (score_list_cx_R2[0][1] > score_list_cx_single[0][1]):
							# remove cx_R2_clust from the parent list of cx
							Remove_ClusterPairConn(cx, cx_R2_clust, RELATION_R2)
							if (DEBUG_LEVEL >= 2):
								fp = open(Output_Text_File, 'a')
								fp.write('\n\n ** Removed parent cluster ' + str(cx_R2_clust) + '  from the current cluster')
								fp.close()
			
			
			# end add - sourya
			#-------------------------------
			
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
										Reachability_Graph_Mat = Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
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
											Reachability_Graph_Mat = Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
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
										Reachability_Graph_Mat = Connect_PossibleR1Reln_ClusterPair(Reachability_Graph_Mat, x, cx)
										flag = True
	
	
				# close the text file
				fp.close()

	return Reachability_Graph_Mat

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

	# open the text file
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')

	for cx in Cluster_Info_Dict:

		"""
		list of clusters belonging to the R2 list of cx
		"""
		cx_R2_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2)
		if (len(cx_R2_list) > 1):
			
			if (DEBUG_LEVEL >= 2):
				fp.write('\n\n ***** In function SelectUniqueParent_Directed ---- Examining cluster -- ' + str(cx))
			
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
			sort the scoring list in ascending order
			for internode count based measure, select the cluster having the lowest average internode count
			"""
			score_list_cx.sort(key=lambda x: x[1])	#, reverse=True)
			if (DEBUG_LEVEL >= 2):
				fp.write('\n --- after sorting the scoring list corresponding to the cluster : ' + str(cx))
				for i in range(len(score_list_cx)):
					fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(score_list_cx[i][0]) + ' score: ' + str(score_list_cx[i][1]))

			"""
			for the internode count based measure, remove all except the first element from the 
			scoring list of the current cluster cx
			since the first element contains the highest internode count measure
			"""
			min_internode_measure = score_list_cx[0][1]
			for i in range(1, len(score_list_cx)):
				if (score_list_cx[i][1] > min_internode_measure):
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
	
	## --------------- old code - sourya
	
	#for i in Cluster_Info_Dict:
		#if (Cluster_Info_Dict[i]._GetExploredStatus() == 0):
			## we check the clusters which have not been explored yet
			#if (valid_node_found == 0):
				#min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
				#min_indeg_node_idx = i
				#valid_node_found = 1
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Minimum indegree cluster so far: ' + str(i))
			#elif (valid_node_found == 1) and (Cluster_Info_Dict[i]._Get_Indegree() < min_indeg):
				#min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
				#min_indeg_node_idx = i
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Minimum indegree cluster: ' + str(i))
			#elif (valid_node_found == 1) and (Cluster_Info_Dict[i]._Get_Indegree() == min_indeg):
				#if (Cluster_Info_Dict[i]._Get_Outdegree() > Cluster_Info_Dict[min_indeg_node_idx]._Get_Outdegree()):
					#min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
					#min_indeg_node_idx = i
					#if (DEBUG_LEVEL >= 2):
						#fp.write('\n Minimum indegree cluster: ' + str(i))
				##elif (valid_node_found == 1) and \
					##(len(Cluster_Info_Dict[i]._GetDistinctPossibleR2List()) == len(Cluster_Info_Dict[min_indeg_node_idx]._GetDistinctPossibleR2List())) \
					##and ((Cluster_Info_Dict[i]._Get_Outdegree() + len(Cluster_Info_Dict[i]._GetPossibleR1List())) > \
						##(Cluster_Info_Dict[min_indeg_node_idx]._Get_Outdegree() + len(Cluster_Info_Dict[min_indeg_node_idx]._GetPossibleR1List()))):    
					##min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
					##min_indeg_node_idx = i
					##if (DEBUG_LEVEL >= 2):
						##fp.write('\n Minimum indegree cluster: ' + str(i))


	#-------------------- new code - sourya

	for i in clust_dict:
		if (clust_dict[i]._GetExploredStatus() == 0):
			# we check the clusters which have not been explored yet
			if (valid_node_found == 0):
				min_indeg = Get_Non_Explored_Indegree(clust_dict, i)
				min_DistPossR2 = Get_Non_Explored_DistinctPossibleR2Count(clust_dict, i)
				min_outdeg = Get_Non_Explored_Outdegree(clust_dict, i)
				min_DistPossR1 = Get_Non_Explored_PossibleR1Count(clust_dict, i)
				min_indeg_node_idx = i
				valid_node_found = 1
				if (DEBUG_LEVEL >= 2):
					fp.write('\n Minimum non explored indegree cluster so far: ' + str(i))
			else:
				"""
				here valid_node_found = 1
				so at least one minimum value is already obtained 
				"""
				curr_indeg = Get_Non_Explored_Indegree(clust_dict, i)
				curr_DistPossR2 = Get_Non_Explored_DistinctPossibleR2Count(clust_dict, i)
				curr_outdeg = Get_Non_Explored_Outdegree(clust_dict, i)
				curr_DistPossR1 = Get_Non_Explored_PossibleR1Count(clust_dict, i)
				if (curr_indeg < min_indeg) or ((curr_indeg == min_indeg) and (curr_DistPossR2 < min_DistPossR2)) \
					or ((curr_indeg == min_indeg) and (curr_DistPossR2 == min_DistPossR2) and (curr_outdeg > min_outdeg)):
					min_indeg = curr_indeg
					min_DistPossR2 = curr_DistPossR2
					min_outdeg = curr_outdeg
					min_DistPossR1 = curr_DistPossR1
					min_indeg_node_idx = i
					if (DEBUG_LEVEL >= 2):
						fp.write('\n Minimum non explored indegree cluster: ' + str(i))

	if (DEBUG_LEVEL >= 2):
		fp.write('\n\n Out of the function --- Extract_Node_Min_Indeg \n\n ')
		fp.close()
	# ---------------------- end new code - sourya
	
	return min_indeg_node_idx

#-----------------------------------------------------  
def Get_Non_Explored_Indegree(clust_dict, i):
	c = 0
	for x in clust_dict[i]._GetClustRelnList(RELATION_R2):
		if (clust_dict[x]._GetExploredStatus() == 0):
			c = c + 1
	return c

def Get_Non_Explored_Outdegree(clust_dict, i):
	c = 0
	for x in clust_dict[i]._GetClustRelnList(RELATION_R1):
		if (clust_dict[x]._GetExploredStatus() == 0):
			c = c + 1
	return c

def Get_Non_Explored_DistinctPossibleR2Count(clust_dict, i):
	c = 0
	for x in clust_dict[i]._GetDistinctPossibleR2List():
		if (clust_dict[x]._GetExploredStatus() == 0):
			c = c + 1
	return c

def Get_Non_Explored_PossibleR1Count(clust_dict, i):
	c = 0
	for x in clust_dict[i]._GetPossibleR1List():
		if (clust_dict[x]._GetExploredStatus() == 0):
			c = c + 1
	return c
#-----------------------------------------------------  
"""
this function separately contains the code for transitive reduction of the dashed edges
"""
def Transitive_Reduction_Dashed(Clust_Possible_R1_Mat, Outfile):
	no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
	
	# open the output text file
	if (DEBUG_LEVEL >= 2):
		fp = open(Outfile, 'a')

	"""
	transitive reduction for the following case: 
	A--->B, B---->C, A---->C ---------- remove A---->C 
	Note: we should not consider the case when B--->C and C--->B simultaneously
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
					# B--->C case
					if (Clust_Possible_R1_Mat[j][k] == 1) and (Clust_Possible_R1_Mat[k][j] == 0):
						# A--->C case
						if (Clust_Possible_R1_Mat[i][k] == 1) and (Clust_Possible_R1_Mat[k][i] == 0):
							# comment - sourya - check
							#Clust_Possible_R1_Mat[i][k] = 0
							Cluster_Info_Dict[clust_i]._RemovePossibleR1(clust_k)
							Cluster_Info_Dict[clust_k]._RemovePossibleR2(clust_i)
							if (DEBUG_LEVEL >= 2):
								fp.write('\n ---- ' + str(clust_i) + ' ----> ' + str(clust_j) + ', ' + str(clust_j) + ' ---->' + str(clust_k) + ', and ' \
									+ str(clust_i) + ' -----> ' + str(clust_k) + ' --- so removing ' + str(clust_i) + ' ----> ' + str(clust_k))

	#"""
	#transitive reduction for the following case: 
	#A--->B / A<--->B, B---->C, A---->C or A<---->C ---------- remove A---->C / A<---->C
	#Note: we should not consider the case when B--->C and C--->B simultaneously
	#"""
	#for j in range(no_of_clusters):
		#clust_j = CURRENT_CLUST_IDX_LIST[j]
		#for i in range(no_of_clusters):
			#if (i == j):
				#continue
			#clust_i = CURRENT_CLUST_IDX_LIST[i]
			## A--->B / A<--->B case
			#if (Clust_Possible_R1_Mat[i][j] == 1):	# and (Clust_Possible_R1_Mat[j][i] == 0):
				#for k in range(no_of_clusters):
					#if (i == k) or (j == k):
						#continue
					#clust_k = CURRENT_CLUST_IDX_LIST[k]
					## B--->C case
					#if (Clust_Possible_R1_Mat[j][k] == 1) and (Clust_Possible_R1_Mat[k][j] == 0):
						## A--->C / A<--->C case
						#if (Clust_Possible_R1_Mat[i][k] == 1):	# and (Clust_Possible_R1_Mat[k][i] == 0):
							## comment - sourya - check
							##Clust_Possible_R1_Mat[i][k] = 0
							#Cluster_Info_Dict[clust_i]._RemovePossibleR1(clust_k)
							#Cluster_Info_Dict[clust_k]._RemovePossibleR2(clust_i)
							#if (DEBUG_LEVEL >= 2):
								#fp.write('\n ---- ' + str(clust_i) + ' ----> / <----> ' + str(clust_j) + ', ' + str(clust_j) + ' ---->' + str(clust_k) + ', and ' \
									#+ str(clust_i) + ' -----> / <----> ' + str(clust_k) + ' --- so removing ' + str(clust_i) + ' ----> / <----> ' + str(clust_k))


	#"""
	#transitive reduction for the following case: 
	#A--->B, B---->C / B<---->C , A---->C or A<---->C ---------- remove A---->C / A<---->C
	#Note: we should not consider the case when A--->B and B--->A simultaneously
	#"""
	#for j in range(no_of_clusters):
		#clust_j = CURRENT_CLUST_IDX_LIST[j]
		#for i in range(no_of_clusters):
			#if (i == j):
				#continue
			#clust_i = CURRENT_CLUST_IDX_LIST[i]
			## A--->B case
			#if (Clust_Possible_R1_Mat[i][j] == 1) and (Clust_Possible_R1_Mat[j][i] == 0):
				#for k in range(no_of_clusters):
					#if (i == k) or (j == k):
						#continue
					#clust_k = CURRENT_CLUST_IDX_LIST[k]
					## B--->C / B<---->C case
					#if (Clust_Possible_R1_Mat[j][k] == 1):	# and (Clust_Possible_R1_Mat[k][j] == 0):
						## A--->C / A<--->C case
						#if (Clust_Possible_R1_Mat[i][k] == 1):	# and (Clust_Possible_R1_Mat[k][i] == 0):
							## comment - sourya - check
							##Clust_Possible_R1_Mat[i][k] = 0
							#Cluster_Info_Dict[clust_i]._RemovePossibleR1(clust_k)
							#Cluster_Info_Dict[clust_k]._RemovePossibleR2(clust_i)
							#if (DEBUG_LEVEL >= 2):
								#fp.write('\n ---- ' + str(clust_i) + ' ----> ' + str(clust_j) + ', ' + str(clust_j) + ' ----> / <---->' + str(clust_k) + ', and ' \
									#+ str(clust_i) + ' -----> / <----> ' + str(clust_k) + ' --- so removing ' + str(clust_i) + ' ----> / <----> ' + str(clust_k))


	# close the output text file
	if (DEBUG_LEVEL >= 2):
		fp.close()

	return

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
							#if (Reachability_Graph_Mat[i][k] == 1):
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
						#if (Reachability_Graph_Mat[j][k] == 1):
						if clust_k in Cluster_Info_Dict[clust_j]._GetClustRelnList(RELATION_R1):
							#if (Reachability_Graph_Mat[i][k] == 1):
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
def CompressDirectedGraph(Reachability_Graph_Mat, Outfile):
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
	target_clust_spec_list = Cluster_Info_Dict[target_clust_idx]._GetSpeciesList()
	if target_taxa_label not in target_clust_spec_list:
		Cluster_Info_Dict[target_clust_idx]._Append_taxa(target_taxa_label)
		# mention the cluster index in the taxa information
		taxa_key = COMPLETE_INPUT_TAXA_LIST.index(target_taxa_label)
		Taxa_Info_Dict[taxa_key]._Set_Clust_Idx_taxa_Part(target_clust_idx)  

