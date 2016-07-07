#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *
import Conflict_Detect
from Conflict_Detect import *

##-----------------------------------------------
#"""
#this function checks whether any given pair of clusters has no relation at all
#"""
#def Check_Cluster_Not_Connected(Reachability_Graph_Mat, outfile):
	#if (DEBUG_LEVEL >= 2):
		#fp = open(outfile, 'a')
		#fp.write('\n\n\n **** Within function --- Check_Cluster_Not_Connected **** \n\n')
		#fp.close()
	
	#no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
	#for i in range(no_of_clusters - 1):
		#for j in range(i+1, no_of_clusters):
			#if (Reachability_Graph_Mat[i][j] == 0) and (Reachability_Graph_Mat[j][i] == 0):
				#clust1 = CURRENT_CLUST_IDX_LIST[i]
				#clust2 = CURRENT_CLUST_IDX_LIST[j]
				#if (DEBUG_LEVEL >= 2):
					#fp = open(outfile, 'a')
					#fp.write('\n\n\n Cluster pair ' + str(clust1) + ' and ' + str(clust2) + ' is not connected ------ ')
				
				#"""
				#obtain the parameters from this cluster pair
				#"""
				#clust_pair_key = (clust1, clust2)
				#R1_freq = Cluster_Pair_Info_Dict[clust_pair_key]._GetFreq(RELATION_R1)
				#R2_freq = Cluster_Pair_Info_Dict[clust_pair_key]._GetFreq(RELATION_R2)
				#R3_freq = Cluster_Pair_Info_Dict[clust_pair_key]._GetFreq(RELATION_R3)
				#R4_freq = Cluster_Pair_Info_Dict[clust_pair_key]._GetFreq(RELATION_R4)
				##R1_score = Cluster_Pair_Info_Dict[clust_pair_key]._GetSupportScore(RELATION_R1)
				##R2_score = Cluster_Pair_Info_Dict[clust_pair_key]._GetSupportScore(RELATION_R2)
				##R3_score = Cluster_Pair_Info_Dict[clust_pair_key]._GetSupportScore(RELATION_R3)
				##R4_score = Cluster_Pair_Info_Dict[clust_pair_key]._GetSupportScore(RELATION_R4)
				#couplet_count = Cluster_Pair_Info_Dict[clust_pair_key]._GetCoupletCount()
				##r1_reln_allowed = Cluster_Pair_Info_Dict[clust_pair_key]._GetRelnAllowed(RELATION_R1)
				##r2_reln_allowed = Cluster_Pair_Info_Dict[clust_pair_key]._GetRelnAllowed(RELATION_R2)
				##pseudo_r4_r1_count = Cluster_Pair_Info_Dict[clust_pair_key]._GetPseudoR4RelnCount(0) 
				##pseudo_r4_r2_count = Cluster_Pair_Info_Dict[clust_pair_key]._GetPseudoR4RelnCount(1)  
				##pseudo_r4_r3_count = Cluster_Pair_Info_Dict[clust_pair_key]._GetPseudoR4RelnCount(2)
				
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n R1_freq: ' + str(R1_freq) + ' R2_freq: ' + str(R2_freq) + ' R4_freq: ' + str(R4_freq) + \
						#' couplet_count: ' + str(couplet_count))	# + ' r1_reln_allowed: ' + str(r1_reln_allowed) + ' r2_reln_allowed: ' + str(r2_reln_allowed)) 
				
				#RELN_LIST = [RELATION_R4, RELATION_R1, RELATION_R2]
				
				#"""
				#if R1 relation or R2 relation is strictly consensus then apply the relation between this cluster pair
				#"""
				#if 1:	#(r1_reln_allowed == False):
					#R1_freq = 0
					#RELN_LIST.remove(RELATION_R1)
				#if 1:	#(r2_reln_allowed == False):
					#R2_freq = 0
					#RELN_LIST.remove(RELATION_R2)
				
				## slightly modified - sourya
				
				#if (R1_freq > R2_freq) and (R1_freq > R4_freq):
					#RELN_LIST.remove(RELATION_R1)
					#RELN_LIST.insert(0, RELATION_R1)
				
				#if (R2_freq > R1_freq) and (R2_freq > R4_freq):
					#RELN_LIST.remove(RELATION_R2)
					#RELN_LIST.insert(0, RELATION_R2)
				
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n RELN_LIST: ' + str(RELN_LIST))
					#fp.close()
				
				#for reln_type in RELN_LIST:	
					#if (Possible_Conflict_Curr_Reln(clust1, clust2, Reachability_Graph_Mat, reln_type, outfile) == 0):
						#Reachability_Graph_Mat = Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(clust1), \
							#CURRENT_CLUST_IDX_LIST.index(clust2), reln_type, clust1, clust2)
						#if (DEBUG_LEVEL >= 2):
							#fp = open(outfile, 'a')
							#fp.write('\n Established the relation: ' + str(reln_type) \
								#+ ' from the cluster ' + str(clust1) + ' to the cluster ' + str(clust2))
							#fp.close()
						#break
			
	#if (DEBUG_LEVEL >= 2):
		#fp = open(outfile, 'a')
		#fp.write('\n\n\n **** Out from the function --- Check_Cluster_Not_Connected **** \n\n')
		#fp.close()
	
	#return Reachability_Graph_Mat

#------------------------------------------------
"""
this function removes reln_type connection from clust1 to clust2
"""
def Remove_ClusterPairConn(clust1, clust2, reln_type):
	Cluster_Info_Dict[clust1]._RemoveRelnInstance(reln_type, clust2)
	Cluster_Info_Dict[clust2]._RemoveRelnInstance(Complementary_Reln(reln_type), clust1)
	return

##----------------------------------------------------
#"""
#this function establishes a possible R1 connection from clust1 to clust2
#"""
#def Connect_PossibleR1Reln_ClusterPair(clust1, clust2):
	#Cluster_Info_Dict[clust1]._AddPossibleR1(clust2)
	#Cluster_Info_Dict[clust2]._AddPossibleR2(clust1)
	#return

##-----------------------------------------------------
#"""
#this function checks and connects two cluster pairs with pseudo R1 / R2 connectivity
#"""
#def CheckPseudoR1RelnClust(clust1, clust2, outfile):
	#r1_reln = False
	#r2_reln = False
	
	#for t1 in Cluster_Info_Dict[clust1]._GetSpeciesList():
		#t1_idx = COMPLETE_INPUT_TAXA_LIST.index(t1)
		#for t2 in Cluster_Info_Dict[clust2]._GetSpeciesList():
			#t2_idx = COMPLETE_INPUT_TAXA_LIST.index(t2)
			#"""
			#formation of the couplet key
			#"""
			#if (t1_idx < t2_idx):
				#target_key = (t1_idx, t2_idx)
				#complement_operation = False
			#else:
				#target_key = (t2_idx, t1_idx)
				#complement_operation = True
				
			#if target_key in TaxaPair_Reln_Dict:
				#supp_tree = TaxaPair_Reln_Dict[target_key]._GetNoSupportTrees()
				#fr4 = TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R4)
				#pr3 = TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(2)
				#if (complement_operation == False):
					#pr1 = TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(0)
					#pr2 = TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(1)
				#else:
					#pr1 = TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(1)
					#pr2 = TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(0)
					
				#"""
				#sourya - previously the number of support trees were restricted as > 1
				#now we have included the condition >= 1 as well
				#"""
				#if (supp_tree >= 1) and (fr4 > 0) and (fr4 == (pr1 + pr3)):
					#r1_reln = True
				#if (supp_tree >= 1) and (fr4 > 0) and (fr4 == (pr2 + pr3)):
					#r2_reln = True
				#if (r1_reln == True) and (r2_reln == True):
					#break
	
	#if (DEBUG_LEVEL > 2):
		#fp = open(outfile, 'a')
		#fp.write('\n In function CheckPseudoR1RelnClust -- clusters : ' + str(clust1) + '  and  ' + str(clust2) + \
			#' r1_reln: ' + str(r1_reln) + ' r2_reln: ' + str(r2_reln))
		#fp.close()
		
	#return r1_reln, r2_reln

##-----------------------------------------------------
#"""
#checks for the possibility of pseudo R1 / R2 relation according to the connectivity of the parent clusters
#"""
#def CheckPseudoR1Reln_ParentBased_Clust(clust1, clust2, outfile):
	
	#clust1_parent_list = Cluster_Info_Dict[clust1]._GetClustRelnList(RELATION_R2)
	#clust2_parent_list = Cluster_Info_Dict[clust2]._GetClustRelnList(RELATION_R2)
	
	#r1_reln = True
	#r2_reln = True
	
	#if (len(clust1_parent_list) > 0) and (len(clust2_parent_list) > 0):
		#clust1_parent = clust1_parent_list[0]
		#clust2_parent = clust2_parent_list[0]
		#if (clust1_parent != clust2_parent):
			#clust2_parent_parent_list = Cluster_Info_Dict[clust2_parent]._GetClustRelnList(RELATION_R2)
			#if (len(clust2_parent_parent_list) == 0):
				#r1_reln = False
			#elif (clust1_parent != clust2_parent_parent_list[0]):
				#r1_reln = False
			#clust1_parent_parent_list = Cluster_Info_Dict[clust1_parent]._GetClustRelnList(RELATION_R2)
			#if (len(clust1_parent_parent_list) == 0):
				#r2_reln = False
			#elif (clust2_parent != clust1_parent_parent_list[0]):
				#r2_reln = False
		
	#return r1_reln, r2_reln

##-----------------------------------------------------
#"""
#this function checks the cluster pairs having R4 relation between them
#for such clusters, we check if there exists possible R1 / R2 relations between them
#and update the cluster information accordingly
#"""
#def Check_PossibleR1R2Reln(Reachability_Graph_Mat, outfile):

	## total number of clusters
	#no_of_clusters = len(CURRENT_CLUST_IDX_LIST)

	#for i in range(no_of_clusters - 1):
		#clust1 = CURRENT_CLUST_IDX_LIST[i]
		#for j in range(i+1, no_of_clusters):
			#clust2 = CURRENT_CLUST_IDX_LIST[j]
			#if (Reachability_Graph_Mat[i][j] == 2) and (Reachability_Graph_Mat[j][i] == 2):

				##clust1_parent_list = Cluster_Info_Dict[clust1]._GetClustRelnList(RELATION_R2)
				##clust2_parent_list = Cluster_Info_Dict[clust2]._GetClustRelnList(RELATION_R2)
				##if (DEBUG_LEVEL > 2):
					##fp = open(outfile, 'a')
					##fp.write('\n\n\n Checking possible R4 (R1/R2) relations between the clusters ' + str(clust1) + ' and ' + str(clust2) + \
						##'  their parent lists: ' + str(clust1_parent_list) + '  and  ' + str(clust2_parent_list))
					##fp.close()

				### first check according to the cluster based connectivity
				##r1_reln_1, r2_reln_1 = CheckPseudoR1Reln_ParentBased_Clust(clust1, clust2, outfile)
				
				## if any of the relations are true
				#if 1:	#(r1_reln_1 == True) or (r2_reln_1 == True):
				
					#"""
					#this cluster pair is related by R4 relation
					#so check whether there exists possible R1 / R2 relations between this pair of clusters
					#"""
					#r1_reln, r2_reln = CheckPseudoR1RelnClust(clust1, clust2, outfile)
				
				#"""
				#we add pseudo R1 / R2 edges, depending on whether it satisfies the given relations
				#"""
				#if (r1_reln == True):	# and (r1_reln_1 == True):
					#Connect_PossibleR1Reln_ClusterPair(clust1, clust2)
					#if (DEBUG_LEVEL > 2):
						#fp = open(outfile, 'a')
						#fp.write('\n Connected possible R4 (R1) relations from the cluster ' + str(clust1) + ' to the cluster ' + str(clust2))
						#fp.close()
				#if (r2_reln == True):	# and (r2_reln_1 == True):
					#Connect_PossibleR1Reln_ClusterPair(clust2, clust1)
					#if (DEBUG_LEVEL > 2):
						#fp = open(outfile, 'a')
						#fp.write('\n Connected possible R4 (R1) relations from the cluster ' + str(clust2) + ' to the cluster ' + str(clust1))
						#fp.close()
	
	#return

##-----------------------------------------------------
#"""
#@parameters:
#clustA and clustB: two input clusters such that clustA->clustB connection is sought
#the function checks whether clustA->clustB connection is basically a R4 relation with pseudo R1 connection
#in such a case, clustB is added in the pseudo R1 list of clustA 
#"""
#def Check_Possible_R1_Reln(clustA, clustB):
	#"""
	#the pseudo R1 relation applies when the cardinality of clustA is 1
	#and the the cardinality of clustB is > 1
	#"""
	#if 1:	#(Cluster_Info_Dict[clustA]._GetCardinality() == 1) and (Cluster_Info_Dict[clustB]._GetCardinality() > 1):
		#"""
		#first of all, R4 relation should be allowed between clustA and clustB
		#"""
		#if (CheckAllowedRelnClusterPair(clustA, clustB, RELATION_R4) == True):
			#"""
			#next, we check whether R4 relation has significant frequency between this cluster pair
			#and also the pseudo R1 - R3 frequencies have significant proportion among the R4 frequency
			#"""
			#if (clustA < clustB):
				#clust_pair_key = (clustA, clustB)
			#else:
				#clust_pair_key = (clustB, clustA)
				
			#if clust_pair_key in Cluster_Pair_Info_Dict:
				#if (Cluster_Pair_Info_Dict[clust_pair_key]._CheckR4RelnSignificant()) and \
					#(Cluster_Pair_Info_Dict[clust_pair_key]._CheckPseudoR4RelnSignificant()):
					#print '\n *** the cluster pair ' + str(clust_pair_key) + '  has significant R4 proportion'

	#return

#-----------------------------------------------------
""" 
this function adds an edge between a pair of clusters (of taxa) 
it also updates the entries of reachability matrix 
"""
def Connect_ClusterPair(Reachability_Graph_Mat, clustA_reach_mat_idx, clustB_reach_mat_idx, reln_type, clustA, clustB):
	if (reln_type == RELATION_R1):
		"""
		adjust the clusters
		"""
		Cluster_Info_Dict[clustA]._AddRelnInstance(RELATION_R1, clustB)
		Cluster_Info_Dict[clustB]._AddRelnInstance(RELATION_R2, clustA)
		"""
		update the reachability matrix
		"""
		Reachability_Graph_Mat[clustA_reach_mat_idx][clustB_reach_mat_idx] = 1
	elif (reln_type == RELATION_R4):
		"""
		adjust the clusters
		"""
		Cluster_Info_Dict[clustA]._AddRelnInstance(RELATION_R4, clustB)
		Cluster_Info_Dict[clustB]._AddRelnInstance(RELATION_R4, clustA)    
		"""
		update the reachability matrix
		"""
		Reachability_Graph_Mat[clustA_reach_mat_idx][clustB_reach_mat_idx] = 2
		Reachability_Graph_Mat[clustB_reach_mat_idx][clustA_reach_mat_idx] = 2
		
	return Reachability_Graph_Mat

##-----------------------------------------------------        
#"""
#Associated function of Solve_MPP_PossibleR1R2 / RemovePossibleR2List_DiffParent
#resolves the inp_score_list if its length is > 1
#by discarding possible clusters

#@parameters:
	#curr_clust_idx: the cluster whose possible R2 relation clusters are checked
	#inp_score_list: contains the scores and cluster indices
	#ascending_order: if True, the scores in the inp_score_list are sorted in the ascending order
										#if False, scores are sorted in the descending order
#"""
#def ResolveScoreList(curr_clust_idx, inp_score_list, outfile, ascending_order=True):
	
	## this is a boolean flag 
	## it is True, only when at least one deletion operation from the inp_score_list is performed
	#flag = False
	
	#if (DEBUG_LEVEL >= 2):
		#fp = open(outfile, 'a')
	
	#if (ascending_order == True):
		#"""
		#here the first element of the score list has the minimum value
		#we note this value
		#we compare all the other score list elements according to this minimum value
		#and delete the clusters having values greater than this minimum value
		#"""
		#min_val = inp_score_list[0][1]
		#for i in range(1, len(inp_score_list)):
			#if (inp_score_list[i][1] > min_val):
				#"""
				#delete this entry from the possible R2 list
				#"""
				#target_delete_clust_idx = inp_score_list[i][0]
				#Cluster_Info_Dict[curr_clust_idx]._RemovePossibleR2(target_delete_clust_idx)
				#Cluster_Info_Dict[target_delete_clust_idx]._RemovePossibleR1(curr_clust_idx)
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Removed possible R2 edge from the cluster ' + str(curr_clust_idx) \
						#+ '  to the cluster: ' + str(target_delete_clust_idx))
				## set the flag variable
				#flag = True
	#else:
		#"""
		#here the first element of the score list has the maximum value
		#we note this value
		#we compare all the other score list elements according to this maximum value
		#and delete the clusters having values less than this maximum value
		#"""
		#max_val = inp_score_list[0][1]
		#for i in range(1, len(inp_score_list)):
			#if (inp_score_list[i][1] < max_val):
				#"""
				#delete this entry from the possible R2 list
				#"""
				#target_delete_clust_idx = inp_score_list[i][0]
				#Cluster_Info_Dict[curr_clust_idx]._RemovePossibleR2(target_delete_clust_idx)
				#Cluster_Info_Dict[target_delete_clust_idx]._RemovePossibleR1(curr_clust_idx)
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Removed possible R2 edge from the cluster ' + str(curr_clust_idx) \
						#+ '  to the cluster: ' + str(target_delete_clust_idx))
				## set the flag variable
				#flag = True
	
	#if (DEBUG_LEVEL >= 2):
		#fp.close()
	
	#return flag

##--------------------------------------------------------
#"""
#this function resolves the clusters having zero indegree
#such that it is either assigned as a descendant or as a possible R1 candidate to an existing cluster
#"""
#def Resolve_Cluster_ZeroIndegree(outfile):
	
	#if (DEBUG_LEVEL >= 2):
		#fp = open(outfile, 'a')
	
	##----------------------------------------------------------
	#"""
	#first check the clusters and try to resolve them using the R1 relation priority measure
	#"""
	#for cx in Cluster_Info_Dict:
		#"""
		#explore only those clusters cx having no in edge (R2 relation)
		#"""
		#cx_R2_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2)
		#if (len(cx_R2_list) == 0):
			#"""
			#explore the distinct possible R2 list
			#i.e. clusters which do not occur in the possible R1 list
			#"""
			#cx_possible_R2_list = Cluster_Info_Dict[cx]._GetPossibleR2List()
			#cx_distinct_possible_R2_list = Cluster_Info_Dict[cx]._GetDistinctPossibleR2List()
			#if (len(cx_possible_R2_list) > 0):
				
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n\n\n ***** In function - Resolve_Cluster_ZeroIndegree - Examining candidate cluster -- ' + str(cx) + \
					#' Its distinct possible R2 list: ' + str(cx_distinct_possible_R2_list) + \
						#'  Its possible R2 list: ' + str(cx_possible_R2_list))
				
				#"""
				#initialize empty list containing the score measure for individual x to current cx
				#"""
				#score_cx = []
				
				#for x in cx_possible_R2_list:	#cx_distinct_possible_R2_list:
					#"""
					#first find the R1 relation priority measure for individual clusters with respect to the cluster cx
					#the function returns two outputs:
					#no of supported couplets and the R1 relation priority sum
					#if the no of supported couplets > 0 then only we consider the score
					#"""
					#cc, dist1 = AvgR1RelnPriority(Cluster_Info_Dict[x]._GetSpeciesList(), Cluster_Info_Dict[cx]._GetSpeciesList())
					#if (cc > 0):
						#subl = [x, dist1]
						#score_cx.append(subl)
						#if (DEBUG_LEVEL >= 2):
							#fp.write('\n R1 relation priority from the cluster: ' + str(x) + '  is : ' + str(dist1))

				#"""
				#sort the scoring list in the descending order
				#"""
				#if (len(score_cx) > 1):
					#score_cx.sort(key=lambda x: x[1], reverse=True)
					
				#"""
				#if the first element of the scoring list contains a positive or zero quantity
				#then we place corresponding x as the parent of cx
				#"""
				#if (len(score_cx) > 0):
					#if (score_cx[0][1] >= 0):
						#x = score_cx[0][0]
						#Cluster_Info_Dict[x]._AddRelnInstance(RELATION_R1, cx)
						#Cluster_Info_Dict[cx]._AddRelnInstance(RELATION_R2, x)
						#if (DEBUG_LEVEL >= 2):
							#fp.write('\n ********* Established R1 relation (based on priority of R1 relation) from the cluster: ' \
								#+ str(x) + '  to the cluster : ' + str(cx))
	
	##----------------------------------------------------------
	
	##for cx in Cluster_Info_Dict:
		##"""
		##explore only those clusters cx having no in edge (R2 relation)
		##"""
		##cx_R2_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2)
		##if (len(cx_R2_list) == 0):
			##"""
			##explore the distinct possible R2 list
			##i.e. clusters which do not occur in the possible R1 list
			##"""
			##cx_possible_R2_list = Cluster_Info_Dict[cx]._GetPossibleR2List()
			##cx_distinct_possible_R2_list = Cluster_Info_Dict[cx]._GetDistinctPossibleR2List()
			##if (len(cx_possible_R2_list) > 0):
				
				##if (DEBUG_LEVEL >= 2):
					##fp.write('\n\n\n ***** In function - Resolve_Cluster_ZeroIndegree - Examining candidate cluster -- ' + str(cx) + \
					##' Its distinct possible R2 list: ' + str(cx_distinct_possible_R2_list) + \
						##'  Its possible R2 list: ' + str(cx_possible_R2_list))
				
				##for x in cx_possible_R2_list:	#cx_distinct_possible_R2_list:
					##"""
					##find the internode count measure for individual clusters with respect to the cluster cx
					##"""
					##dist = FindAvgInternodeCount(Cluster_Info_Dict[x]._GetSpeciesList(), Cluster_Info_Dict[cx]._GetSpeciesList(), 2, 1)
					##if (DEBUG_LEVEL >= 2):
						##fp.write('\n Internode distance from the cluster: ' + str(x) + '  is : ' + str(dist))

	#if (DEBUG_LEVEL >= 2):
		#fp.close()
		
	#return

##-----------------------------------------------------        
#"""
#this function computes the level difference (in terms of the ancestor / descendant relationship)
#from the src_clust to the dest_clust
#"""
#def ComputeClusterParentDistance(src_clust, dest_clust):
	#dist = 0
	#while (1):
		#if (src_clust == dest_clust):
			#break
		#dist_clust_parent_list = Cluster_Info_Dict[dest_clust]._GetClustRelnList(RELATION_R2)
		#if (len(dist_clust_parent_list) == 0):
			#return -1
		#else:
			#dist = dist + 1
			#dest_clust = dist_clust_parent_list[0]
	
	#return dist

##-----------------------------------------------------        
#"""
#Resolve the possible R2 lists for individual clusters cx
#Let the following cases are true:
	#1) cx_parent -> cx
	#2) x ---> cx
	#3) x_parent -> x
#We have to select whether we retain cx_parent as the ancestor of cx
#or we enforce x ---> cx as the primary connection to traverse cx
#"""
#def RemovePossibleR2List_DiffParent(outfile, Temporary_Cluster_Reln_List):

	#for cx in Cluster_Info_Dict:
		#"""
		#explore only those clusters cx such that there exists at least one cx_parent for cx
		#and also one x for cx
		#"""
		#cx_R2_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2)
		#cx_possible_R2_list = Cluster_Info_Dict[cx]._GetPossibleR2List()
		
		## add - sourya
		#if (len(cx_R2_list) == 0) and (len(cx_possible_R2_list) == 1):
			#Cluster_Info_Dict[cx]._AddFinalPossibleR2(cx_possible_R2_list[0])
		## end add - sourya
		
		#if (len(cx_R2_list) > 0) and (len(cx_possible_R2_list) > 0):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n\n\n ***** In function - RemovePossibleR2List_DiffParent - Examining candidate cluster -- ' + str(cx) + \
					#' Its R2 list: ' + str(cx_R2_list) + '  Its possible R2 list: ' + str(cx_possible_R2_list))

			#cx_parent = cx_R2_list[0]
			
			#cx_parent_internode_count = FindAvgInternodeCount(Cluster_Info_Dict[cx_parent]._GetSpeciesList(), \
				#Cluster_Info_Dict[cx]._GetSpeciesList(), 2, 1)
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n *** Internode count between cx_parent: ' + str(cx_parent) + ' to the current cluster : ' + str(cx_parent_internode_count))
			
			#"""
			#define an empty list whose individual elements consist of two fields:
			#1) x such that x-->cx
			#2) cx_parent = x_parent
			#otherwise, all other x information gets deleted
			#"""
			#score_cx = []
			
			#"""
			#define a scoring list for all x such that x has no parent
			#"""
			#score_x_no_parent = []

			#for x in cx_possible_R2_list:
				#if (Cluster_Info_Dict[x]._Get_Indegree() > 0):
					
					#x_parent = (Cluster_Info_Dict[x]._GetClustRelnList(RELATION_R2))[0]
					#dist = FindAvgInternodeCount(Cluster_Info_Dict[x]._GetSpeciesList(), Cluster_Info_Dict[cx]._GetSpeciesList(), 2, 1)
					#parent_dist = ComputeClusterParentDistance(cx_parent, x_parent)
					
					#if (DEBUG_LEVEL >= 2):
						#fp.write('\n cx_parent: ' + str(cx_parent) + ' x: ' + str(x) + ' x_parent: ' + str(x_parent))
						#fp.write('\n Internode count between the cluster ' + str(cx) + ' and ' + str(x) + '  is: ' + str(dist) + \
							#'   and the parent distance:  ' + str(parent_dist))
					
					#if (parent_dist < 0):
						##if (cx_parent != x_parent):
						#if (DEBUG_LEVEL >= 2):
							#fp.write('\n The cluster x: ' + str(x) + '  is not a descendant of cx_parent: ' + str(cx_parent))
					#else:
						#"""
						#here cx_parent = x_parent
						#insert the following three fields:
						#1) information of the cluster x
						#2) internode count "dist"
						#3) "parent_dist"
						#"""
						#subl = [x, dist, parent_dist]
						#score_cx.append(subl)
				
				#else:
					#dist = FindAvgInternodeCount(Cluster_Info_Dict[x]._GetSpeciesList(), Cluster_Info_Dict[cx]._GetSpeciesList(), 2, 1)
					#if (DEBUG_LEVEL >= 2):
						#fp.write('\n The cluster x: ' + str(x) + '  has no parent')
						#fp.write('\n Internode count between the cluster ' + str(cx) + ' and ' + str(x) + '  is: ' + str(dist))
					#subl = [x, dist]
					#score_x_no_parent.append(subl)
				
			#"""
			#sort the scoring list containing only the possible R2 relations in descending order, 
			#with respect to the distance between x_parent and cx_parent
			#"""
			#if (len(score_cx) > 1):
				#score_cx.sort(key=lambda x: x[1], reverse=True)

			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n --- Sorted scoring list (internode count distance) in descending order : ' + str(score_cx))
				
			#if (len(score_x_no_parent) > 1):
				#score_x_no_parent.sort(key=lambda x: x[1], reverse=True)

			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n --- Sorted list (internode count distance) of nodes having no parent -- in descending order : ' + str(score_x_no_parent))

			#if (DEBUG_LEVEL >= 2):
				#fp.close()

			## comment - sourya
			##if (len(score_cx) > 1):
				##ResolveScoreList(cx, score_cx, outfile, False)
			
			###--------------------------------------
			### add - sourya
			
			##if (len(score_cx) > 0):
				##"""
				##the first element should be the target possible R2 cluster
				##let this element be denoted as cz
				##so cz-->cx
				##"""
				##cz = score_cx[0][0]
				##"""
				##insert cz as the final possible R2 list candidate
				##"""
				##Cluster_Info_Dict[cx]._AddFinalPossibleR2(cz)
				##"""
				##we check if cx --> cz also
				##if not, then we delete the cx_parent from the parent list of cx 
				##"""
				##if cz not in Cluster_Info_Dict[cx]._GetPossibleR1List():
					##Cluster_Info_Dict[cx]._RemoveRelnInstance(RELATION_R2, cx_parent)
					##if (DEBUG_LEVEL >= 2):
						##fp = open(outfile, 'a')
						##fp.write('\n *** Here the cluster ' + str(cz) + '-->' + str(cx) + \
							##'  but the converse is not true --- removing ' + str(cx_parent) + ' from the in edge list of ' + str(cx))
						##fp.close()
				##else:
					##if (DEBUG_LEVEL >= 2):
						##fp = open(outfile, 'a')
						##fp.write('\n *** Here the cluster ' + str(cz) + '-->' + str(cx) + \
							##'  and the converse also holds --- ' + str(cx) + '-->' + str(cz) + '  so do not manipulate the parent list ')
						##fp.close()

			### end add - sourya
			###--------------------------------------

			##--------------------------------------
			## add - sourya
			
			#if (len(score_cx) > 0):
				#"""
				#the first element should be the target possible R2 cluster
				#let this element be denoted as cz
				#so cz-->cx
				#"""
				#cz = score_cx[0][0]
				#cz_internode_count = score_cx[0][1]
				#parent_dist = score_cx[0][2]
				
				#"""
				#we add the following condition:
				#only if the cz_internode_count exceeds the cx_parent_internode_count
				#then only we update the parent information
				#"""
				#if (cz_internode_count > cx_parent_internode_count) and (parent_dist >= 0):

					#"""
					#insert cz as the final possible R2 list candidate
					#"""
					#Cluster_Info_Dict[cx]._AddFinalPossibleR2(cz)
					#if (DEBUG_LEVEL >= 2):
						#fp = open(outfile, 'a')
						#fp.write('\n *** Established the connection ' + str(cz) + '-->' + str(cx))
						#fp.close()
					#if cz in Cluster_Info_Dict[cx]._GetPossibleR1List():
						#Cluster_Info_Dict[cz]._AddFinalPossibleR2(cx)
						#if (DEBUG_LEVEL >= 2):
							#fp = open(outfile, 'a')
							#fp.write('\n *** Also Established the connection ' + str(cx) + '-->' + str(cz))
							#fp.close()
					#"""
					#two conditions:
					#1) If parent_dist = 0 then parent of cz = parent of cz
					#in such a case, along with cz --> cx, we check whether cx --> cz also
					#if cx --> cz does not hold, we delete the cx_parent from the parent list of cx 
					
					#2) if parent_dist > 0, cx_parent -> cz_parent
					#so, we delete cx_parent from the parent list of cx 
					#"""
					#if (parent_dist == 0):
						#if cz not in Cluster_Info_Dict[cx]._GetPossibleR1List():
							#"""
							#note these clusters and their relations in the Temporary_Cluster_Reln_List
							#it will be used to remove cx_parent -> cx connection
							#"""
							#subl = [cx, cx_parent, RELATION_R2]
							#Temporary_Cluster_Reln_List.append(subl)
							#if (DEBUG_LEVEL >= 2):
								#fp = open(outfile, 'a')
								#fp.write('\n *** Here parent_dist = 0 ---- the cluster ' + str(cz) + '-->' + str(cx) + \
									#'  but the converse is not true --- removing ' + str(cx_parent) + ' from the in edge list of ' + str(cx))
								#fp.close()
						#else:
							#if (DEBUG_LEVEL >= 2):
								#fp = open(outfile, 'a')
								#fp.write('\n *** Here parent_dist = 0 ---- the cluster ' + str(cz) + '-->' + str(cx) + \
									#'  and the converse also holds --- ' + str(cx) + '-->' + str(cz) + '  so do not manipulate the parent list ')
								#fp.close()
					#else:
						#"""
						#note these clusters and their relations in the Temporary_Cluster_Reln_List
						#it will be used to remove cx_parent -> cx connection
						#"""
						#subl = [cx, cx_parent, RELATION_R2]
						#Temporary_Cluster_Reln_List.append(subl)
						#if (DEBUG_LEVEL >= 2):
							#fp = open(outfile, 'a')
							#fp.write('\n *** Here parent_dist > 0 ---- removing ' + str(cx_parent) + ' from the in edge list of ' + str(cx))
							#fp.close()

			## end add - sourya
			##--------------------------------------
			
	#return

##-----------------------------------------------------        
#"""
#This function checks the clusters having length of possible R2 list > 1
#In that case, it resolves the possible R2 list to have only one element
#the selection is carried out using the internode count based scoring mechanism
#"""
#def Solve_MPP_PossibleR1R2(Output_Text_File):
	##--------------------------------------------------------------
	#"""
	#now we check the clusters cx such that there exists at least two different clusters 
	#belonging to the possible R2 list of cx
	#"""
	#for cx in Cluster_Info_Dict:
		#"""
		#explore the possible R2 and R1 lists of the current cluster cx
		#"""
		#cx_R2_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2)
		#cx_possible_R2_list = Cluster_Info_Dict[cx]._GetPossibleR2List()
		#cx_possible_R1_list = Cluster_Info_Dict[cx]._GetPossibleR1List()
		#cx_FinalPossibleR2List = Cluster_Info_Dict[cx]._GetFinalPossibleR2List()
		
		## open the text file
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n\n\n ***** In function - Solve_MPP_PossibleR1R2 - cluster -- ' + str(cx) + \
				#'  cx_R2_list: ' + str(cx_R2_list) + '  cx_possible_R2_list: ' + str(cx_possible_R2_list) + \
					#'  cx_FinalPossibleR2List:  ' + str(cx_FinalPossibleR2List))
			#fp.close()
		
		## condition add - sourya
		#if (len(cx_R2_list) == 0) and (len(cx_possible_R2_list) == 1) and (len(cx_FinalPossibleR2List) == 0):
			#Cluster_Info_Dict[cx]._AddFinalPossibleR2(cx_possible_R2_list[0])
		## end condition - sourya
		
		##if (len(cx_possible_R2_list) > 1):
		#if (len(cx_R2_list) == 0) and (len(cx_possible_R2_list) > 1) and (len(cx_FinalPossibleR2List) == 0):

			## open the text file
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n\n\n ***** In function - Solve_MPP_PossibleR1R2 - Examining cluster -- ' + str(cx))
			
			#"""
			#initialize two scoring lists corresponding to the cluster cx
			#one will contain the scores of clusters cy such that cy <---> cx
			#other will contain the scores of clusters cy such that cy ---> cx
			#"""
			#score_list_cx_single = []
			#score_list_cx_both = []
			#"""
			#explore all candidate clusters cz belonging to cx_possible_R2_list
			#"""
			#for cz in cx_possible_R2_list:
				#cz_score = FindAvgInternodeCount(Cluster_Info_Dict[cz]._GetSpeciesList(), Cluster_Info_Dict[cx]._GetSpeciesList(), 2, 1)
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n --- element (possible R4 (R2) reln): ' + str(cz) + ' Internode count score: ' + str(cz_score))
				#temp_subl = [cz, cz_score]
				#if cz in cx_possible_R1_list:
					#score_list_cx_both.append(temp_subl)
				#else:
					#score_list_cx_single.append(temp_subl)
		
			#"""
			#sort score_list_cx_both in the ascending order
			#Note: Here we sort the list with respect to minimum internode count
			#"""
			#score_list_cx_both.sort(key=lambda x: x[1])	#, reverse=True)
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n --- Sorted score_list_cx_both (ascending order -- internode count distance) : ' + str(score_list_cx_both))
		
			#"""
			#sort score_list_cx_single in the ascending order
			#Note: Here we sort the list with respect to minimum internode count
			#"""
			#score_list_cx_single.sort(key=lambda x: x[1])	#, reverse=True)
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n --- Sorted score_list_cx_single (ascending order -- internode count distance) : ' + str(score_list_cx_single))
			
			## close the text file
			#if (DEBUG_LEVEL >= 2):
				#fp.close()

			##"""
			##process score_list_cx_both so that the first element in the list
			##containing the minimum internode count measure, gets remained
			##"""
			##if (len(score_list_cx_both) > 1):
				##ResolveScoreList(cx, score_list_cx_both, Output_Text_File)
			
			##"""
			##process score_list_cx_single so that the first element in the list
			##containing the minimum internode count measure, gets remained
			##"""
			##if (len(score_list_cx_single) > 1):
				##ResolveScoreList(cx, score_list_cx_single, Output_Text_File)
			
			#"""
			#if score_list_cx_single is empty, then the cluster should not have any possible R2 cluster
			#otherwise, scan through the sorted score_list_cx_both 
			#until an element in it contains either a parent cluster or an entry in its final allowed R2 list
			#"""
			#if (len(score_list_cx_single) > 0):
				#flag = False
				#if (len(score_list_cx_both) > 0):
					#for i in range(len(score_list_cx_both)):
						#if (flag == True):
							#break
						#x = score_list_cx_both[i][0]
						#x_parent_list = Cluster_Info_Dict[x]._GetClustRelnList(RELATION_R2)
						#x_allowed_R2_list = Cluster_Info_Dict[x]._GetFinalPossibleR2List()
						#if (len(x_parent_list) > 0) or (len(x_allowed_R2_list) > 0):
							#flag = True
						#"""
						#add this cluster in the final possible R2 list of the current cluster cx
						#"""
						#Cluster_Info_Dict[cx]._AddFinalPossibleR2(x)
			
				#if (flag == False):
					#if (len(score_list_cx_single) > 0):
						#for i in range(len(score_list_cx_single)):
							#if (flag == True):
								#break
							#x = score_list_cx_single[i][0]
							#x_parent_list = Cluster_Info_Dict[x]._GetClustRelnList(RELATION_R2)
							#x_allowed_R2_list = Cluster_Info_Dict[x]._GetFinalPossibleR2List()
							#if (len(x_parent_list) > 0) or (len(x_allowed_R2_list) > 0):
								#flag = True
							#"""
							#add this cluster in the final possible R2 list of the current cluster cx
							#"""
							#Cluster_Info_Dict[cx]._AddFinalPossibleR2(x)
					
	#return

##-----------------------------------------------------        
#"""
#solves no parent problem - checks the clusters having no parent at all 
#(empty R2 relation list) and empty possible R2 relation list as well
#in such a case, it finds a candidate cluster which will be placed as a R2 or 
#possible R2 relation candidate
#"""
#def Solve_NPP_NoPossibleR2(Reachability_Graph_Mat, Output_Text_File):

	#if (DEBUG_LEVEL >= 2):
		#fp = open(Output_Text_File, 'a')
		#fp.write('\n\n\n ***** Within function  --  Solve_NPP_NoPossibleR2 ******* ')
		#fp.close()

	##--------------------------------------------------------------
	#"""
	#first we check the clusters cx having following properties:
	#1) No cluster cy exists such that cy -> cx holds
	#2) There exists no clusters c1 such that c1--->cx 
	#"""
	#for cx in Cluster_Info_Dict:
		#"""
		#clusters having no indegree (no cluster is connected with R2 relation)
		#"""
		#if (Cluster_Info_Dict[cx]._Get_Indegree() == 0):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n ***** Examining cluster with zero indegree -- ' + str(cx))
				#fp.close()
			#if (DEBUG_LEVEL > 2):
				#Cluster_Info_Dict[cx]._PrintClusterInfo(cx, Output_Text_File)    
			
			#"""
			#possible R2 list of the cluster cx --- that is, collection of cz such that cz--->cx holds
			#"""
			#cx_possible_R2_list = Cluster_Info_Dict[cx]._GetPossibleR2List()
			#"""
			#here we check only those cases such that no cz exists
			#"""
			#if (len(cx_possible_R2_list) == 0):

				## open the text file
				#fp = open(Output_Text_File, 'a')
				
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n ***** The cluster has no element in possible R2 list as well -- ')
				
				#"""
				#it is a boolean flag, indicating that cx has no in edge or possible R2 cluster
				#"""
				#flag = False
				
				#"""
				#first find clusters x having the following property:
				#let cx->y (Relation R1, if holds)
				#and x->y (parent of y, excluding cx)
				#for such x, if x and cx are not related by any means, we include x->cx
				#"""
				#cx_R1_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R1)
				#if (len(cx_R1_list) > 0):
					#for y in cx_R1_list:
						#y_R2_list = Cluster_Info_Dict[y]._GetClustRelnList(RELATION_R2)
						#if (len(y_R2_list) > 0):
							#for x in y_R2_list:
								#if (x != cx):
									#if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(cx)] == 0) and \
										#(Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(cx)][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
										#if (DEBUG_LEVEL >= 2):
											#fp.write('\n We have ' + str(cx) + '->' + str(y) + '  and  ' + str(x) + '->' + str(y) + \
												#'  ===>>> So inserting ' + str(x) + '-> ' + str(cx) + '  since they are not related ')
										## establish x->cx
										#Reachability_Graph_Mat = Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
											#CURRENT_CLUST_IDX_LIST.index(cx), RELATION_R1, x, cx)
										## set the flag as well
										#flag = True
				
				#"""
				#if the flag is True at this point, then some cluster is assigned as the parent of cx
				#otherwise we explore the clusters x having the following property:
				#let cx--->y (Relation R1, if holds)
				#and x->y (parent of y, excluding cx)
				#for such x, if x and cx are not related by any means, we include x->cx
				#"""
				#if (flag == False):
					#cx_possible_R1_list = Cluster_Info_Dict[cx]._GetPossibleR1List()
					#if (len(cx_possible_R1_list) > 0):
						#for y in cx_possible_R1_list:
							#y_R2_list = Cluster_Info_Dict[y]._GetClustRelnList(RELATION_R2)
							#if (len(y_R2_list) > 0):
								#for x in y_R2_list:
									#if (x != cx):
										#if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(cx)] == 0) and \
											#(Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(cx)][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
											## establish x->cx
											#if (DEBUG_LEVEL >= 2):
												#fp.write('\n We have ' + str(cx) + '-->' + str(y) + '  and  ' + str(x) + '->' + str(y) + \
													#'  ===>>> So inserting ' + str(x) + '-> ' + str(cx) + '  since they are not related ')
											## establish x->cx
											#Reachability_Graph_Mat = Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
												#CURRENT_CLUST_IDX_LIST.index(cx), RELATION_R1, x, cx)
											#flag = True
				
				#"""
				#if the flag is True at this point, then some cluster is assigned as the parent of cx
				#otherwise, again we explore all the clusters such that cx--->y
				#this time we explore the clusters such that x--->y (and x is not cx)
				#we establish the relation x--->cx
				#"""
				#if (flag == False):
					#cx_possible_R1_list = Cluster_Info_Dict[cx]._GetPossibleR1List()
					#for y in cx_possible_R1_list:
						#y_possible_R2_list = Cluster_Info_Dict[y]._GetPossibleR2List()
						#if (len(y_possible_R2_list) > 0):
							#for x in y_possible_R2_list:
								#if (x != cx):
									#if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(cx)] == 0) and \
										#(Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(cx)][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
										#if (DEBUG_LEVEL >= 2):
											#fp.write('\n We have ' + str(cx) + '-->' + str(y) + '  and  ' + str(x) + '--->' + str(y) + \
												#'  ===>>> So inserting ' + str(x) + '----> ' + str(cx) + '  since they are not related ')
										## establish x--->cx
										#Reachability_Graph_Mat = Connect_PossibleR1Reln_ClusterPair(Reachability_Graph_Mat, x, cx)
										#flag = True
	
	
				## close the text file
				#fp.close()

	#if (DEBUG_LEVEL >= 2):
		#fp = open(Output_Text_File, 'a')
		#fp.write('\n\n\n ***** Out from the function  --  Solve_NPP_NoPossibleR2 ******* ')
		#fp.close()

	#return Reachability_Graph_Mat

#-----------------------------------------------------        
"""
this function solves multiple parent problem (C2)
by uniquely selecting one particular parent

Measures used: 
Support score measure (main)
internode count based scoring mechanism
"""
def SelectUniqueParent_Directed_SupportScore(ReachMat, Output_Text_File):

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')

	for cx in Cluster_Info_Dict:
		"""
		list of clusters cz such that cz->cx is satisfied
		"""
		cx_R2_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2)
		if (len(cx_R2_list) > 1):
			if (DEBUG_LEVEL >= 2):
				fp.write('\n\n ***** In function SelectUniqueParent_Directed_Part1 ---- Examining cluster -- ' + str(cx))
			
			"""
			initialize one scoring list corresponding to the cluster cx
			the scoring list contains scores corresponding to the relations cz->cx for all possible clusters cz
			"""
			score_list_cx = []
			for cz in cx_R2_list:
				"""
				second scoring measure : average internode count between cz and cx
				"""
				cz_score_internode = FindAvgDistanceMeasure(Cluster_Info_Dict[cz]._GetSpeciesList(), \
					Cluster_Info_Dict[cx]._GetSpeciesList(), 3, 2, 1)
				"""
				first scoring measure : support score of cz->cx connection 
				"""
				cz_score = GetRelnScore_ClusterPair(cx, cz, RELATION_R2)
				"""
				create a list containing the cluster cz information
				and also the support score and internode count measure
				"""
				temp_subl = [cz, cz_score, cz_score_internode]
				score_list_cx.append(temp_subl)
				if (DEBUG_LEVEL >= 2):
					fp.write('\n --- cluster (R2 reln): ' + str(cz) + '  Support score: ' + str(cz_score) + \
						'  Internode count score: ' + str(cz_score_internode))
			
			"""
			here we find the cluster having the maximum support score
			tie cases are broken with respect to lower internode count
			or in further, with respect to larger R1 relation list
			"""
			"""
			index of the score_list_cx containing the max score cluster
			"""
			max_score_clust_idx = 0
			for i in range(1, len(score_list_cx)):
				if (score_list_cx[i][1] > score_list_cx[max_score_clust_idx][1]):
					# higher support score
					max_score_clust_idx = i
				elif (score_list_cx[i][1] == score_list_cx[max_score_clust_idx][1]):
					# equal support score
					if (score_list_cx[i][2] < score_list_cx[max_score_clust_idx][2]):
						# lower internode count
						max_score_clust_idx = i
					elif (score_list_cx[i][2] == score_list_cx[max_score_clust_idx][2]):
						# equal internode count
						if (len(Cluster_Info_Dict[score_list_cx[i][0]]._GetClustRelnList(RELATION_R1)) > \
							len(Cluster_Info_Dict[score_list_cx[max_score_clust_idx][0]]._GetClustRelnList(RELATION_R1))):
							# higher length of the R1 cluster list
							max_score_clust_idx = i
			
			"""
			cluster corresponding to the maximum support score
			"""
			max_score_clust = score_list_cx[max_score_clust_idx][0]
			
			if (DEBUG_LEVEL >= 2):
				fp.write('\n Initially selected max_score_clust_idx: ' + str(max_score_clust_idx) + '  max_score_clust: ' + str(max_score_clust))
			
			#--------------------------------------------------------
			# this portion can be added - sourya
			#--------------------------------------------------------
			"""
			first we check whether there is a cluster "clust2"
			such that the R1 and R2 relation lists of clust2 are subsets of corresponding in max_score_clust
			"""
			subset_clust_list = []
			for i in range(len(score_list_cx)):
				if (i == max_score_clust_idx):
					continue
				if (CheckSubsetClust(max_score_clust, score_list_cx[i][0]) == True):
					templ = [score_list_cx[i][0], score_list_cx[i][1], score_list_cx[i][2]]
					subset_clust_list.append(templ)
					if (DEBUG_LEVEL >= 2):
						fp.write('\n Here the cluster ' + str(score_list_cx[i][0]) + \
							' has R1 & R2 lists which are subsets of the corresponding in ' + str(max_score_clust) + \
								'  appending in the list of subset_clust_list')
			
			if (DEBUG_LEVEL >= 2):
				fp.write('\n Formed subset_clust_list: ' + str(subset_clust_list))
			
			"""
			if the subset_clust_list is not empty:
			sort the list according to the support score measure (descending order)
			and select the cluster with the highest support score as the parent of the current cluster
			"""
			if (len(subset_clust_list) > 0):
				subset_clust_list.sort(key=lambda x: x[1], reverse=True)
				max_score_clust_new = subset_clust_list[0][0]
				if (DEBUG_LEVEL >= 2):
					fp.write('\n max_score_clust_new: ' + str(max_score_clust_new))

				"""
				first establish a directed edge connection from the cluster "max_score_clust"
				to the cluster "max_score_clust_new"
				"""
				Cluster_Info_Dict[max_score_clust]._AddRelnInstance(RELATION_R1, max_score_clust_new)
				Cluster_Info_Dict[max_score_clust_new]._AddRelnInstance(RELATION_R2, max_score_clust)
				
				"""
				then remove all the clusters except the "max_score_clust_new"
				from the parent list of cx
				"""
				for i in range(len(score_list_cx)):
					target_delete_clust = score_list_cx[i][0]
					if (target_delete_clust != max_score_clust_new):
						Remove_ClusterPairConn(cx, target_delete_clust, RELATION_R2)
						if (DEBUG_LEVEL >= 2):
							fp.write('\n Removed connection ' + str(target_delete_clust) + '->' + str(cx))
			
				#--------------------------------------------------------
				# end
				#--------------------------------------------------------
			else:
				"""
				otherwise, we have already selected max_score_clust_idx and max_score_clust
				remove other clusters from the parent list
				"""
				for i in range(len(score_list_cx)):
					target_delete_clust = score_list_cx[i][0]
					if (target_delete_clust != max_score_clust):
						Remove_ClusterPairConn(cx, target_delete_clust, RELATION_R2)
						if (DEBUG_LEVEL >= 2):
							fp.write('\n Removed connection ' + str(target_delete_clust) + '->' + str(cx))

	# close the text file
	if (DEBUG_LEVEL >= 2):
		fp.close()
		
	return 

#-----------------------------------------------------        
"""
this function solves multiple parent problem (C2)
by uniquely selecting one particular parent

Measures used: 
internode count based scoring mechanism (main)
"""
def SelectUniqueParent_Directed_Internode(ReachMat, Output_Text_File):
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')

	for cx in Cluster_Info_Dict:
		"""
		list of clusters cz such that cz->cx is satisfied
		"""
		cx_R2_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2)
		if (len(cx_R2_list) > 1):
			if (DEBUG_LEVEL >= 2):
				fp.write('\n\n ***** In function SelectUniqueParent_Directed_Part1 ---- Examining cluster -- ' + str(cx))
			
			"""
			initialize one scoring list corresponding to the cluster cx
			the scoring list contains scores corresponding to the relations cz->cx for all possible clusters cz
			"""
			score_list_cx = []
			for cz in cx_R2_list:
				"""
				second scoring measure : average internode count between cz and cx
				"""
				cz_score_internode = FindAvgDistanceMeasure(Cluster_Info_Dict[cz]._GetSpeciesList(), \
					Cluster_Info_Dict[cx]._GetSpeciesList(), 3, 2, 1)
				"""
				create a list containing the cluster cz information
				and also the support score and internode count measure
				"""
				temp_subl = [cz, cz_score_internode]
				score_list_cx.append(temp_subl)
				if (DEBUG_LEVEL >= 2):
					fp.write('\n --- cluster (R2 reln): ' + str(cz) + '  Internode count score: ' + str(cz_score_internode))
			
			"""
			here we find the cluster having the minimum internode count
			tie cases are broken with respect to larger R1 relation list
			"""
			min_internode_clust_idx = 0
			for i in range(1, len(score_list_cx)):
				if (score_list_cx[i][1] < score_list_cx[min_internode_clust_idx][1]):
					# lower internode count
					min_internode_clust_idx = i
				elif (score_list_cx[i][1] == score_list_cx[min_internode_clust_idx][1]):
					# equal internode count
					if (len(Cluster_Info_Dict[score_list_cx[i][0]]._GetClustRelnList(RELATION_R1)) > \
						len(Cluster_Info_Dict[score_list_cx[min_internode_clust_idx][0]]._GetClustRelnList(RELATION_R1))):
						# higher length of the R1 cluster list
						min_internode_clust_idx = i
			
			"""
			cluster corresponding to the minimum internode count
			"""
			min_internode_clust = score_list_cx[min_internode_clust_idx][0]
			
			if (DEBUG_LEVEL >= 2):
				fp.write('\n Initially selected min_internode_clust_idx: ' + str(min_internode_clust_idx) \
					+ '  min_internode_clust: ' + str(min_internode_clust))
			
			#--------------------------------------------------------
			# this portion can be added - sourya
			#--------------------------------------------------------
			"""
			first we check whether there is a cluster "clust2"
			such that the R1 and R2 relation lists of clust2 are subsets of corresponding in min_internode_clust
			"""
			subset_clust_list = []
			for i in range(len(score_list_cx)):
				if (i == min_internode_clust_idx):
					continue
				if (CheckSubsetClust(min_internode_clust, score_list_cx[i][0]) == True):
					templ = [score_list_cx[i][0], score_list_cx[i][1]]
					subset_clust_list.append(templ)
					if (DEBUG_LEVEL >= 2):
						fp.write('\n Here the cluster ' + str(score_list_cx[i][0]) + \
							' has R1 & R2 lists which are subsets of the corresponding in ' + str(min_internode_clust) + \
								'  appending in the list of subset_clust_list')
			
			if (DEBUG_LEVEL >= 2):
				fp.write('\n Formed subset_clust_list: ' + str(subset_clust_list))
			
			"""
			if the subset_clust_list is not empty:
			sort the list according to the lower internode count measure (ascending order)
			and select the cluster of the first index
			"""
			if (len(subset_clust_list) > 0):
				subset_clust_list.sort(key=lambda x: x[1])
				min_internode_clust_new = subset_clust_list[0][0]
				if (DEBUG_LEVEL >= 2):
					fp.write('\n min_internode_clust_new: ' + str(min_internode_clust_new))
				
				"""
				first establish a directed edge connection from the cluster "min_internode_clust"
				to the cluster "min_internode_clust_new"
				"""
				Cluster_Info_Dict[min_internode_clust]._AddRelnInstance(RELATION_R1, min_internode_clust_new)
				Cluster_Info_Dict[min_internode_clust_new]._AddRelnInstance(RELATION_R2, min_internode_clust)
				
				"""
				then remove all the clusters except the "min_internode_clust_new"
				from the parent list of cx
				"""
				for i in range(len(score_list_cx)):
					target_delete_clust = score_list_cx[i][0]
					if (target_delete_clust != min_internode_clust_new):
						Remove_ClusterPairConn(cx, target_delete_clust, RELATION_R2)
						if (DEBUG_LEVEL >= 2):
							fp.write('\n Removed connection ' + str(target_delete_clust) + '->' + str(cx))
				
				#--------------------------------------------------------
				# end
				#--------------------------------------------------------
			else:
				"""
				otherwise, we have already selected min_internode_clust_idx and min_internode_clust
				remove other clusters from the parent list
				"""
				for i in range(len(score_list_cx)):
					target_delete_clust = score_list_cx[i][0]
					if (target_delete_clust != min_internode_clust):
						Remove_ClusterPairConn(cx, target_delete_clust, RELATION_R2)
						if (DEBUG_LEVEL >= 2):
							fp.write('\n Removed connection ' + str(target_delete_clust) + '->' + str(cx))

	# close the text file
	if (DEBUG_LEVEL >= 2):
		fp.close()
		
	return 


##-----------------------------------------------------        
#"""
#this function returns True if subl1 is lower (should be placed before subl2 in a list sorted in ascending order) than subl2
#otherwise, it returns False
#subl has following items:
#1) cx - cluster index
#2) R1 freq
#3) Transitive R1 freq
#4) Support score
#"""
#def Lower_Elem(subl1, subl2):
	#cl1 = subl1[0]
	#r1_freq1 = subl1[1]
	#trans_r1_freq1 = subl1[2]
	#score1 = subl1[3]
	
	#cl2 = subl2[0]
	#r1_freq2 = subl2[1]
	#trans_r1_freq2 = subl2[2]
	#score2 = subl2[3]
	
	#if (score1 >= 0) and (score2 >= 0):
		#if (score1 < score2):
			#return True
		#elif (score1 > score2):
			#return False
	
	#if (trans_r1_freq1 < trans_r1_freq2):
		#return True
	#elif (trans_r1_freq1 > trans_r1_freq2):
		#return False
	
	#if (r1_freq1 < r1_freq2):
		#return True
	#elif (r1_freq1 > r1_freq2):
		#return False

	#return False

##-----------------------------------------------------        
#"""
#this function solves multiple parent problem (C2)
#by uniquely selecting one particular parent
#the selection is carried out using one of the following measure
	#internode count based scoring mechanism
#This function is applicable for directed in edge (Relation R2) only
#"""
#def SelectUniqueParent_Directed(Output_Text_File):

	## open the text file
	#if (DEBUG_LEVEL >= 2):
		#fp = open(Output_Text_File, 'a')

	#for cx in Cluster_Info_Dict:
		#"""
		#list of clusters belonging to the R2 list of cx
		#"""
		#cx_R2_list = Cluster_Info_Dict[cx]._GetClustRelnList(RELATION_R2)
		#if (len(cx_R2_list) > 1):
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n\n ***** In function SelectUniqueParent_Directed ---- Examining cluster -- ' + str(cx))
			
			#"""
			#initialize one scoring list corresponding to the cluster cx
			#"""
			#score_list_cx = []
			#"""
			#explore all the in edges (relation R2)
			#"""
			#for cz in cx_R2_list:
				#"""
				#obtain the R1 freq + transitive freq + support score measures 
				#"""
				#cz_cx_R1_freq = GetRelnFreq_ClusterPair(cz, cx, RELATION_R1)
				#cz_cx_R1_transitive_freq = GetRelnTransitiveFreq_ClusterPair(cz, cx, RELATION_R1)
				#cz_cx_R1_score = GetRelnScore_ClusterPair(cz, cx, RELATION_R1)
				#temp_subl = [cz, cz_cx_R1_freq, cz_cx_R1_transitive_freq, cz_cx_R1_score]
				#score_list_cx.append(temp_subl)
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n --- element (R2 reln): ' + str(cz) + \
						#'  R1 freq: ' + str(cz_cx_R1_freq) + '  Transitive R1 freq: ' + str(cz_cx_R1_transitive_freq) + \
							#'  Support score: ' + str(cz_cx_R1_score))
			
			#"""
			#now we sort the score_list_cx, according to the above mentioned parameters
			#we employ insertion sort
			#"""
			#for j in range(1, len(score_list_cx)):
				#"""
				#insert score_list_cx[j] into the sorted sequence score_list_cx[0...j-1]
				#"""
				#i = j - 1
				#while (i >= 0) and (Lower_Elem(score_list_cx[j], score_list_cx[i]) == True):
					#i = i - 1
				#"""
				#remove score_list_cx[j] and insert it in the (i+1)th position
				#"""
				#elem = score_list_cx.pop(j)
				#score_list_cx.insert(i+1, elem)
			
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n --- after sorting the scoring list (by insertion sort) corresponding to the cluster : ' + str(cx))
				#for i in range(len(score_list_cx)):
					#fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(score_list_cx[i][0]) + \
						#' R1 freq: ' + str(score_list_cx[i][1]) + '  Transitive R1 freq: ' + str(score_list_cx[i][2]) + \
							#'  Support score: ' + str(score_list_cx[i][3]))
			
			#"""
			#as the scoring list is based on ascending order, we have to select the last element of this list
			#as the parent of the target cluster
			#"""
			#for i in range(len(score_list_cx) - 1):
				#target_delete_clust_idx = score_list_cx[i][0]
				#Remove_ClusterPairConn(cx, target_delete_clust_idx, RELATION_R2)
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Removed out edge from the cluster ' + str(target_delete_clust_idx) + '  to the cluster: ' + str(cx))

	## close the text file
	#if (DEBUG_LEVEL >= 2):
		#fp.close()
		
	#return 

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
				if (Cluster_Info_Dict[i]._Get_Outdegree() > Cluster_Info_Dict[min_indeg_node_idx]._Get_Outdegree()):
					min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
					min_indeg_node_idx = i
					if (DEBUG_LEVEL >= 2):
						fp.write('\n Minimum indegree cluster: ' + str(i))
				#elif (valid_node_found == 1) and \
					#(len(Cluster_Info_Dict[i]._GetDistinctPossibleR2List()) == len(Cluster_Info_Dict[min_indeg_node_idx]._GetDistinctPossibleR2List())) \
					#and ((Cluster_Info_Dict[i]._Get_Outdegree() + len(Cluster_Info_Dict[i]._GetPossibleR1List())) > \
						#(Cluster_Info_Dict[min_indeg_node_idx]._Get_Outdegree() + len(Cluster_Info_Dict[min_indeg_node_idx]._GetPossibleR1List()))):    
					#min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
					#min_indeg_node_idx = i
					#if (DEBUG_LEVEL >= 2):
						#fp.write('\n Minimum indegree cluster: ' + str(i))


	##-------------------- new code - sourya

	#for i in clust_dict:
		#if (clust_dict[i]._GetExploredStatus() == 0):
			## we check the clusters which have not been explored yet
			#if (valid_node_found == 0):
				#min_indeg = Get_Non_Explored_Indegree(clust_dict, i)
				#min_DistPossR2 = Get_Non_Explored_DistinctPossibleR2Count(clust_dict, i)
				#min_outdeg = Get_Non_Explored_Outdegree(clust_dict, i)
				#min_DistPossR1 = Get_Non_Explored_PossibleR1Count(clust_dict, i)
				#min_indeg_node_idx = i
				#valid_node_found = 1
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Minimum non explored indegree cluster so far: ' + str(i))
			#else:
				#"""
				#here valid_node_found = 1
				#so at least one minimum value is already obtained 
				#"""
				#curr_indeg = Get_Non_Explored_Indegree(clust_dict, i)
				#curr_DistPossR2 = Get_Non_Explored_DistinctPossibleR2Count(clust_dict, i)
				#curr_outdeg = Get_Non_Explored_Outdegree(clust_dict, i)
				#curr_DistPossR1 = Get_Non_Explored_PossibleR1Count(clust_dict, i)
				#if (curr_indeg < min_indeg) or ((curr_indeg == min_indeg) and (curr_DistPossR2 < min_DistPossR2)) \
					#or ((curr_indeg == min_indeg) and (curr_DistPossR2 == min_DistPossR2) and (curr_outdeg > min_outdeg)):
					#min_indeg = curr_indeg
					#min_DistPossR2 = curr_DistPossR2
					#min_outdeg = curr_outdeg
					#min_DistPossR1 = curr_DistPossR1
					#min_indeg_node_idx = i
					#if (DEBUG_LEVEL >= 2):
						#fp.write('\n Minimum non explored indegree cluster: ' + str(i))

	#if (DEBUG_LEVEL >= 2):
		#fp.write('\n\n Out of the function --- Extract_Node_Min_Indeg \n\n ')
		#fp.close()
	## ---------------------- end new code - sourya
	
	return min_indeg_node_idx

##-----------------------------------------------------  
#def Get_Non_Explored_Indegree(clust_dict, i):
	#c = 0
	#for x in clust_dict[i]._GetClustRelnList(RELATION_R2):
		#if (clust_dict[x]._GetExploredStatus() == 0):
			#c = c + 1
	#return c

#def Get_Non_Explored_Outdegree(clust_dict, i):
	#c = 0
	#for x in clust_dict[i]._GetClustRelnList(RELATION_R1):
		#if (clust_dict[x]._GetExploredStatus() == 0):
			#c = c + 1
	#return c

#def Get_Non_Explored_DistinctPossibleR2Count(clust_dict, i):
	#c = 0
	#for x in clust_dict[i]._GetDistinctPossibleR2List():
		#if (clust_dict[x]._GetExploredStatus() == 0):
			#c = c + 1
	#return c

#def Get_Non_Explored_PossibleR1Count(clust_dict, i):
	#c = 0
	#for x in clust_dict[i]._GetPossibleR1List():
		#if (clust_dict[x]._GetExploredStatus() == 0):
			#c = c + 1
	#return c
##-----------------------------------------------------  
#"""
#this function separately contains the code for transitive reduction of the dashed edges (possible R1 / R2 relations)
#"""
#def Transitive_Reduction_Dashed(Clust_Possible_R1_Mat, Outfile):
	#no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
	
	## open the output text file
	#if (DEBUG_LEVEL >= 2):
		#fp = open(Outfile, 'a')

	#"""
	#transitive reduction for the following case: 
	#A--->B, B---->C / B<---->C, A---->C / A<---->C  ---------- remove A---->C / A<---->C
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
						## A--->C / A<---->C case
						#if (Clust_Possible_R1_Mat[i][k] == 1):	# and (Clust_Possible_R1_Mat[k][i] == 0):
							## comment - sourya - check
							##Clust_Possible_R1_Mat[i][k] = 0
							#Cluster_Info_Dict[clust_i]._RemovePossibleR1(clust_k)
							#Cluster_Info_Dict[clust_k]._RemovePossibleR2(clust_i)
							## add - sourya
							#Cluster_Info_Dict[clust_i]._RemovePossibleR2(clust_k)
							#Cluster_Info_Dict[clust_k]._RemovePossibleR1(clust_i)
							## end add - sourya
							#if (DEBUG_LEVEL >= 2):
								#fp.write('\n ---- ' + str(clust_i) + ' ----> ' + str(clust_j) + ', ' \
									#+ str(clust_j) + ' ----> / <----> ' + str(clust_k) + ', and ' \
									#+ str(clust_i) + ' -----> / <-----> ' + str(clust_k) \
										#+ ' --- so removing ' + str(clust_i) + ' ----> / <-----> ' + str(clust_k))

	#"""
	#transitive reduction for the following case: 
	#A--->B / A<--->B, B---->C, A---->C / A<---->C ---------- remove A---->C / A<---->C
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
							## add - sourya
							#Cluster_Info_Dict[clust_i]._RemovePossibleR2(clust_k)
							#Cluster_Info_Dict[clust_k]._RemovePossibleR1(clust_i)
							## end add - sourya
							#if (DEBUG_LEVEL >= 2):
								#fp.write('\n ---- ' + str(clust_i) + ' ----> / <----> ' + str(clust_j) \
									#+ ', ' + str(clust_j) + ' ---->' + str(clust_k) + ', and ' \
									#+ str(clust_i) + ' -----> / <----> ' + str(clust_k) \
										#+ ' --- so removing ' + str(clust_i) + ' ----> / <----> ' + str(clust_k))

	## close the output text file
	#if (DEBUG_LEVEL >= 2):
		#fp.close()

	#return

##-----------------------------------------------------  
#""" 
#this function performs transitive reduction of a graph for the dashed edges (possible R1 / R2 relations)
#"""
#def CompressDAG_Dashed(Clust_Possible_R1_Mat, Outfile):
	#no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
	
	## open the output text file
	#if (DEBUG_LEVEL >= 2):
		#fp = open(Outfile, 'a')

	## comment - sourya

	##"""
	##transitive reduction for the following case: 
	##A->B, B---->C, A->C ---------- remove A->C
	##"""
	##for j in range(no_of_clusters):
		##clust_j = CURRENT_CLUST_IDX_LIST[j]
		##for i in range(no_of_clusters):
			##if (i == j):
				##continue
			##clust_i = CURRENT_CLUST_IDX_LIST[i]
			### A->B case
			##if clust_j in Cluster_Info_Dict[clust_i]._GetClustRelnList(RELATION_R1):
				##for k in range(no_of_clusters):
					##if (i == k) or (j == k):
						##continue
					##clust_k = CURRENT_CLUST_IDX_LIST[k]
					### B---->C case (and also ensuring that C---->B does not hold)
					##if clust_k in Cluster_Info_Dict[clust_j]._GetPossibleR1List():
						##if clust_j not in Cluster_Info_Dict[clust_k]._GetPossibleR1List():
							### A->C case
							##if clust_k in Cluster_Info_Dict[clust_i]._GetClustRelnList(RELATION_R1):
								### remove the A->C connection
								##Remove_ClusterPairConn(clust_i, clust_k, RELATION_R1)
								##if (DEBUG_LEVEL >= 2):
									##fp.write('\n ---- ' + str(clust_i) + ' -> ' + str(clust_j) + ', ' \
										##+ str(clust_j) + ' ----> ' + str(clust_k) + ', and ' \
										##+ str(clust_i) + ' -> ' + str(clust_k) \
											##+ ' --- so removing ' + str(clust_i) + ' -> ' + str(clust_k))


	##"""
	##transitive reduction for the following case: 
	##A---->B, B->C, A->C ---------- remove A->C
	##"""
	##for j in range(no_of_clusters):
		##clust_j = CURRENT_CLUST_IDX_LIST[j]
		##for i in range(no_of_clusters):
			##if (i == j):
				##continue
			##clust_i = CURRENT_CLUST_IDX_LIST[i]
			### A---->B case (and also ensuring that B---->A does not hold)
			##if clust_j in Cluster_Info_Dict[clust_i]._GetPossibleR1List():
				##if clust_i not in Cluster_Info_Dict[clust_j]._GetPossibleR1List():
					##for k in range(no_of_clusters):
						##if (i == k) or (j == k):
							##continue
						##clust_k = CURRENT_CLUST_IDX_LIST[k]
						### B->C case
						##if clust_k in Cluster_Info_Dict[clust_j]._GetClustRelnList(RELATION_R1):
							### A->C case
							##if clust_k in Cluster_Info_Dict[clust_i]._GetClustRelnList(RELATION_R1):
								### remove the A->C connection
								##Remove_ClusterPairConn(clust_i, clust_k, RELATION_R1)
								##if (DEBUG_LEVEL >= 2):
									##fp.write('\n ---- ' + str(clust_i) + ' ----> ' \
										##+ str(clust_j) + ', ' + str(clust_j) + ' -> ' + str(clust_k) + ', and ' \
										##+ str(clust_i) + ' -> ' + str(clust_k) \
											##+ ' --- so removing ' + str(clust_i) + ' -> ' + str(clust_k))

	## end comment - sourya

	##------------------------------------
	## add - sourya
	#"""
	#transitive reduction for the following case: 
	#A->B, B---->C / B<---->C, A--->C / A<--->C ---------- remove A--->C / A<--->C
	#"""
	#for j in range(no_of_clusters):
		#clust_j = CURRENT_CLUST_IDX_LIST[j]
		#for i in range(no_of_clusters):
			#if (i == j):
				#continue
			#clust_i = CURRENT_CLUST_IDX_LIST[i]
			## A->B case
			#if clust_j in Cluster_Info_Dict[clust_i]._GetClustRelnList(RELATION_R1):
				#for k in range(no_of_clusters):
					#if (i == k) or (j == k):
						#continue
					#clust_k = CURRENT_CLUST_IDX_LIST[k]
					## B---->C / B<---->C case
					#if clust_k in Cluster_Info_Dict[clust_j]._GetPossibleR1List():
						## A--->C / A<--->C case
						#if clust_k in Cluster_Info_Dict[clust_i]._GetPossibleR1List():
							## remove A--->C / A<--->C connection
							#Cluster_Info_Dict[clust_i]._RemovePossibleR1(clust_k)
							#Cluster_Info_Dict[clust_k]._RemovePossibleR2(clust_i)
							## add - sourya
							#Cluster_Info_Dict[clust_i]._RemovePossibleR2(clust_k)
							#Cluster_Info_Dict[clust_k]._RemovePossibleR1(clust_i)
							## end add - sourya
							#if (DEBUG_LEVEL >= 2):
								#fp.write('\n ---- ' + str(clust_i) + ' -> ' \
									#+ str(clust_j) + ', ' + str(clust_j) + ' ---> / <---> ' + str(clust_k) + ', and ' \
									#+ str(clust_i) + ' ---> / <---> ' + str(clust_k) \
										#+ ' --- so removing ' + str(clust_i) + ' ---> / <---> ' + str(clust_k))

	#"""
	#transitive reduction for the following case: 
	#A--->B / A<--->B, B->C, A--->C / A<--->C ---------- remove A--->C / A<--->C
	#"""
	#for j in range(no_of_clusters):
		#clust_j = CURRENT_CLUST_IDX_LIST[j]
		#for i in range(no_of_clusters):
			#if (i == j):
				#continue
			#clust_i = CURRENT_CLUST_IDX_LIST[i]
			## A---->B / A<--->B case
			#if clust_j in Cluster_Info_Dict[clust_i]._GetPossibleR1List():
				#for k in range(no_of_clusters):
					#if (i == k) or (j == k):
						#continue
					#clust_k = CURRENT_CLUST_IDX_LIST[k]
					## B->C case
					#if clust_k in Cluster_Info_Dict[clust_j]._GetClustRelnList(RELATION_R1):
						## A--->C / A<--->C case
						#if clust_k in Cluster_Info_Dict[clust_i]._GetPossibleR1List():
							## remove A--->C / A<--->C connection
							#Cluster_Info_Dict[clust_i]._RemovePossibleR1(clust_k)
							#Cluster_Info_Dict[clust_k]._RemovePossibleR2(clust_i)
							## add - sourya
							#Cluster_Info_Dict[clust_i]._RemovePossibleR2(clust_k)
							#Cluster_Info_Dict[clust_k]._RemovePossibleR1(clust_i)
							## end add - sourya
							#if (DEBUG_LEVEL >= 2):
								#fp.write('\n ---- ' + str(clust_i) + ' ---> / <----> ' \
									#+ str(clust_j) + ', ' + str(clust_j) + ' -> ' + str(clust_k) + ', and ' \
									#+ str(clust_i) + ' ---> / <---> ' + str(clust_k) \
										#+ ' --- so removing ' + str(clust_i) + ' ---> / <---> ' + str(clust_k))

	## end add - sourya
	##------------------------------------

	## close the output text file
	#if (DEBUG_LEVEL >= 2):
		#fp.close()

	#return

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

