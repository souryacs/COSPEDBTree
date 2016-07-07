#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

##--------------------------------------------
#"""
#this function checks all the cluster pairs having at least R1 / R2 relation between them
#"""
#def Check_Transitive_R1R2_Allowed_Connections(outfile):
	
	#if (DEBUG_LEVEL >= 2):
		#fp = open(outfile, 'a')
	
	#for cx in Cluster_Pair_Info_Dict:
		#cx_allowed_reln_list = Cluster_Pair_Info_Dict[cx]._GetPossibleRelnList()
		#"""
		#check if R1 relation is allowed between this cluster pair
		#in such a case, find the following transitive connectivity:
		#Let A = cx[0] and B = cx[1]
		#R1(A,B) means A->B
		#accumulate the frequency of D->A and D->B if both are allowed
		#also accumulate the frequency of A->E and B->E if both are allowed
		#"""
		#A = cx[0]
		#B = cx[1]

		#if RELATION_R1 in cx_allowed_reln_list:
			
			#cx_r1_freq = Cluster_Pair_Info_Dict[cx]._GetFreq(RELATION_R1)
			#Cluster_Pair_Info_Dict[cx]._AddTransitiveFreq(RELATION_R1, cx_r1_freq)
			
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n\n *** Cluster pair: ' + str(cx) + '  has R1 as its allowed relation --- with frequency: ' + str(cx_r1_freq))
			
			#Common_R2_Clust_list = list(set(Cluster_Info_Dict[A]._GetInitialList(RELATION_R2)).intersection(set(Cluster_Info_Dict[B]._GetInitialList(RELATION_R2))))
			
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n Its Common_R2_Clust_list (corresponds to the list of D) ' + str(Common_R2_Clust_list)) 
			
			#for D in Common_R2_Clust_list:
				#D_A_R1_freq = GetRelnFreq_ClusterPair(D, A, RELATION_R1)
				#D_B_R1_freq = GetRelnFreq_ClusterPair(D, B, RELATION_R1)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveFreq(RELATION_R1, D_A_R1_freq)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveFreq(RELATION_R1, D_B_R1_freq)
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Adding R1 frequencies of values ' + str(D_A_R1_freq) + '  and  ' + str(D_B_R1_freq) + \
						#'  corresponding to the cluster ' + str(D))
			
			#Common_R1_Clust_list = list(set(Cluster_Info_Dict[A]._GetInitialList(RELATION_R1)).intersection(set(Cluster_Info_Dict[B]._GetInitialList(RELATION_R1))))
			
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n Its Common_R1_Clust_list (corresponds to the list of E) ' + str(Common_R1_Clust_list)) 
			
			#for E in Common_R1_Clust_list:
				#A_E_R1_freq = GetRelnFreq_ClusterPair(A, E, RELATION_R1)
				#B_E_R1_freq = GetRelnFreq_ClusterPair(B, E, RELATION_R1)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveFreq(RELATION_R1, A_E_R1_freq)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveFreq(RELATION_R1, B_E_R1_freq)
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Adding R1 frequencies of values ' + str(A_E_R1_freq) + '  and  ' + str(B_E_R1_freq) + \
						#'  corresponding to the cluster ' + str(E))
					
			#"""
			#now print the transitively implied R1 frequency of this cluster pair
			#"""
			#cx_trans_R1_freq = Cluster_Pair_Info_Dict[cx]._GetTransitiveFreq(RELATION_R1)
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n ===>>> Cluster pair: ' + str(cx) + '  original R1 freq: ' + str(cx_r1_freq) + \
					#'  transitively implied R1 freq: ' + str(cx_trans_R1_freq))
			
		#"""
		#check if R2 relation is allowed between this cluster pair
		#in such a case, find the following transitive connectivity:
		#Let A = cx[0] and B = cx[1]
		#R2(A,B) means A<-B
		#accumulate the frequency of D<-A and D<-B if both are allowed
		#also accumulate the frequency of A<-E and B<-E if both are allowed
		#"""
		#if RELATION_R2 in cx_allowed_reln_list:
			
			#cx_r2_freq = Cluster_Pair_Info_Dict[cx]._GetFreq(RELATION_R2)
			#Cluster_Pair_Info_Dict[cx]._AddTransitiveFreq(RELATION_R2, cx_r2_freq)
			
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n\n *** Cluster pair: ' + str(cx) + '  has R2 as its allowed relation --- with frequency: ' + str(cx_r2_freq))
			
			#Common_R1_Clust_list = list(set(Cluster_Info_Dict[A]._GetInitialList(RELATION_R1)).intersection(set(Cluster_Info_Dict[B]._GetInitialList(RELATION_R1))))
			
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n Its Common_R1_Clust_list (corresponds to the list of D) ' + str(Common_R1_Clust_list)) 
			
			#for D in Common_R1_Clust_list:
				#D_A_R2_freq = GetRelnFreq_ClusterPair(D, A, RELATION_R2)
				#D_B_R2_freq = GetRelnFreq_ClusterPair(D, B, RELATION_R2)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveFreq(RELATION_R2, D_A_R2_freq)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveFreq(RELATION_R2, D_B_R2_freq)
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Adding R2 frequencies of values ' + str(D_A_R2_freq) + '  and  ' + str(D_B_R2_freq) + \
						#'  corresponding to the cluster ' + str(D))
			
			#Common_R2_Clust_list = list(set(Cluster_Info_Dict[A]._GetInitialList(RELATION_R2)).intersection(set(Cluster_Info_Dict[B]._GetInitialList(RELATION_R2))))

			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n Its Common_R2_Clust_list (corresponds to the list of E) ' + str(Common_R2_Clust_list)) 

			#for E in Common_R2_Clust_list:
				#A_E_R2_freq = GetRelnFreq_ClusterPair(A, E, RELATION_R2)
				#B_E_R2_freq = GetRelnFreq_ClusterPair(B, E, RELATION_R2)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveFreq(RELATION_R2, A_E_R2_freq)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveFreq(RELATION_R2, B_E_R2_freq)
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Adding R2 frequencies of values ' + str(A_E_R2_freq) + '  and  ' + str(B_E_R2_freq) + \
						#'  corresponding to the cluster ' + str(E))
	
			#"""
			#now print the transitively implied R2 frequency of this cluster pair
			#"""
			#cx_trans_R2_freq = Cluster_Pair_Info_Dict[cx]._GetTransitiveFreq(RELATION_R2)
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n ===>>> Cluster pair: ' + str(cx) + '  original R2 freq: ' + str(cx_r2_freq) + \
					#'  transitively implied R2 freq: ' + str(cx_trans_R2_freq))
	
	#if (DEBUG_LEVEL >= 2):	
		#fp.close()
	
	#return

##--------------------------------------------
#"""
#this function checks all the cluster pairs having at least R1 / R2 relation between them
#"""
#def Check_Transitive_R1R2_Allowed_Connections(outfile):
	
	#if (DEBUG_LEVEL >= 2):
		#fp = open(outfile, 'a')
	
	#for cx in Cluster_Pair_Info_Dict:
		#cx_allowed_reln_list = Cluster_Pair_Info_Dict[cx]._GetPossibleRelnList()
		#"""
		#check if R1 relation is allowed between this cluster pair
		#in such a case, find the following transitive connectivity:
		#Let A = cx[0] and B = cx[1]
		#R1(A,B) means A->B
		#accumulate the support score of D->A and D->B if both are allowed
		#also accumulate the support score of A->E and B->E if both are allowed
		#"""
		#A = cx[0]
		#B = cx[1]

		#if RELATION_R1 in cx_allowed_reln_list:
			
			#cx_r1_freq = Cluster_Pair_Info_Dict[cx]._GetFreq(RELATION_R1)
			#cx_r1_score = Cluster_Pair_Info_Dict[cx]._GetSupportScore(RELATION_R1)
			#Cluster_Pair_Info_Dict[cx]._AddTransitiveScore(RELATION_R1, cx_r1_score)
			
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n\n *** Cluster pair: ' + str(cx) + '  has R1 as its allowed relation --- with support score: ' + str(cx_r1_score))
			
			#Common_R2_Clust_list = list(set(Cluster_Info_Dict[A]._GetInitialList(RELATION_R2)).intersection(set(Cluster_Info_Dict[B]._GetInitialList(RELATION_R2))))
			
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n Its Common_R2_Clust_list (corresponds to the list of D) ' + str(Common_R2_Clust_list)) 
			
			#for D in Common_R2_Clust_list:
				#D_A_R1_score = GetRelnScore_ClusterPair(D, A, RELATION_R1)
				#D_B_R1_score = GetRelnScore_ClusterPair(D, B, RELATION_R1)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveScore(RELATION_R1, D_A_R1_score)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveScore(RELATION_R1, D_B_R1_score)
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Adding R1 support scores of values ' + str(D_A_R1_score) + '  and  ' + str(D_B_R1_score) + \
						#'  corresponding to the cluster ' + str(D))
			
			#Common_R1_Clust_list = list(set(Cluster_Info_Dict[A]._GetInitialList(RELATION_R1)).intersection(set(Cluster_Info_Dict[B]._GetInitialList(RELATION_R1))))
			
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n Its Common_R1_Clust_list (corresponds to the list of E) ' + str(Common_R1_Clust_list)) 
			
			#for E in Common_R1_Clust_list:
				#A_E_R1_score = GetRelnScore_ClusterPair(A, E, RELATION_R1)
				#B_E_R1_score = GetRelnScore_ClusterPair(B, E, RELATION_R1)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveScore(RELATION_R1, A_E_R1_score)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveScore(RELATION_R1, B_E_R1_score)
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Adding R1 support scores of values ' + str(A_E_R1_score) + '  and  ' + str(B_E_R1_score) + \
						#'  corresponding to the cluster ' + str(E))
					
			#"""
			#now print the transitively implied R1 frequency of this cluster pair
			#"""
			#cx_trans_R1_score = Cluster_Pair_Info_Dict[cx]._GetTransitiveScore(RELATION_R1)
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n ===>>> Cluster pair: ' + str(cx) + '  original R1 support score: ' + str(cx_r1_score) + \
					#'  transitively implied R1 support score: ' + str(cx_trans_R1_score))
			
			#"""
			#now insert the transitive support score and the R1 frequency in the support score queue
			#"""
			#subl = [A, B, RELATION_R1, cx_r1_freq, cx_trans_R1_score]
			#Queue_Score_Cluster_Pair.append(subl)
			
		#"""
		#check if R2 relation is allowed between this cluster pair
		#in such a case, find the following transitive connectivity:
		#Let A = cx[0] and B = cx[1]
		#R2(A,B) means A<-B
		#accumulate the frequency of D<-A and D<-B if both are allowed
		#also accumulate the frequency of A<-E and B<-E if both are allowed
		#"""
		#if RELATION_R2 in cx_allowed_reln_list:
			
			#cx_r2_freq = Cluster_Pair_Info_Dict[cx]._GetFreq(RELATION_R2)
			#cx_r2_score = Cluster_Pair_Info_Dict[cx]._GetSupportScore(RELATION_R2)
			#Cluster_Pair_Info_Dict[cx]._AddTransitiveScore(RELATION_R2, cx_r2_score)
			
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n\n *** Cluster pair: ' + str(cx) + '  has R2 as its allowed relation --- with support score: ' + str(cx_r2_score))
			
			#Common_R1_Clust_list = list(set(Cluster_Info_Dict[A]._GetInitialList(RELATION_R1)).intersection(set(Cluster_Info_Dict[B]._GetInitialList(RELATION_R1))))
			
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n Its Common_R1_Clust_list (corresponds to the list of D) ' + str(Common_R1_Clust_list)) 
			
			#for D in Common_R1_Clust_list:
				#D_A_R2_score = GetRelnScore_ClusterPair(D, A, RELATION_R2)
				#D_B_R2_score = GetRelnScore_ClusterPair(D, B, RELATION_R2)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveScore(RELATION_R2, D_A_R2_score)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveScore(RELATION_R2, D_B_R2_score)
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Adding R2 support scores of values ' + str(D_A_R2_score) + '  and  ' + str(D_B_R2_score) + \
						#'  corresponding to the cluster ' + str(D))
			
			#Common_R2_Clust_list = list(set(Cluster_Info_Dict[A]._GetInitialList(RELATION_R2)).intersection(set(Cluster_Info_Dict[B]._GetInitialList(RELATION_R2))))

			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n Its Common_R2_Clust_list (corresponds to the list of E) ' + str(Common_R2_Clust_list)) 

			#for E in Common_R2_Clust_list:
				#A_E_R2_score = GetRelnScore_ClusterPair(A, E, RELATION_R2)
				#B_E_R2_score = GetRelnScore_ClusterPair(B, E, RELATION_R2)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveScore(RELATION_R2, A_E_R2_score)
				#Cluster_Pair_Info_Dict[cx]._AddTransitiveScore(RELATION_R2, B_E_R2_score)
				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n Adding R2 support scores of values ' + str(A_E_R2_score) + '  and  ' + str(B_E_R2_score) + \
						#'  corresponding to the cluster ' + str(E))
	
			#"""
			#now print the transitively implied R2 frequency of this cluster pair
			#"""
			#cx_trans_R2_score = Cluster_Pair_Info_Dict[cx]._GetTransitiveScore(RELATION_R2)
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n ===>>> Cluster pair: ' + str(cx) + '  original R2 support score: ' + str(cx_r2_score) + \
					#'  transitively implied R2 support score: ' + str(cx_trans_R2_score))
	
			#"""
			#now insert the transitive support score and the R2 frequency in the support score queue
			#"""
			#subl = [A, B, RELATION_R2, cx_r2_freq, cx_trans_R2_score]
			#Queue_Score_Cluster_Pair.append(subl)
	
	#if (DEBUG_LEVEL >= 2):	
		#fp.close()
	
	#return


#----------------------------------------------------------------
"""
this function checks the cluster pairs having at least one of the relations R1 or R2 as their allowed relation
and forms a queue of support scores using those relations
"""
def Form_Conflict_Support_Score_Queue():
	for cx in Cluster_Pair_Info_Dict:
		cx_allowed_reln_list = Cluster_Pair_Info_Dict[cx]._GetPossibleRelnList()
		if RELATION_R1 in cx_allowed_reln_list:
			#subl = [cx[0], cx[1], RELATION_R1, Cluster_Pair_Info_Dict[cx]._GetFreq(RELATION_R1), \
				#Cluster_Pair_Info_Dict[cx]._GetTransitiveFreq(RELATION_R1), \
					#Cluster_Pair_Info_Dict[cx]._GetSupportScore(RELATION_R1)]
			subl = [cx[0], cx[1], RELATION_R1, Cluster_Pair_Info_Dict[cx]._GetFreq(RELATION_R1), \
				Cluster_Pair_Info_Dict[cx]._GetSupportScore(RELATION_R1)]
			Queue_Score_Cluster_Pair.append(subl)
			
		if RELATION_R2 in cx_allowed_reln_list:
			#subl = [cx[0], cx[1], RELATION_R2, Cluster_Pair_Info_Dict[cx]._GetFreq(RELATION_R2), \
				#Cluster_Pair_Info_Dict[cx]._GetTransitiveFreq(RELATION_R2), \
					#Cluster_Pair_Info_Dict[cx]._GetSupportScore(RELATION_R2)]
			subl = [cx[0], cx[1], RELATION_R2, Cluster_Pair_Info_Dict[cx]._GetFreq(RELATION_R2), \
				Cluster_Pair_Info_Dict[cx]._GetSupportScore(RELATION_R2)]
			Queue_Score_Cluster_Pair.append(subl)
			
		# modify - sourya
		# here we add relation R4 only if it is the only relation allowed between this cluster pair
		if (RELATION_R4 in cx_allowed_reln_list):	# and (len(cx_allowed_reln_list) > 1):
			if (RELATION_R1 not in cx_allowed_reln_list) and (RELATION_R2 not in cx_allowed_reln_list):
				subl = [cx[0], cx[1], RELATION_R4, Cluster_Pair_Info_Dict[cx]._GetFreq(RELATION_R4), \
					Cluster_Pair_Info_Dict[cx]._GetSupportScore(RELATION_R4)]
				Queue_Score_Cluster_Pair.append(subl)
		
	return

##----------------------------------------------------------------
#"""
#this function checks the cluster pairs having only R4 relation as their allowed relation
#and forms a queue of support scores using those relations
#"""
#def Form_NonConflict_Support_Score_Queue(outfile):
	
	#if (DEBUG_LEVEL > 2):
		#fp = open(outfile, 'a')
	
	#for cx in Cluster_Pair_Info_Dict:
		#cx_allowed_reln_list = Cluster_Pair_Info_Dict[cx]._GetPossibleRelnList()
		#if (len(cx_allowed_reln_list) == 1) and (RELATION_R4 in cx_allowed_reln_list):
			#clust1 = cx[0]
			#clust2 = cx[1]
			#subl = [clust1, clust2, RELATION_R4, Cluster_Pair_Info_Dict[cx]._GetFreq(RELATION_R4), \
				#Cluster_Pair_Info_Dict[cx]._GetSupportScore(RELATION_R4)]
			#Queue_Score_Cluster_Pair_NonConflict.append(subl)
			#if (DEBUG_LEVEL > 2):
				#fp.write('\n Cluster pair: ' + str(cx) + '  has only R4 as its allowed relation -- adding it in Queue_Score_Cluster_Pair_NonConflict')

	#if (DEBUG_LEVEL > 2):
		#fp.close()
	
	#return

#----------------------------------------------------------------
"""
this function fills the support score values for individual cluster pairs
also accumulates the frequency values for individual relations (except the relation R3)
"""
def Initialize_Cluster_Pair(outfile):
	
	no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
	
	for i in range(no_of_clusters - 1):
		clust1 = CURRENT_CLUST_IDX_LIST[i]
		taxa_list1 = Cluster_Info_Dict[clust1]._GetSpeciesList()
		for j in range(i+1, no_of_clusters):
			clust2 = CURRENT_CLUST_IDX_LIST[j]
			taxa_list2 = Cluster_Info_Dict[clust2]._GetSpeciesList()
			R1_freq, R2_freq, R3_freq, R4_freq, couplet_count, pseudo_r4_r1_count, pseudo_r4_r2_count, \
						pseudo_r4_r3_count = GetFreqScore_ClusterPair(taxa_list1, taxa_list2)
			
			if (DEBUG_LEVEL > 2):
				fp = open(outfile, 'a')
				fp.write('\n **** Examine cluster pair ' + str(clust1) + ' and ' + str(clust2) + ' R1_freq: ' + str(R1_freq) + \
					' R2_freq: ' + str(R2_freq) + ' R3_freq: ' + str(R3_freq) + ' R4_freq: ' + str(R4_freq))
				fp.close()
			
			"""
			initialize this cluster pair only if they are related by one of the 
			four relations r1 to r4 in at least one input tree
			"""
			support_tree_count = R1_freq + R2_freq + R3_freq + R4_freq
			if (support_tree_count > 0):
			
				"""
				initialize the cluster pair key and corresponding relation instance
				"""
				clust_pair_key = (clust1, clust2)
				if clust_pair_key not in Cluster_Pair_Info_Dict:
					Cluster_Pair_Info_Dict.setdefault(clust_pair_key, Cluster_Pair(R1_freq, R2_freq, R3_freq, R4_freq, \
						couplet_count, pseudo_r4_r1_count, pseudo_r4_r2_count, pseudo_r4_r3_count))
				
				"""
				frequency of the consensus relation
				"""
				max_freq = max(R1_freq, R2_freq, R3_freq, R4_freq)
				
				"""
				add R4 relation if it has positive frequency
				"""
				if (R4_freq > 0):
					if 1:	#(R4_freq == max_freq) or (R4_freq >= (PERCENT_MAX_FREQ1 * max_freq)):
						Cluster_Pair_Info_Dict[clust_pair_key]._AddPossibleReln(RELATION_R4)
				
				"""
				add R1 or R2 relations only if they are either consensus 
				or having a significant frequency compared to the consensus relation
				Note: We maintain two different allowed relation list
				1) _AddAllPossibleRelnList function inserts any significant relation
				2) _AddPossibleReln inserts R1 / R2 relation only if it has greater frequency than corresponding R2 / R1 relation
				"""
				if (R1_freq > 0):
					if ((R1_freq >= PERCENT_MAX_FREQ1 * max_freq) or ((R1_freq + R3_freq) >= (PERCENT_MAX_FREQ2 * max_freq))):
						if (R1_freq >= R2_freq):
							Cluster_Pair_Info_Dict[clust_pair_key]._AddPossibleReln(RELATION_R1)
				
				if (R2_freq > 0):
					if ((R2_freq >= PERCENT_MAX_FREQ1 * max_freq) or ((R2_freq + R3_freq) >= (PERCENT_MAX_FREQ1 * max_freq))):
						if (R2_freq >= R1_freq):
							Cluster_Pair_Info_Dict[clust_pair_key]._AddPossibleReln(RELATION_R2)
				
				#"""
				#extract the allowed relation list between this pair of cluster
				#"""
				#poss_reln_list = Cluster_Pair_Info_Dict[clust_pair_key]._GetPossibleRelnList()

				## add - sourya
				#"""
				#here we check if the relations R1, R2, and R4 are all included in the allowed relation list
				#in such a case, we check the following conditions
				#"""
				#if (len(poss_reln_list) == 3):
					#min_freq = min(R1_freq, R2_freq, R3_freq, R4_freq)
					#if (R4_freq == min_freq):
						#"""
						#frequency of R4 relation is the lowest among all the relations
						#so, discard it from the set of allowed relations
						#"""
						#Cluster_Pair_Info_Dict[clust_pair_key]._RemovePossibleReln(RELATION_R4)
					#else:
						#"""
						#here R4 relation is surely included
						#and in fact, discard both R1 and R2 relations from the set of allowed relations
						#since both the relations have identical frequencies
						#"""
						#Cluster_Pair_Info_Dict[clust_pair_key]._RemovePossibleReln(RELATION_R1)
						#Cluster_Pair_Info_Dict[clust_pair_key]._RemovePossibleReln(RELATION_R2)
				## end add - sourya
				
				"""
				again re-extract the allowed relation list between this cluster pair
				"""
				poss_reln_list = Cluster_Pair_Info_Dict[clust_pair_key]._GetPossibleRelnList()
				if RELATION_R1 in poss_reln_list:
					Cluster_Info_Dict[clust1]._AppendInitialList(RELATION_R1, clust2)
					Cluster_Info_Dict[clust2]._AppendInitialList(RELATION_R2, clust1)
				if RELATION_R2 in poss_reln_list:
					Cluster_Info_Dict[clust1]._AppendInitialList(RELATION_R2, clust2)
					Cluster_Info_Dict[clust2]._AppendInitialList(RELATION_R1, clust1)
				
				#"""
				#we differentiate between two cases:
				#1) When R4 relation is the only allowed relation between this pair of cluster
				#2) Otherwise, one of the relations R1 / R2 is allowed
				#"""
				#if (len(poss_reln_list) == 1) and (RELATION_R4 in poss_reln_list):
					#subl = [clust1, clust2, RELATION_R4, R4_freq, Cluster_Pair_Info_Dict[clust_pair_key]._GetSupportScore(RELATION_R4)]
					#Queue_Score_Cluster_Pair_NonConflict.append(subl)
				#else:
					#"""
					#we include the R1 / R2 relation information in the Queue_Score_Cluster_Pair
					#provided they belong to the allowed relation list between this pair of cluster
					#"""
					#if RELATION_R1 in poss_reln_list:
						#subl = [clust1, clust2, RELATION_R1, R1_freq, Cluster_Pair_Info_Dict[clust_pair_key]._GetSupportScore(RELATION_R1)]
						#Queue_Score_Cluster_Pair.append(subl)
					#if RELATION_R2 in poss_reln_list:
						#subl = [clust1, clust2, RELATION_R2, R2_freq, Cluster_Pair_Info_Dict[clust_pair_key]._GetSupportScore(RELATION_R2)]
						#Queue_Score_Cluster_Pair.append(subl)
					##if RELATION_R4 in poss_reln_list:
						##subl = [clust1, clust2, RELATION_R4, R4_freq]
						##Queue_Score_Cluster_Pair.append(subl)
			
	return

#----------------------------------------------------------------
"""
here we fill the support score queues with the couplet relations and the support score values
from which the DAG will be constructed
"""
def Fill_Support_Score_Queues_Couplet_Based():
	
	""" 
	here, we process all the couplets
	individual couplets lead to individual cluster pairs
	for individual relations between a couplet, corresponding edge between the cluster pair is established
	"""  
	for l in TaxaPair_Reln_Dict:
		r1_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(RELATION_R1)
		r2_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(RELATION_R2)
		r3_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(RELATION_R3)
		r4_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(RELATION_R4)
		
		"""
		single_edge_exist: if TRUE, means that only one type of relation 
		is supported (with respect to input trees) between this couplet
		detection of it during setting the priority values of different relations
		basically, we are looking for the consensus relation
		"""
		single_edge_exist_list = TaxaPair_Reln_Dict[l]._CheckNonConflictingCouplet()
		single_edge_exist = single_edge_exist_list[0]
		consensus_reln_type = single_edge_exist_list[1]
		
		"""
		maximum frequency among all 4 relations
		"""
		max_freq = TaxaPair_Reln_Dict[l]._GetConsensusFreq()

		#"""
		#first set the priority values of this couplet
		#"""
		#TaxaPair_Reln_Dict[l]._SetConnPrVal()

		#""" 
		#calculate the support score and priority measures for individual couplets
		#and for individual relations
		#it is the product of frequency and the priority measures 
		#"""
		#TaxaPair_Reln_Dict[l]._SetCostMetric()
		
		"""
		for other types of couplets, not included in above criterion
		"""
		for reln_type in [RELATION_R4, RELATION_R1, RELATION_R2, RELATION_R3]:
			curr_reln_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(reln_type)
			if (curr_reln_freq > 0):
				if (reln_type == RELATION_R4):
					"""
					non zero frequency of R4 relation is considered
					"""
					TaxaPair_Reln_Dict[l]._AddAllowedReln(reln_type)
				elif (reln_type == RELATION_R1) or (reln_type == RELATION_R2):
					if (r4_freq == 0) or (curr_reln_freq >= (0.7 * max_freq)):	# or (curr_reln_freq == max_freq_nonR3):
					#if (curr_reln_freq >= (0.7 * max_freq)) or (curr_reln_freq == max_freq_nonR3):
						"""
						either zero frequency of R4 relation, or sufficient frequency of R1/R2 relations are considered
						for their inclusion in the 
						"""
						TaxaPair_Reln_Dict[l]._AddAllowedReln(reln_type)
				else:	#if (reln_type == RELATION_R3):
					if (len(TaxaPair_Reln_Dict[l]._GetAllowedRelnList()) <= 1) \
						and (TaxaPair_Reln_Dict[l]._CheckTargetRelnConsensus(reln_type) == True):
						"""
						first add relation R3 in the allowed relation list
						"""
						TaxaPair_Reln_Dict[l]._AddAllowedReln(reln_type)
						"""
						differentiate between two cases
						1) When R3 is the sole relation between this couplet
						2) When R3 is the majority consensus relation between this couplet
						"""
						if (consensus_reln_type == RELATION_R3) and (single_edge_exist == 1):
							"""
							we also add R1, R2, R4 relations as the allowed relations
							"""
							TaxaPair_Reln_Dict[l]._AddAllowedReln(RELATION_R1)
							TaxaPair_Reln_Dict[l]._AddAllowedReln(RELATION_R2)
							TaxaPair_Reln_Dict[l]._AddAllowedReln(RELATION_R4)
							"""
							add the support score relation in the queue
							"""
							reln_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(consensus_reln_type)
							#support_score = TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(consensus_reln_type)
							sublist = [l[0], l[1], consensus_reln_type, reln_freq, reln_freq]	#, support_score]
							Queue_Score_R3_SingleReln.append(sublist)
						elif (TaxaPair_Reln_Dict[l]._CheckTargetRelnMajorityConsensus(RELATION_R3) == True) and (r3_freq > 1):
							"""
							we impose one condition that the no of trees bearing this R3 relation should be > 1
							add the support score relation in the queue
							"""
							#support_score = TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(RELATION_R3)
							sublist = [l[0], l[1], RELATION_R3, r3_freq, r3_freq]	#, support_score]
							Queue_Score_R3_MajCons.append(sublist)
	
	return

#------------------------------------------------------
""" this code section implements the max priority queue """
#------------------------------------------------------
# parent node of current node
def Parent(idx):
	return int((idx-1) / 2)

# left child node of current node  
def Left(idx):
	return (2*idx+1)

# right child node of current node
def Right(idx):
	return (2*idx+2)

#------------------------------------------------------------------------
# version 1X - 15th Sep, 2015 
#-------------------------
"""
this function is used for sorting the support score queue (priority queue structure)

Parameters:
1) inp_queue: Input priority queue
	can be either Nonconflicting support score queue or the conflicting support score queue
2) i, j denotes indices of two different elements within the priority queue
these two elements are to be compared
3) inp_no = can be either 1 or 2
	if 1- denotes the Nonconflicting support score queue and corresponding configuration
	if 2 - denotes the conflicting support score queue and corresponding configuration

output:
returns index with the higher support score measure
"""
def Higher_Score_Value(inp_queue, i, j, inp_no=1):
	key1 = (inp_queue[i][0], inp_queue[i][1])
	reln1 = inp_queue[i][2]
	freq1 = inp_queue[i][3]

	key2 = (inp_queue[j][0], inp_queue[j][1])
	reln2 = inp_queue[j][2]
	freq2 = inp_queue[j][3]

	if (inp_no == 1):
		score1 = inp_queue[i][4]
		score2 = inp_queue[j][4]
	else:
		trans_freq1 = inp_queue[i][4]
		trans_freq2 = inp_queue[j][4]
		score1 = inp_queue[i][5]
		score2 = inp_queue[j][5]
	
	if (inp_no == 1):
		"""
		function when non-conflicting support score queue is used
		"""
		# case A - if both scores are strictly positive
		if (score1 > 0) or (score2 > 0):
			if (score1 < score2):
				return j
			elif (score1 > score2):
				return i
			else:	#if (score1 == score2):
				if (freq1 < freq2):
					return j      
				elif (freq2 < freq1):
					return i
				#elif ((reln1 == RELATION_R1) or (reln1 == RELATION_R2)) and ((reln2 != RELATION_R1) and (reln2 != RELATION_R2)):
					#return i
				#elif ((reln2 == RELATION_R1) or (reln2 == RELATION_R2)) and ((reln1 != RELATION_R1) and (reln1 != RELATION_R2)):
					#return j
				elif (reln1 == RELATION_R4) and (reln2 != RELATION_R4):
					return i
				elif (reln1 != RELATION_R4) and (reln2 == RELATION_R4):
					return j
		else:
			"""
			here one or both scores are either zero or negative
			"""
			if (freq1 < freq2):
				return j      
			elif (freq2 < freq1):
				return i     
			#elif ((reln1 == RELATION_R1) or (reln1 == RELATION_R2)) and ((reln2 != RELATION_R1) and (reln2 != RELATION_R2)):
				#return i
			#elif ((reln2 == RELATION_R1) or (reln2 == RELATION_R2)) and ((reln1 != RELATION_R1) and (reln1 != RELATION_R2)):
				#return j
			elif (reln1 == RELATION_R4) and (reln2 != RELATION_R4):
				return i
			elif (reln1 != RELATION_R4) and (reln2 == RELATION_R4):
				return j
		
		# default return condition
		return i
		
	else:
		"""
		function when conflicting support score queue is used
		"""
		# case A - if both scores are strictly positive
		if (score1 > 0) or (score2 > 0):
			if (score1 < score2):
				return j
			elif (score1 > score2):
				return i
			else:	#if (score1 == score2):
				if (trans_freq1 < trans_freq2):
					return j
				elif (trans_freq1 > trans_freq2):
					return i
				elif (freq1 < freq2):
					return j      
				elif (freq2 < freq1):
					return i
				elif ((reln1 == RELATION_R1) or (reln1 == RELATION_R2)) and ((reln2 != RELATION_R1) and (reln2 != RELATION_R2)):
					return i
				elif ((reln2 == RELATION_R1) or (reln2 == RELATION_R2)) and ((reln1 != RELATION_R1) and (reln1 != RELATION_R2)):
					return j
				#elif (reln1 == RELATION_R4) and (reln2 != RELATION_R4):
					#return i
				#elif (reln1 != RELATION_R4) and (reln2 == RELATION_R4):
					#return j
		else:
			"""
			here one or both scores are either zero or negative
			"""
			if (trans_freq1 < trans_freq2):
				return j
			elif (trans_freq1 > trans_freq2):
				return i
			elif (freq1 < freq2):
				return j      
			elif (freq2 < freq1):
				return i     
			elif ((reln1 == RELATION_R1) or (reln1 == RELATION_R2)) and ((reln2 != RELATION_R1) and (reln2 != RELATION_R2)):
				return i
			elif ((reln2 == RELATION_R1) or (reln2 == RELATION_R2)) and ((reln1 != RELATION_R1) and (reln1 != RELATION_R2)):
				return j
			#elif (reln1 == RELATION_R4) and (reln2 != RELATION_R4):
				#return i
			#elif (reln1 != RELATION_R4) and (reln2 == RELATION_R4):
				#return j
		
		# default return condition
		return i
	
#------------------------------------------------------------------------
"""
this function exchanges two elements in the heap
"""
def Exchange_Elem(inp_queue, i, j, inp_no=1):
	"""
	here r corresponds to different fields (of the sublists)
	of the entries in the support score queue
	as individual queue entry consists of 5 different fields
	the copy operation should take care of the 5 fields
	"""
	if (inp_no == 1):
		total_elem = 5
	else:
		total_elem = 6
	
	# swap individual fields
	for r in range(total_elem):	#(4):
		temp_key = inp_queue[i][r]
		inp_queue[i][r] = inp_queue[j][r]
		inp_queue[j][r] = temp_key

	return

#-----------------------------------------------
"""
maintain max heap property
note: heap_size may not be the actual length of the queue
but the working length (on which the remaining sorting operation will take place)
"""
def Max_Heapify(inp_queue, idx, heap_size, inp_no=1):
	l = Left(idx)
	r = Right(idx)
	if (l < heap_size) and (Higher_Score_Value(inp_queue, idx, l, inp_no) == l):
		largest_idx = l
	else:
		largest_idx = idx
	if (r < heap_size) and (Higher_Score_Value(inp_queue, largest_idx, r, inp_no) == r):
		largest_idx = r
	if (largest_idx != idx):
		Exchange_Elem(inp_queue, idx, largest_idx, inp_no)
		Max_Heapify(inp_queue, largest_idx, heap_size, inp_no)

#-----------------------------------------------
"""
extract the current maximum and also pop the element from the heap
"""
def Heap_Extract_Max(inp_queue, inp_no=1):
	if (len(inp_queue) < 1):
		print 'underflow of max priority queue'
	"""
	1st element is the maximum
	"""
	max_elem = list(inp_queue[0])
	"""
	replace the first element with the last element of the queue
	"""
	Exchange_Elem(inp_queue, 0, len(inp_queue) - 1, inp_no)
	"""
	delete the last element of the queue
	"""
	del inp_queue[len(inp_queue) - 1]
	heap_size = len(inp_queue)
	"""
	call the max_heapify function to maintain the heap property
	0 is the starting index of the list storing the heap structure
	"""
	Max_Heapify(inp_queue, 0, heap_size, inp_no)
	return max_elem

#-----------------------------------------------
"""
this function builds the priority queue (max heap property)
"""
def Build_Max_Heap(inp_queue, inp_no=1):
	heap_size = len(inp_queue)
	for idx in range(int(len(inp_queue) / 2), -1, -1):
		Max_Heapify(inp_queue, idx, heap_size, inp_no)

#-----------------------------------------------
"""
this is the heap sort algorithm
parameters:
1) inp_queue: Input priority queue
	can be either Nonconflicting support score queue or the conflicting support score queue
2) inp_no = can be either 1 or 2 (default 1)
	if 1- denotes the Nonconflicting support score queue and corresponding configuration
	if 2 - denotes the conflicting support score queue and corresponding configuration

"""
def Sort_Priority_Queue(inp_queue, inp_no=1):
	Build_Max_Heap(inp_queue, inp_no)
	heap_size = len(inp_queue)

