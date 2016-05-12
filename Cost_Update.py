#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#----------------------------------------------------------------
"""
this function fills the support score values for individual cluster pairs
also accumulates the frequency values for individual relations (except the relation R3)
"""
def Fill_Queue_Support_Score_Cluster_Pair(outfile):
	
	#if (DEBUG_LEVEL >= 2):
		#fp = open(outfile, 'a')
		#fp.write('\n\n\n **** Within function --- Fill_Queue_Support_Score_Cluster_Pair **** \n\n')

	no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
	
	for i in range(no_of_clusters - 1):
		clust1 = CURRENT_CLUST_IDX_LIST[i]
		taxa_list1 = Cluster_Info_Dict[clust1]._GetSpeciesList()
		for j in range(i+1, no_of_clusters):
			clust2 = CURRENT_CLUST_IDX_LIST[j]
			taxa_list2 = Cluster_Info_Dict[clust2]._GetSpeciesList()
			R1_freq, R2_freq, R4_freq, R1_score, R2_score, R4_score, couplet_count, \
				r1_reln_allowed, r2_reln_allowed = GetFreqScore_ClusterPair(taxa_list1, taxa_list2)
			
			"""
			now for this cluster pair, insert the support score and frequency values (along with the 
			cluster pair information) to the support score queue
			"""
			if (R1_freq > 0) and (r1_reln_allowed == True):
				sublist = [clust1, clust2, RELATION_R1, R1_freq, R1_score]	# ((R1_score * 1.0) / couplet_count)]
				Queue_Score_Cluster_Pair.append(sublist)
			
			if (R2_freq > 0) and (r2_reln_allowed == True):
				sublist = [clust1, clust2, RELATION_R2, R2_freq, R2_score]	# ((R2_score * 1.0) / couplet_count)]
				Queue_Score_Cluster_Pair.append(sublist)
			
			if (R4_freq > 0):
				sublist = [clust1, clust2, RELATION_R4, R4_freq, R4_score]	# ((R4_score * 1.0) / couplet_count)]
				Queue_Score_Cluster_Pair.append(sublist)

	#---------------------------------------------------------------------
	#if (DEBUG_LEVEL >= 2):
		#fp.write('\n\n\n **** Out from the function --- Fill_Queue_Support_Score_Cluster_Pair **** \n\n')
		#fp.close()
	
	return

#----------------------------------------------------------------
"""
here we fill the support score queues with the couplet relations and the support score values
from which the DAG will be constructed
"""
def Fill_Support_Score_Queues_Couplet_Based_Reln_R3():
	
	""" 
	here, we process all the couplets
	individual couplets lead to individual cluster pairs
	for individual relations between a couplet, corresponding edge between the cluster pair is established
	"""  
	for l in TaxaPair_Reln_Dict:
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

		"""
		first set the priority values of this couplet
		"""
		TaxaPair_Reln_Dict[l]._SetConnPrVal()

		""" 
		calculate the support score and priority measures for individual couplets
		and for individual relations
		it is the product of frequency and the priority measures 
		"""
		TaxaPair_Reln_Dict[l]._SetCostMetric()
		
		for reln_type in [RELATION_R1, RELATION_R2, RELATION_R4, RELATION_R3]:
			curr_reln_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(reln_type)
			if (curr_reln_freq > 0):
				if (reln_type == RELATION_R4):
					TaxaPair_Reln_Dict[l]._AddAllowedReln(reln_type)

				if (reln_type == RELATION_R1) or (reln_type == RELATION_R2):
					if (TaxaPair_Reln_Dict[l]._GetEdgeWeight(RELATION_R4) == 0) \
						or (TaxaPair_Reln_Dict[l]._GetEdgeWeight(reln_type) >= (0.7 * max_freq)):
						TaxaPair_Reln_Dict[l]._AddAllowedReln(reln_type)

				if (reln_type == RELATION_R3):
					if (len(TaxaPair_Reln_Dict[l]._GetAllowedRelnList()) <= 1) and (TaxaPair_Reln_Dict[l]._CheckTargetRelnConsensus(RELATION_R3) == True):
						TaxaPair_Reln_Dict[l]._AddAllowedReln(reln_type)
						"""
						only put the support score corresponding to the relation R3
						the support score queue has 5 elements - modification
						"""
						if ((consensus_reln_type == RELATION_R3) and (single_edge_exist == 1)):
							# add - sourya
							# we also add R1, R2 relations as the allowed relations
							TaxaPair_Reln_Dict[l]._AddAllowedReln(RELATION_R1)
							TaxaPair_Reln_Dict[l]._AddAllowedReln(RELATION_R2)
							# end add - sourya
							reln_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(consensus_reln_type)
							support_score = TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(consensus_reln_type)
							sublist = [l[0], l[1], consensus_reln_type, reln_freq, support_score]
							Queue_Score_R3_SingleReln.append(sublist)
						else:
							if (TaxaPair_Reln_Dict[l]._CheckTargetRelnMajorityConsensus(RELATION_R3) == True):
								reln_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(RELATION_R3)
								if (reln_freq > 1):
									support_score = TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(RELATION_R3)
									sublist = [l[0], l[1], RELATION_R3, reln_freq, support_score]
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
this function is used for cost based sorting of the inp_queue in the unweighted supertree
this function compares two elements of the heap
and returns the index with the higher support score measure
"""
def Higher_Score_Value(inp_queue, i, j):
	key1 = (inp_queue[i][0], inp_queue[i][1])
	reln1 = inp_queue[i][2]
	freq1 = inp_queue[i][3]
	score1 = inp_queue[i][4]
	key2 = (inp_queue[j][0], inp_queue[j][1])
	reln2 = inp_queue[j][2]
	freq2 = inp_queue[j][3]
	score2 = inp_queue[j][4]
	
	# case A - if both scores are strictly positive
	if (score1 > 0) or (score2 > 0):
		if (score1 < score2):
			return j
		elif (score1 > score2):
			return i
		else:	#if (score1 == score2):
			"""
			for tie case of cost
			"""
			if (freq1 < freq2):
				return j      
			elif (freq2 < freq1):
				return i
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
		elif (reln1 == RELATION_R4) and (reln2 != RELATION_R4):
			return i
		elif (reln1 != RELATION_R4) and (reln2 == RELATION_R4):
			return j
    
	return i
#------------------------------------------------------------------------
  
"""
this function exchanges two elements in the heap
"""
def Exchange_Elem(inp_queue, i, j):
	"""
	here r corresponds to different fields (of the sublists)
	of the entries in the support score queue
	as individual queue entry consists of 5 different fields
	the copy operation should take care of the 5 fields
	"""
	for r in range(5):
		temp_key = inp_queue[i][r]
		inp_queue[i][r] = inp_queue[j][r]
		inp_queue[j][r] = temp_key

"""
maintain max heap property
note: heap_size may not be the actual length of the queue
but the working length (on which the remaining sorting operation will take place)
"""
def Max_Heapify(inp_queue, idx, heap_size):
	l = Left(idx)
	r = Right(idx)
	if (l < heap_size) and (Higher_Score_Value(inp_queue, idx, l) == l):
		largest_idx = l
	else:
		largest_idx = idx
	if (r < heap_size) and (Higher_Score_Value(inp_queue, largest_idx, r) == r):
		largest_idx = r
	if (largest_idx != idx):
		Exchange_Elem(inp_queue, idx, largest_idx)
		Max_Heapify(inp_queue, largest_idx, heap_size)

"""
extract the current maximum and also pop the element from the heap
"""
def Heap_Extract_Max(inp_queue):
	if (len(inp_queue) < 1):
		print 'underflow of max priority queue'
	"""
	1st element is the maximum
	"""
	max_elem = list(inp_queue[0])
	"""
	replace the first element with the last element of the queue
	"""
	Exchange_Elem(inp_queue, 0, len(inp_queue) - 1)
	"""
	delete the last element of the queue
	"""
	del inp_queue[len(inp_queue) - 1]
	heap_size = len(inp_queue)
	"""
	call the max_heapify function to maintain the heap property
	0 is the starting index of the list storing the heap structure
	"""
	Max_Heapify(inp_queue, 0, heap_size)
	return max_elem
  
"""
this function builds the priority queue (max heap property)
"""
def Build_Max_Heap(inp_queue):
	heap_size = len(inp_queue)
	for idx in range(int(len(inp_queue) / 2), -1, -1):
		Max_Heapify(inp_queue, idx, heap_size)
    
"""
this is the heap sort algorithm
"""
def Sort_Priority_Queue(inp_queue):
	Build_Max_Heap(inp_queue)
	heap_size = len(inp_queue)

