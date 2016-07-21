#!/usr/bin/env python

import Header
from Header import *
import Cluster_Manage
from Cluster_Manage import *

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
		
	return

#-----------------------------------------------------
""" 
this function updates the transitive closure of the cluster of nodes
on inclusion of a new edge between a pair of clusters 
"""
def TransClosUpd(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, reln_type, outfile):
	if (reln_type == RELATION_R1) or (reln_type == RELATION_R4):
		src_reach_mat_idx = nodeA_reach_mat_idx
		dest_reach_mat_idx = nodeB_reach_mat_idx
	elif (reln_type == RELATION_R2):
		src_reach_mat_idx = nodeB_reach_mat_idx
		dest_reach_mat_idx = nodeA_reach_mat_idx
	else:
		return
		
	src_taxa_clust_idx = CURRENT_CLUST_IDX_LIST[src_reach_mat_idx]
	dest_taxa_clust_idx = CURRENT_CLUST_IDX_LIST[dest_reach_mat_idx]
		
	if (reln_type == RELATION_R1) or (reln_type == RELATION_R2):
		"""
		for A->B connection
		if D->A exists
		then establish D->B
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R2):
			if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][dest_reach_mat_idx] == 0):
				Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
					dest_reach_mat_idx, RELATION_R1, x, dest_taxa_clust_idx)
				if (DEBUG_LEVEL >= 2):
					fp = open(outfile, 'a')
					fp.write('\n ===>> (cluster pair): ' + str(x)  + ' and ' + str(CURRENT_CLUST_IDX_LIST[dest_reach_mat_idx]) + \
						'  Connected by reln R1')      
					fp.close()
		
		"""
		for A->B connection
		if B->E exists
		then establish A->E
		"""
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
			if (Reachability_Graph_Mat[src_reach_mat_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
				Connect_ClusterPair(Reachability_Graph_Mat, src_reach_mat_idx, \
					CURRENT_CLUST_IDX_LIST.index(x), RELATION_R1, src_taxa_clust_idx, x)
				if (DEBUG_LEVEL >= 2):
					fp = open(outfile, 'a')
					fp.write('\n ===>> (cluster pair): ' + str(CURRENT_CLUST_IDX_LIST[src_reach_mat_idx])  + ' and ' + str(x) + \
						'  Connected by reln R1')      
					fp.close()

		"""
		for A->B connection
		if D->A and B->E exists
		then establish D->E  
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R2):
			for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
				if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 0):  
					Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
						CURRENT_CLUST_IDX_LIST.index(y), RELATION_R1, x, y)  
					if (DEBUG_LEVEL >= 2):
						fp = open(outfile, 'a')
						fp.write('\n ===>> (cluster pair): ' + str(x)  + ' and ' + str(y) + \
							'  Connected by reln R1')      
						fp.close()

		"""
		for A->B connection
		if D><A exists
		then establish D><B
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R4):
			if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][dest_reach_mat_idx] == 0):
				Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
					dest_reach_mat_idx, RELATION_R4, x, dest_taxa_clust_idx)
				if (DEBUG_LEVEL >= 2):
					fp = open(outfile, 'a')
					fp.write('\n ===>> (cluster pair): ' + str(x)  + ' and ' + str(CURRENT_CLUST_IDX_LIST[dest_reach_mat_idx]) + \
						'  Connected by reln R4')      
					fp.close()
		
		"""
		for A->B connection
		if D><A exists
		then for all B->E
		establish D><E
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetClustRelnList(RELATION_R4):
			for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetClustRelnList(RELATION_R1):
				if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 0):
					Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
						CURRENT_CLUST_IDX_LIST.index(y), RELATION_R4, x, y)  
					if (DEBUG_LEVEL >= 2):
						fp = open(outfile, 'a')
						fp.write('\n ===>> (cluster pair): ' + str(x)  + ' and ' + str(y) + \
							'  Connected by reln R4')      
						fp.close()
		
	else:
		"""
		the target relation is R4
		"""
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
		for any couplet (x,y) where x in src_clust_out_neighb, and y in dest_clust_out_neighb,
		(x, y) will be related via the relation R4
		"""
		for x in src_clust_out_neighb:
			for y in dest_clust_out_neighb:
				if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 0):
					Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
						CURRENT_CLUST_IDX_LIST.index(y), RELATION_R4, x, y)  
					if (DEBUG_LEVEL >= 2):
						fp = open(outfile, 'a')
						fp.write('\n ===>> (cluster pair): ' + str(x)  + ' and ' + str(y) + \
							'  Connected by reln R4')      
						fp.close()

	return

#-----------------------------------------------------
""" 
this function merges two clusters 
basically the cluster src_clust_idx will be merged to the dest_clust_idx
so all the entries concerning src_clust_idx will now point to the dest_clust_idx
"""
def Merge_Clusters(Reachability_Graph_Mat, dest_clust_idx, src_clust_idx, \
			dest_clust_reach_mat_idx, src_clust_reach_mat_idx, outfile):
	
	""" 
	first update the reachability matrix entries 
	originally all the Reachability_Graph_Mat entries concerning src_clust_reach_mat_idx 
	will now point to the dest_clust_reach_mat_idx 
	also update the dest cluster out edge and in edge lists 
	"""
	"""
	Important - Note -
	In copying the out / in / no edge information
	we can overwrite the no edge with a definite out / in edge
	"""
	for x in Cluster_Info_Dict[src_clust_idx]._GetClustRelnList(RELATION_R1):
		Cluster_Info_Dict[x]._RemoveRelnInstance(RELATION_R2, src_clust_idx)
		if (Reachability_Graph_Mat[dest_clust_reach_mat_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
			Connect_ClusterPair(Reachability_Graph_Mat, dest_clust_reach_mat_idx, \
				CURRENT_CLUST_IDX_LIST.index(x), RELATION_R1, dest_clust_idx, x)  

	for x in Cluster_Info_Dict[src_clust_idx]._GetClustRelnList(RELATION_R2):
		Cluster_Info_Dict[x]._RemoveRelnInstance(RELATION_R1, src_clust_idx)
		if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][dest_clust_reach_mat_idx] == 0):
			Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
				dest_clust_reach_mat_idx, RELATION_R1, x, dest_clust_idx)  

	for x in Cluster_Info_Dict[src_clust_idx]._GetClustRelnList(RELATION_R4):
		Cluster_Info_Dict[x]._RemoveRelnInstance(RELATION_R4, src_clust_idx)
		if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][dest_clust_reach_mat_idx] == 0) and \
			(Reachability_Graph_Mat[dest_clust_reach_mat_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
			Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
				dest_clust_reach_mat_idx, RELATION_R4, x, dest_clust_idx)      
	
	# then adjust the taxa in the other cluster
	for tax in Cluster_Info_Dict[src_clust_idx]._GetSpeciesList():
		Append_Cluster_Taxa_Label(dest_clust_idx, tax)
  
	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')
		fp.write('\n ===>> Successfully merged the cluster ' + str(src_clust_idx) + '  to the cluster ' + str(dest_clust_idx))      
		fp.close()
  
	return
    
#-----------------------------------------------------
""" 
this function updates the reachability graph 
on the basis of input edge type between input 2 taxa 
"""
def AdjustReachGraph(Reachability_Graph_Mat, nodeA_clust_idx, nodeB_clust_idx, reln_type, outfile):    
  
	nodeA_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(nodeA_clust_idx)
	nodeB_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(nodeB_clust_idx)
		
	if (reln_type == RELATION_R3):
		"""
		keep the minimum cluster index intact
		merge two clusters
		"""
		if (nodeA_clust_idx > nodeB_clust_idx):
			Merge_Clusters(Reachability_Graph_Mat, nodeB_clust_idx, \
				nodeA_clust_idx, nodeB_reach_mat_idx, nodeA_reach_mat_idx, outfile)
			"""
			delete the index of nodeA_clust_idx from the reachability matrix
			"""
			Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeA_reach_mat_idx), axis=0)	# delete the row
			Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeA_reach_mat_idx), axis=1)	# delete the column
			"""
			delete the entry of nodeA_clust_idx from the CURRENT_CLUST_IDX_LIST
			"""
			CURRENT_CLUST_IDX_LIST.remove(nodeA_clust_idx)
			"""
			also remove the cluster key from the dictionary
			"""
			Cluster_Info_Dict.pop(nodeA_clust_idx, None)          
		else:
			Merge_Clusters(Reachability_Graph_Mat, nodeA_clust_idx, nodeB_clust_idx, \
				nodeA_reach_mat_idx, nodeB_reach_mat_idx, outfile)    
			"""
			delete the index of nodeB_clust_idx from the reachability matrix
			"""
			Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeB_reach_mat_idx), axis=0)	# delete the row
			Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeB_reach_mat_idx), axis=1)	# delete the column
			"""
			delete the entry of nodeB_clust_idx from the CURRENT_CLUST_IDX_LIST
			"""
			CURRENT_CLUST_IDX_LIST.remove(nodeB_clust_idx)
			"""
			also remove the cluster key from the dictionary
			"""
			Cluster_Info_Dict.pop(nodeB_clust_idx, None)    
		#-----------------------------
	elif (reln_type == RELATION_R1):
		"""
		connect the pair of clusters, along with updating the reachability matrix
		"""
		Connect_ClusterPair(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, \
			RELATION_R1, nodeA_clust_idx, nodeB_clust_idx)
		"""
		now perform the transitive closure on the derived reachability matrix  
		"""
		TransClosUpd(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, RELATION_R1, outfile)
	elif (reln_type == RELATION_R2):
		"""
		connect the pair of clusters, along with updating the reachability matrix
		"""
		Connect_ClusterPair(Reachability_Graph_Mat, nodeB_reach_mat_idx, nodeA_reach_mat_idx, \
			RELATION_R1, nodeB_clust_idx, nodeA_clust_idx)
		"""
		now perform the transitive closure on the derived reachability matrix  
		"""
		TransClosUpd(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, RELATION_R2, outfile)    
	else:	#reln_type == RELATION_R4:
		"""
		connect the pair of clusters, along with updating the reachability matrix
		"""
		Connect_ClusterPair(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, \
			RELATION_R4, nodeA_clust_idx, nodeB_clust_idx)
		"""
		now perform the transitive closure on the derived reachability matrix  
		"""
		TransClosUpd(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, RELATION_R4, outfile)    
			
	return Reachability_Graph_Mat

#-----------------------------------------------------
""" 
this function merges two clusters, by simply copying the contents in the cluster "src_clust_idx"
to the cluster "dest_clust_idx"
"""
def Copy_Cluster_Content(dest_clust_idx, src_clust_idx, outfile):
	"""
	copy the set of taxa belonging in the cluster "src_clust_idx"
	to the cluster "dest_clust_idx"
	"""
	for tax in Cluster_Info_Dict[src_clust_idx]._GetSpeciesList():
		Append_Cluster_Taxa_Label(dest_clust_idx, tax)
  
	"""
	delete the entry of src_clust_idx from the CURRENT_CLUST_IDX_LIST
	"""
	CURRENT_CLUST_IDX_LIST.remove(src_clust_idx)
	"""
	also remove the "src_clust_idx" cluster key from the dictionary
	"""
	Cluster_Info_Dict.pop(src_clust_idx, None)          
  
	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')
		fp.write('\n ===>> Successfully merged the cluster ' + str(src_clust_idx) + '  to the cluster ' + str(dest_clust_idx))      
		fp.close()

	return
