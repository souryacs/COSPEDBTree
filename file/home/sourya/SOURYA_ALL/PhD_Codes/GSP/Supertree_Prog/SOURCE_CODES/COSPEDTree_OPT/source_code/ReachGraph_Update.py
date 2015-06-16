#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

##-----------------------------------------------------
""" this function adds an edge between a pair of clusters (of taxa) 
it also updates the entries of reachability matrix """
def Connect_ClusterPair(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, edge_type, nodeA_clust_idx, nodeB_clust_idx):
  if (edge_type == DIRECTED_OUT_EDGE):
    # adjust the clusters
    Cluster_Info_Dict[nodeA_clust_idx]._AddOutEdge(nodeB_clust_idx)
    Cluster_Info_Dict[nodeB_clust_idx]._AddInEdge(nodeA_clust_idx)
    # update the reachability matrix
    Reachability_Graph_Mat[nodeA_reach_mat_idx][nodeB_reach_mat_idx] = 1
  elif (edge_type == NO_EDGE):
    # adjust the clusters
    Cluster_Info_Dict[nodeA_clust_idx]._AddNoEdge(nodeB_clust_idx)
    Cluster_Info_Dict[nodeB_clust_idx]._AddNoEdge(nodeA_clust_idx)    
    # update the reachability matrix
    Reachability_Graph_Mat[nodeA_reach_mat_idx][nodeB_reach_mat_idx] = 2
    Reachability_Graph_Mat[nodeB_reach_mat_idx][nodeA_reach_mat_idx] = 2
        
##-----------------------------------------------------
""" this function updates the transitive closure of the cluster of nodes
on inclusion ogf a new edge between a pair of clusters """
def TransClosUpd(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, taxaA_label, taxaB_label, edge_type):
  if (edge_type == DIRECTED_OUT_EDGE) or (edge_type == NO_EDGE):
    src_reach_mat_idx = nodeA_reach_mat_idx
    dest_reach_mat_idx = nodeB_reach_mat_idx
    src_taxa_label = taxaA_label
    dest_taxa_label = taxaB_label
  elif (edge_type == DIRECTED_IN_EDGE):
    src_reach_mat_idx = nodeB_reach_mat_idx
    dest_reach_mat_idx = nodeA_reach_mat_idx
    src_taxa_label = taxaB_label
    dest_taxa_label = taxaA_label    
  else:
    return
    
  src_taxa_clust_idx = CURRENT_CLUST_IDX_LIST[src_reach_mat_idx]
  dest_taxa_clust_idx = CURRENT_CLUST_IDX_LIST[dest_reach_mat_idx]
    
  if (edge_type == DIRECTED_OUT_EDGE) or (edge_type == DIRECTED_IN_EDGE):
    # for A->B connection
    # if D->A exists
    # then establish D->B
    for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetInEdgeList():
      if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][dest_reach_mat_idx] == 0):
	Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), dest_reach_mat_idx, DIRECTED_OUT_EDGE, x, dest_taxa_clust_idx)
    
    # for A->B connection
    # if B->E exists
    # then establish A->E
    for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
      if (Reachability_Graph_Mat[src_reach_mat_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
	Connect_ClusterPair(Reachability_Graph_Mat, src_reach_mat_idx, CURRENT_CLUST_IDX_LIST.index(x), DIRECTED_OUT_EDGE, src_taxa_clust_idx, x)

    # for A->B connection
    # if D->A and B->E exists
    # then establish D->E  
    for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetInEdgeList():
      for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
	if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 0):  
	  Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), CURRENT_CLUST_IDX_LIST.index(y), DIRECTED_OUT_EDGE, x, y)  

    # for A->B connection
    # if D><A exists
    # then establish D><B
    for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetNoEdgeList():
      if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][dest_reach_mat_idx] == 0):
	Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), dest_reach_mat_idx, NO_EDGE, x, dest_taxa_clust_idx)
	      
    # for A->B connection
    # if D><A exists
    # then for all B->E
    # establish D><E
    for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetNoEdgeList():
      for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
	if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 0):
	  Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), CURRENT_CLUST_IDX_LIST.index(y), NO_EDGE, x, y)  
    
  else:
    # construct the out neighborhood of src_cluster
    # it will contain the cluster itself and all the other clusters connected via out edges from this cluster
    src_clust_out_neighb = []
    src_clust_out_neighb.append(src_taxa_clust_idx)
    src_clust_out_neighb.extend(Cluster_Info_Dict[src_taxa_clust_idx]._GetOutEdgeList())
    # construct the out neighborhood of dest cluster
    # it will contain the cluster itself and all the other clusters connected via out edges from this cluster    
    dest_clust_out_neighb = []
    dest_clust_out_neighb.append(dest_taxa_clust_idx)
    dest_clust_out_neighb.extend(Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList())
    
    for x in src_clust_out_neighb:
      for y in dest_clust_out_neighb:
	if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 0):
	  Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), CURRENT_CLUST_IDX_LIST.index(y), NO_EDGE, x, y)  

##-----------------------------------------------------
""" this function merges two clusters 
basically the cluster src_clust_idx will be merged to the dest_clust_idx
so all the entries concerning src_clust_idx will now point to the dest_clust_idx"""
def Merge_Clusters(Reachability_Graph_Mat, dest_taxa_label, src_taxa_label, dest_clust_idx, src_clust_idx, \
		    dest_clust_reach_mat_idx, src_clust_reach_mat_idx):
  """ first update the reachability matrix entries 
  originally all the Reachability_Graph_Mat entries concerning src_clust_reach_mat_idx 
  will now point to the dest_clust_reach_mat_idx 
  also update the dest cluster out edge and in edge lists """
  
  for x in Cluster_Info_Dict[src_clust_idx]._GetOutEdgeList():
    Cluster_Info_Dict[x]._RemoveInEdge(src_clust_idx)
    #Cluster_Info_Dict[src_clust_idx]._RemoveOutEdge(x)
    if (Reachability_Graph_Mat[dest_clust_reach_mat_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
      Connect_ClusterPair(Reachability_Graph_Mat, dest_clust_reach_mat_idx, CURRENT_CLUST_IDX_LIST.index(x), DIRECTED_OUT_EDGE, dest_clust_idx, x)  

  for x in Cluster_Info_Dict[src_clust_idx]._GetInEdgeList():
    Cluster_Info_Dict[x]._RemoveOutEdge(src_clust_idx)
    #Cluster_Info_Dict[src_clust_idx]._RemoveInEdge(x)
    if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][dest_clust_reach_mat_idx] == 0):
      Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), dest_clust_reach_mat_idx, DIRECTED_OUT_EDGE, x, dest_clust_idx)  

  for x in Cluster_Info_Dict[src_clust_idx]._GetNoEdgeList():
    Cluster_Info_Dict[x]._RemoveNoEdge(src_clust_idx)
    #Cluster_Info_Dict[src_clust_idx]._RemoveNoEdge(x)
    if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][dest_clust_reach_mat_idx] == 0):
      Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), dest_clust_reach_mat_idx, NO_EDGE, x, dest_clust_idx)      
          
  # then adjust the taxa in the other cluster
  for tax in Cluster_Info_Dict[src_clust_idx]._GetSpeciesList():
    Append_Cluster_Taxa_Label(dest_clust_idx, tax)
    
##-----------------------------------------------------
""" this function updates the reachability graph 
on the basis of input edge type between input 2 taxa """
def AdjustReachGraph(Reachability_Graph_Mat, taxaA_label, taxaB_label, edge_type):    
  
  nodeA_clust_idx = Taxa_Info_Dict[taxaA_label]._Get_Taxa_Part_Clust_Idx()
  nodeB_clust_idx = Taxa_Info_Dict[taxaB_label]._Get_Taxa_Part_Clust_Idx()

  # perform cluster merging operation (content shifting plus transitive closure)
  nodeA_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(nodeA_clust_idx)
  nodeB_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(nodeB_clust_idx)
  
  # establish the final edge connection between these two taxa
  #Connect_TaxaPair(taxaA_label, taxaB_label, edge_type)
  
  # update tge reachability matrix information
  if (edge_type == BI_DIRECTED_EDGE):      
    #-----------------------------
    # keep the minimum cluster index intact
    # merge two clusters
    if (nodeA_clust_idx > nodeB_clust_idx):
      Merge_Clusters(Reachability_Graph_Mat, taxaB_label, taxaA_label, nodeB_clust_idx, nodeA_clust_idx, nodeB_reach_mat_idx, nodeA_reach_mat_idx)
      # delete the index of nodeA_clust_idx from the reachability matrix
      Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeA_reach_mat_idx), axis=0)	# delete the row
      Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeA_reach_mat_idx), axis=1)	# delete the column
      # delete the entry of nodeA_clust_idx from the CURRENT_CLUST_IDX_LIST
      CURRENT_CLUST_IDX_LIST.remove(nodeA_clust_idx)
      # also remove the cluster key from the dictionary
      Cluster_Info_Dict.pop(nodeA_clust_idx, None)          
    else:
      Merge_Clusters(Reachability_Graph_Mat, taxaA_label, taxaB_label, nodeA_clust_idx, nodeB_clust_idx, nodeA_reach_mat_idx, nodeB_reach_mat_idx)    
      # delete the index of nodeB_clust_idx from the reachability matrix
      Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeB_reach_mat_idx), axis=0)	# delete the row
      Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeB_reach_mat_idx), axis=1)	# delete the column
      # delete the entry of nodeB_clust_idx from the CURRENT_CLUST_IDX_LIST
      CURRENT_CLUST_IDX_LIST.remove(nodeB_clust_idx)
      # also remove the cluster key from the dictionary
      Cluster_Info_Dict.pop(nodeB_clust_idx, None)    
    #-----------------------------
  elif (edge_type == DIRECTED_OUT_EDGE):
    # connect the pair of clusters, along with updating the reachability matrix
    Connect_ClusterPair(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, \
			DIRECTED_OUT_EDGE, nodeA_clust_idx, nodeB_clust_idx)
    # now perform the transitive closure on the derived reachability matrix  
    TransClosUpd(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, taxaA_label, taxaB_label, DIRECTED_OUT_EDGE)
  elif (edge_type == DIRECTED_IN_EDGE):
    # connect the pair of clusters, along with updating the reachability matrix
    Connect_ClusterPair(Reachability_Graph_Mat, nodeB_reach_mat_idx, nodeA_reach_mat_idx, \
			DIRECTED_OUT_EDGE, nodeB_clust_idx, nodeA_clust_idx)
    # now perform the transitive closure on the derived reachability matrix  
    TransClosUpd(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, taxaA_label, taxaB_label, DIRECTED_IN_EDGE)    
  else:	#edge_type == NO_EDGE:
    Connect_ClusterPair(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, \
			NO_EDGE, nodeA_clust_idx, nodeB_clust_idx)
    # now perform the transitive closure on the derived reachability matrix  
    TransClosUpd(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, taxaA_label, taxaB_label, NO_EDGE)    
      
  return Reachability_Graph_Mat
