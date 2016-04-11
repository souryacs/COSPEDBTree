#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

##-----------------------------------------------------
#"""
#this function deletes the support score queue entry
#@param: 
#src_clust_idx, dest_clust_idx: indices of clusters containing single taxon each
#inp_reln: input relation type
#"""
#def Del_Support_Score_Entry(src_clust_idx, dest_clust_idx, inp_reln, outfile):
	#if (DEBUG_LEVEL >= 2):
		#fp = open(outfile, 'a')
	
	#t1 = COMPLETE_INPUT_TAXA_LIST[src_clust_idx]
	#t2 = COMPLETE_INPUT_TAXA_LIST[dest_clust_idx]
	#key1 = (t1, t2)
	#key2 = (t2, t1)
	
	#if key1 in TaxaPair_Reln_Dict:
		#queue_entry = [t1, t2, inp_reln, TaxaPair_Reln_Dict[key1]._GetEdgeCost_ConnReln(inp_reln)]
		#if queue_entry in Cost_List_Taxa_Pair_Multi_Reln:
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n ==> Deleting the support score queue entry: ' + str(queue_entry))
			#Cost_List_Taxa_Pair_Multi_Reln.remove(queue_entry)
	
	#if key2 in TaxaPair_Reln_Dict:
		#queue_entry = [t2, t1, Complementary_Reln(inp_reln), TaxaPair_Reln_Dict[key2]._GetEdgeCost_ConnReln(Complementary_Reln(inp_reln))]
		#if queue_entry in Cost_List_Taxa_Pair_Multi_Reln:
			#if (DEBUG_LEVEL >= 2):
				#fp.write('\n ==> Deleting the support score queue entry: ' + str(queue_entry))
			#Cost_List_Taxa_Pair_Multi_Reln.remove(queue_entry)
	
	#if (DEBUG_LEVEL >= 2):
		#fp.close()
	
	#return

#-----------------------------------------------------
#"""
#this function processes individual clusters (where one cluster corresponds to a single taxa)
#and processes its edge (relation) lists
#it removes redunant (unnecessary) relations from those lists 
#"""
#def DelRedundantReln(Output_Text_File):
	#"""
	#for a cluster cx, if there exists cx---cy and also cx->cy / cx<-cy / cx=cy edges
	#then keep only cx--cy (R4 relation) edge, and remove all other edges
	#"""
	#for cl in Cluster_Info_Dict:
		#"""
		#list of cluster indices cy such that both cx---cy and cx->cy are present
		#"""
		#common_idx = [v for v in Cluster_Info_Dict[cl]._GetNoEdgeList() if v in Cluster_Info_Dict[cl]._GetOutEdgeList()]
		#if (len(common_idx) > 0):
			#for x in common_idx:
				#"""
				#remove cx->cy information
				#"""
				#Cluster_Info_Dict[cl]._RemoveOutEdge(x)
				#Cluster_Info_Dict[x]._RemoveInEdge(cl)
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n ==> (cluster level) As both ' + str(cl) + '---' + str(x) + '  and  ' \
						#+ str(cl) + '-->' + str(x) + ' is true - deleting ' + str(cl) + '-->' + str(x))
					#fp.close()
				#"""
				#Delete the corresponding support score entry
				#"""
				#Del_Support_Score_Entry(cl, x, RELATION_R1, Output_Text_File)
				
				###----------------------------------------------------------
				### add - sourya
				##"""
				##search clusters cz such that cy-->cz / cy==cz exists and also either cx-->cz or cx==cz exists
				##in such a case, remove cx-->cz or cx==cz connection
				##"""
				##x_ol = set(Cluster_Info_Dict[x]._GetOutEdgeList())
				##x_el = set(Cluster_Info_Dict[x]._GetEqEdgeList())
				##cl_ol = set(Cluster_Info_Dict[cl]._GetOutEdgeList())
				##cl_el = set(Cluster_Info_Dict[cl]._GetEqEdgeList())
				##cz_list = list((x_ol.union(x_el)).intersection(cl_ol.union(cl_el)))
				##if (len(cz_list) > 0):
					##for cz in cz_list:
						##if cz in Cluster_Info_Dict[cl]._GetOutEdgeList():
							##"""
							##remove cx->cz information
							##"""
							##Cluster_Info_Dict[cl]._RemoveOutEdge(cz)
							##Cluster_Info_Dict[cz]._RemoveInEdge(cl)
							##if (DEBUG_LEVEL >= 2):
								##fp = open(Output_Text_File, 'a')
								##fp.write('\n ==> After such deletion - ' + str(x) + '--> / == ' + str(cz) + '  and  ' \
									##+ str(cl) + '-->' + str(cz) + ' - deleting ' + str(cl) + '-->' + str(cz))
								##fp.close()

							##"""
							##Delete the corresponding support score entry
							##"""
							##Del_Support_Score_Entry(cl, cz, RELATION_R1, Output_Text_File)
							
						##if cz in Cluster_Info_Dict[cl]._GetEqEdgeList():
							##"""
							##remove cx==cz information
							##"""
							##Cluster_Info_Dict[cl]._RemoveEqEdge(cz)
							##Cluster_Info_Dict[cz]._RemoveEqEdge(cl)
							##if (DEBUG_LEVEL >= 2):
								##fp = open(Output_Text_File, 'a')
								##fp.write('\n ==> After such deletion - ' + str(x) + '--> / == ' + str(cz) + '  and  ' \
									##+ str(cl) + '==' + str(cz) + ' - deleting ' + str(cl) + '==' + str(cz))
								##fp.close()
			
							##"""
							##Delete the corresponding support score entry
							##"""
							##Del_Support_Score_Entry(cl, cz, RELATION_R3, Output_Text_File)
				
				### end add - sourya
				###-------------------------------------------------------------------
				
		#"""
		#list of cluster indices cy such that both cx---cy and cx<-cy are present
		#"""
		#common_idx = [v for v in Cluster_Info_Dict[cl]._GetNoEdgeList() if v in Cluster_Info_Dict[cl]._GetInEdgeList()]
		#if (len(common_idx) > 0):
			#for x in common_idx:
				#"""
				#remove cx<-cy information
				#"""
				#Cluster_Info_Dict[cl]._RemoveInEdge(x)
				#Cluster_Info_Dict[x]._RemoveOutEdge(cl)
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n ==> (cluster level) As both ' + str(cl) + '---' + str(x) + '  and  ' \
						#+ str(cl) + '<--' + str(x) + ' is true - deleting ' + str(cl) + '<--' + str(x))
					#fp.close()
				#"""
				#Delete the corresponding support score entry
				#"""
				#Del_Support_Score_Entry(cl, x, RELATION_R2, Output_Text_File)
		
				###-----------------------------------------------------------
				### add - sourya
				
				##"""
				##search clusters cz such that cy<--cz / cy==cz exists and also either cx<--cz or cx==cz exists
				##in such a case, remove cx<--cz or cx==cz connection
				##"""
				##x_il = set(Cluster_Info_Dict[x]._GetInEdgeList())
				##x_el = set(Cluster_Info_Dict[x]._GetEqEdgeList())
				##cl_il = set(Cluster_Info_Dict[cl]._GetInEdgeList())
				##cl_el = set(Cluster_Info_Dict[cl]._GetEqEdgeList())
				##cz_list = list((x_il.union(x_el)).intersection(cl_il.union(cl_el)))
				##if (len(cz_list) > 0):
					##for cz in cz_list:
						##if cz in Cluster_Info_Dict[cl]._GetInEdgeList():
							##"""
							##remove cx<--cz information
							##"""
							##Cluster_Info_Dict[cl]._RemoveInEdge(cz)
							##Cluster_Info_Dict[cz]._RemoveOutEdge(cl)
							##if (DEBUG_LEVEL >= 2):
								##fp = open(Output_Text_File, 'a')
								##fp.write('\n ==> After such deletion - ' + str(x) + '<-- / == ' + str(cz) + '  and  ' \
									##+ str(cl) + '<--' + str(cz) + ' - deleting ' + str(cl) + '<--' + str(cz))
								##fp.close()

							##"""
							##Delete the corresponding support score entry
							##"""
							##Del_Support_Score_Entry(cl, cz, RELATION_R2, Output_Text_File)
							
						##if cz in Cluster_Info_Dict[cl]._GetEqEdgeList():
							##"""
							##remove cx==cz information
							##"""
							##Cluster_Info_Dict[cl]._RemoveEqEdge(cz)
							##Cluster_Info_Dict[cz]._RemoveEqEdge(cl)
							##if (DEBUG_LEVEL >= 2):
								##fp = open(Output_Text_File, 'a')
								##fp.write('\n ==> After such deletion - ' + str(x) + '<-- / == ' + str(cz) + '  and  ' \
									##+ str(cl) + '==' + str(cz) + ' - deleting ' + str(cl) + '==' + str(cz))
								##fp.close()
			
							##"""
							##Delete the corresponding support score entry
							##"""
							##Del_Support_Score_Entry(cl, cz, RELATION_R3, Output_Text_File)
				
				### end add - sourya
				###-----------------------------------------------------------
		
		#"""
		#list of cluster indices cy such that both cx---cy and cx==cy are present
		#"""
		#common_idx = [v for v in Cluster_Info_Dict[cl]._GetNoEdgeList() if v in Cluster_Info_Dict[cl]._GetEqEdgeList()]
		#if (len(common_idx) > 0):
			#for x in common_idx:
				#"""
				#remove cx==cy information
				#"""
				#Cluster_Info_Dict[cl]._RemoveEqEdge(x)
				#Cluster_Info_Dict[x]._RemoveEqEdge(cl)
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n ==> (cluster level) As both ' + str(cl) + '---' + str(x) + '  and  ' \
						#+ str(cl) + '==' + str(x) + ' is true - deleting ' + str(cl) + '==' + str(x))
					#fp.close()
				#"""
				#Delete the corresponding support score entry
				#"""
				#Del_Support_Score_Entry(cl, x, RELATION_R3, Output_Text_File)
		
				###---------------------------------------------------------
				### add - sourya
				
				##"""
				##search clusters cz such that cy-->cz / cy==cz exists and also either cx-->cz or cx==cz exists
				##in such a case, remove cx-->cz or cx==cz connection
				##"""
				##x_ol = set(Cluster_Info_Dict[x]._GetOutEdgeList())
				##x_el = set(Cluster_Info_Dict[x]._GetEqEdgeList())
				##cl_ol = set(Cluster_Info_Dict[cl]._GetOutEdgeList())
				##cl_el = set(Cluster_Info_Dict[cl]._GetEqEdgeList())
				##cz_list = list((x_ol.union(x_el)).intersection(cl_ol.union(cl_el)))
				##if (len(cz_list) > 0):
					##for cz in cz_list:
						##if cz in Cluster_Info_Dict[cl]._GetOutEdgeList():
							##"""
							##remove cx->cz information
							##"""
							##Cluster_Info_Dict[cl]._RemoveOutEdge(cz)
							##Cluster_Info_Dict[cz]._RemoveInEdge(cl)
							##if (DEBUG_LEVEL >= 2):
								##fp = open(Output_Text_File, 'a')
								##fp.write('\n ==> After such deletion - ' + str(x) + '--> / == ' + str(cz) + '  and  ' \
									##+ str(cl) + '-->' + str(cz) + ' - deleting ' + str(cl) + '-->' + str(cz))
								##fp.close()

							##"""
							##Delete the corresponding support score entry
							##"""
							##Del_Support_Score_Entry(cl, cz, RELATION_R1, Output_Text_File)
							
						##if cz in Cluster_Info_Dict[cl]._GetEqEdgeList():
							##"""
							##remove cx==cz information
							##"""
							##Cluster_Info_Dict[cl]._RemoveEqEdge(cz)
							##Cluster_Info_Dict[cz]._RemoveEqEdge(cl)
							##if (DEBUG_LEVEL >= 2):
								##fp = open(Output_Text_File, 'a')
								##fp.write('\n ==> After such deletion - ' + str(x) + '--> / == ' + str(cz) + '  and  ' \
									##+ str(cl) + '==' + str(cz) + ' - deleting ' + str(cl) + '==' + str(cz))
								##fp.close()
			
							##"""
							##Delete the corresponding support score entry
							##"""
							##Del_Support_Score_Entry(cl, cz, RELATION_R3, Output_Text_File)
		
				##"""
				##search clusters cz such that cy<--cz / cy==cz exists and also either cx<--cz or cx==cz exists
				##in such a case, remove cx<--cz or cx==cz connection
				##"""
				##x_il = set(Cluster_Info_Dict[x]._GetInEdgeList())
				##x_el = set(Cluster_Info_Dict[x]._GetEqEdgeList())
				##cl_il = set(Cluster_Info_Dict[cl]._GetInEdgeList())
				##cl_el = set(Cluster_Info_Dict[cl]._GetEqEdgeList())
				##cz_list = list((x_il.union(x_el)).intersection(cl_il.union(cl_el)))
				##if (len(cz_list) > 0):
					##for cz in cz_list:
						##if cz in Cluster_Info_Dict[cl]._GetInEdgeList():
							##"""
							##remove cx<--cz information
							##"""
							##Cluster_Info_Dict[cl]._RemoveInEdge(cz)
							##Cluster_Info_Dict[cz]._RemoveOutEdge(cl)
							##if (DEBUG_LEVEL >= 2):
								##fp = open(Output_Text_File, 'a')
								##fp.write('\n ==> After such deletion - ' + str(x) + '<-- / == ' + str(cz) + '  and  ' \
									##+ str(cl) + '<--' + str(cz) + ' - deleting ' + str(cl) + '<--' + str(cz))
								##fp.close()

							##"""
							##Delete the corresponding support score entry
							##"""
							##Del_Support_Score_Entry(cl, cz, RELATION_R2, Output_Text_File)
							
						##if cz in Cluster_Info_Dict[cl]._GetEqEdgeList():
							##"""
							##remove cx==cz information
							##"""
							##Cluster_Info_Dict[cl]._RemoveEqEdge(cz)
							##Cluster_Info_Dict[cz]._RemoveEqEdge(cl)
							##if (DEBUG_LEVEL >= 2):
								##fp = open(Output_Text_File, 'a')
								##fp.write('\n ==> After such deletion - ' + str(x) + '<-- / == ' + str(cz) + '  and  ' \
									##+ str(cl) + '==' + str(cz) + ' - deleting ' + str(cl) + '==' + str(cz))
								##fp.close()
			
							##"""
							##Delete the corresponding support score entry
							##"""
							##Del_Support_Score_Entry(cl, cz, RELATION_R3, Output_Text_File)
				
				### end add - sourya
				###--------------------------------------------------------------
		
	#"""
	#for a cluster cx, if there exists cx==cy and also cx->cy / cx<-cy edges
	#then delete cx==cy (R3 relation) edge
	#"""
	#for cl in Cluster_Info_Dict:
		#"""
		#list of cluster indices cy such that both cx==cy and cx->cy are present
		#"""
		#common_idx = [v for v in Cluster_Info_Dict[cl]._GetEqEdgeList() if v in Cluster_Info_Dict[cl]._GetOutEdgeList()]
		#if (len(common_idx) > 0):
			#for x in common_idx:
				#"""
				#remove cx==cy information
				#"""
				#Cluster_Info_Dict[cl]._RemoveEqEdge(x)
				#Cluster_Info_Dict[x]._RemoveEqEdge(cl)
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n ==> (cluster level) As both ' + str(cl) + '-->' + str(x) + '  and  ' \
						#+ str(cl) + '==' + str(x) + ' is true - deleting ' + str(cl) + '==' + str(x))
					#fp.close()
				#"""
				#Delete the corresponding support score entry
				#"""
				#Del_Support_Score_Entry(cl, x, RELATION_R3, Output_Text_File)
				
				###---------------------------------------------------------
				### add - sourya
				
				##"""
				##search clusters cz such that cy==cz exists and also cx==cz exists
				##in such a case, remove cx==cz connection
				##"""
				##x_el = set(Cluster_Info_Dict[x]._GetEqEdgeList())
				##cl_el = set(Cluster_Info_Dict[cl]._GetEqEdgeList())
				##cz_list = list(x_el.intersection(cl_el))
				##if (len(cz_list) > 0):
					##for cz in cz_list:
						##"""
						##remove cx==cz information
						##"""
						##Cluster_Info_Dict[cl]._RemoveEqEdge(cz)
						##Cluster_Info_Dict[cz]._RemoveEqEdge(cl)
						##if (DEBUG_LEVEL >= 2):
							##fp = open(Output_Text_File, 'a')
							##fp.write('\n ==> After such deletion - ' + str(x) + '--> / == ' + str(cz) + '  and  ' \
								##+ str(cl) + '==' + str(cz) + ' - deleting ' + str(cl) + '==' + str(cz))
							##fp.close()
		
						##"""
						##Delete the corresponding support score entry
						##"""
						##Del_Support_Score_Entry(cl, cz, RELATION_R3, Output_Text_File)
		
				### end add - sourya
				###--------------------------------------------------------------
				
		#"""
		#list of cluster indices cy such that both cx==cy and cx<-cy are present
		#"""
		#common_idx = [v for v in Cluster_Info_Dict[cl]._GetEqEdgeList() if v in Cluster_Info_Dict[cl]._GetInEdgeList()]
		#if (len(common_idx) > 0):
			#for x in common_idx:
				#"""
				#remove cx==cy information
				#"""
				#Cluster_Info_Dict[cl]._RemoveEqEdge(x)
				#Cluster_Info_Dict[x]._RemoveEqEdge(cl)
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n ==> (cluster level) As both ' + str(cl) + '<--' + str(x) + '  and  ' \
						#+ str(cl) + '==' + str(x) + ' is true - deleting ' + str(cl) + '==' + str(x))
					#fp.close()
				#"""
				#Delete the corresponding support score entry
				#"""
				#Del_Support_Score_Entry(cl, x, RELATION_R3, Output_Text_File)

				###---------------------------------------------------------
				### add - sourya
				
				##"""
				##search clusters cz such that cy==cz exists and also cx==cz exists
				##in such a case, remove cx==cz connection
				##"""
				##x_el = set(Cluster_Info_Dict[x]._GetEqEdgeList())
				##cl_el = set(Cluster_Info_Dict[cl]._GetEqEdgeList())
				##cz_list = list(x_el.intersection(cl_el))
				##if (len(cz_list) > 0):
					##for cz in cz_list:
						##"""
						##remove cx==cz information
						##"""
						##Cluster_Info_Dict[cl]._RemoveEqEdge(cz)
						##Cluster_Info_Dict[cz]._RemoveEqEdge(cl)
						##if (DEBUG_LEVEL >= 2):
							##fp = open(Output_Text_File, 'a')
							##fp.write('\n ==> After such deletion - ' + str(x) + '--> / == ' + str(cz) + '  and  ' \
								##+ str(cl) + '==' + str(cz) + ' - deleting ' + str(cl) + '==' + str(cz))
							##fp.close()
		
						##"""
						##Delete the corresponding support score entry
						##"""
						##Del_Support_Score_Entry(cl, cz, RELATION_R3, Output_Text_File)
		
				### end add - sourya
				###--------------------------------------------------------------

	#return

##----------------------------------------
#"""
#this function is called after transitive closure of the cluster connectivity
#few relations from the support score queue 
#"""
#def DelRelnMore(Output_Text_File):
	#"""
	#if there exists A->C, A--B and B->C
	#"""




##------------------------------------------------------------
#"""
#this function applies transitive closure operation for the cluster based connectivity
#"""
#def Transitive_Closure_Cluster_Connectivity(Output_Text_File):
	#no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
	#"""
	#transitive closure operation
	#"""

	#"""
	#case A - if A->B and B->C then A->C, provided A--C does not exist
	#case B - if A->B and B==C then A->C, provided A--C does not exist
	#"""
	#for i in range(no_of_clusters):
		#clust_i = CURRENT_CLUST_IDX_LIST[i]
		#for j in range(no_of_clusters):
			#clust_j = CURRENT_CLUST_IDX_LIST[j]
			#if (clust_i == clust_j):
				#continue
			#if clust_j in Cluster_Info_Dict[clust_i]._GetOutEdgeList():
				#"""
				#A->B holds with respect to input gene trees
				#"""
				#for k in range(no_of_clusters):
					#clust_k = CURRENT_CLUST_IDX_LIST[k]
					#if (clust_i == clust_k) or (clust_j == clust_k):
						#continue
					#if (clust_k in Cluster_Info_Dict[clust_j]._GetOutEdgeList()): 
						#"""
						#B->C holds with respect to input gene trees
						#"""
						#if (DEBUG_LEVEL >= 2):
							#fp = open(Output_Text_File, 'a')
							#fp.write('\n ==> (cluster level) Both ' + str(clust_i) + '-->' + str(clust_j) + '  and  ' \
								#+ str(clust_j) + '-->' + str(clust_k) + ' is true ')
							#fp.close()

						#""" 
						#we check whether A--C does not exist
						#in such a case we establish A->C relation
						#"""
						#if clust_k not in Cluster_Info_Dict[clust_i]._GetNoEdgeList():
							#Cluster_Info_Dict[clust_i]._AddOutEdge(clust_k)
							#Cluster_Info_Dict[clust_k]._AddInEdge(clust_i)
							#if (DEBUG_LEVEL >= 2):
								#fp = open(Output_Text_File, 'a')
								#fp.write('\n ==> As ' + str(clust_i) + '--' + str(clust_k) + ' does not exist - we set  ' \
									#+ str(clust_i) + '-->' + str(clust_k))
								#fp.close()
							
						#"""
						#otherwise, if A--C exists before, we remove A->B and also delete the entry from the support score queue
						#"""
						#if clust_k in Cluster_Info_Dict[clust_i]._GetNoEdgeList():
							#Cluster_Info_Dict[clust_i]._RemoveOutEdge(clust_j)
							#Cluster_Info_Dict[clust_j]._RemoveInEdge(clust_i)
							#if (DEBUG_LEVEL >= 2):
								#fp = open(Output_Text_File, 'a')
								#fp.write('\n ==> As ' + str(clust_i) + '--' + str(clust_k) + ' exists previously - we delete  ' \
									#+ str(clust_i) + '-->' + str(clust_j))
								#fp.close()
							#"""
							#Delete the corresponding support score entry
							#"""
							#Del_Support_Score_Entry(clust_i, clust_j, RELATION_R1, Output_Text_File)
							
					#if (clust_k in Cluster_Info_Dict[clust_j]._GetEqEdgeList()):
						#"""
						#B==C holds with respect to input gene trees
						#"""
						#if (DEBUG_LEVEL >= 2):
							#fp = open(Output_Text_File, 'a')
							#fp.write('\n ==> (cluster level) Both ' + str(clust_i) + '-->' + str(clust_j) + '  and  ' \
								#+ str(clust_j) + '==' + str(clust_k) + ' is true ')
							#fp.close()
						
						#""" 
						#we check whether A--C does not exist
						#in such a case we establish A->C relation
						#"""
						#if clust_k not in Cluster_Info_Dict[clust_i]._GetNoEdgeList():
							#Cluster_Info_Dict[clust_i]._AddOutEdge(clust_k)
							#Cluster_Info_Dict[clust_k]._AddInEdge(clust_i)
							#if (DEBUG_LEVEL >= 2):
								#fp = open(Output_Text_File, 'a')
								#fp.write('\n ==> As ' + str(clust_i) + '--' + str(clust_k) + ' does not exist - we set  ' \
									#+ str(clust_i) + '-->' + str(clust_k))
								#fp.close()
							
						#"""
						#otherwise, if A--C exists before, we remove A->B and also delete the entry from the support score queue
						#"""
						#if clust_k in Cluster_Info_Dict[clust_i]._GetNoEdgeList():
							#Cluster_Info_Dict[clust_i]._RemoveOutEdge(clust_j)
							#Cluster_Info_Dict[clust_j]._RemoveInEdge(clust_i)
							#if (DEBUG_LEVEL >= 2):
								#fp = open(Output_Text_File, 'a')
								#fp.write('\n ==> As ' + str(clust_i) + '--' + str(clust_k) + ' exists previously - we delete  ' \
									#+ str(clust_i) + '-->' + str(clust_j))
								#fp.close()
							#"""
							#Delete the corresponding support score entry
							#"""
							#Del_Support_Score_Entry(clust_i, clust_j, RELATION_R1, Output_Text_File)
						
				
	#"""
	#case C - if A==B and B->C then A->C
	#"""
	#for i in range(no_of_clusters):
		#clust_i = CURRENT_CLUST_IDX_LIST[i]
		#for j in range(no_of_clusters):
			#clust_j = CURRENT_CLUST_IDX_LIST[j]
			#if (clust_i == clust_j):
				#continue
			#if clust_j in Cluster_Info_Dict[clust_i]._GetEqEdgeList():
				#"""
				#A==B holds with respect to input gene trees
				#"""
				#for k in range(no_of_clusters):
					#clust_k = CURRENT_CLUST_IDX_LIST[k]
					#if (clust_i == clust_k) or (clust_j == clust_k):
						#continue
					#if clust_k in Cluster_Info_Dict[clust_j]._GetOutEdgeList():
						#"""
						#B->C holds with respect to input gene trees
						#"""
						#if (DEBUG_LEVEL >= 2):
							#fp = open(Output_Text_File, 'a')
							#fp.write('\n ==> (cluster level) Both ' + str(clust_i) + '==' + str(clust_j) + '  and  ' \
								#+ str(clust_j) + '-->' + str(clust_k) + ' is true ')
							#fp.close()
						
						#""" 
						#we check whether A--C does not exist
						#in such a case we establish A->C relation
						#"""
						#if clust_k not in Cluster_Info_Dict[clust_i]._GetNoEdgeList():
							#Cluster_Info_Dict[clust_i]._AddOutEdge(clust_k)
							#Cluster_Info_Dict[clust_k]._AddInEdge(clust_i)
							#if (DEBUG_LEVEL >= 2):
								#fp = open(Output_Text_File, 'a')
								#fp.write('\n ==> As ' + str(clust_i) + '--' + str(clust_k) + ' does not exist - we set  ' \
									#+ str(clust_i) + '-->' + str(clust_k))
								#fp.close()
							
						#"""
						#otherwise, if A--C exists before, we remove A==B and also delete the entry from the support score queue
						#"""
						#if clust_k in Cluster_Info_Dict[clust_i]._GetNoEdgeList():
							#Cluster_Info_Dict[clust_i]._RemoveEqEdge(clust_j)
							#Cluster_Info_Dict[clust_j]._RemoveEqEdge(clust_i)
							#if (DEBUG_LEVEL >= 2):
								#fp = open(Output_Text_File, 'a')
								#fp.write('\n ==> As ' + str(clust_i) + '--' + str(clust_k) + ' exists previously - we delete  ' \
									#+ str(clust_i) + '==' + str(clust_j))
								#fp.close()
							#"""
							#Delete the corresponding support score entry
							#"""
							#Del_Support_Score_Entry(clust_i, clust_j, RELATION_R3, Output_Text_File)
						

	###----------------------------------------------
	##"""
	##transitive reduction - but first create a copy of the cluster dictionary obtained so far
	##we search the clusters using this copies, but modify the original cluster dictionary
	##"""
	##Cluster_Info_Dict_Copy = Cluster_Info_Dict.copy()
	
	##"""
	##now apply transitive reduction to the cluster contents
	##"""

	##"""
	##case A - if A->B, B->C, and A->C exists, remove A->C
	##"""
	##for clust_i in Cluster_Info_Dict_Copy:
		##clust_i_out_edge_list = Cluster_Info_Dict_Copy[clust_i]._GetOutEdgeList()
		##if (len(clust_i_out_edge_list) > 0):
			##for clust_j in clust_i_out_edge_list:
				##clust_j_out_edge_list = Cluster_Info_Dict_Copy[clust_j]._GetOutEdgeList()
				##clust_k_list = [for v in clust_j_out_edge_list if v in clust_i_out_edge_list]
				##if (len(clust_k_list) > 0):
					##for clust_k in clust_k_list:
						##Cluster_Info_Dict[clust_i]._RemoveOutEdge(clust_k)
						##Cluster_Info_Dict[clust_k]._RemoveInEdge(clust_i)
				
			
	##"""
	##case B - if A==B, B->C, and A->C exists, remove A->C
	##"""
	##for clust_i in Cluster_Info_Dict_Copy:
		##clust_i_out_edge_list = Cluster_Info_Dict_Copy[clust_i]._GetOutEdgeList()
		##if (len(clust_i_out_edge_list) > 0):
			##for clust_j in clust_i_out_edge_list:
				##clust_j_out_edge_list = Cluster_Info_Dict_Copy[clust_j]._GetOutEdgeList()
				##clust_k_list = [for v in clust_j_out_edge_list if v in clust_i_out_edge_list]
				##if (len(clust_k_list) > 0):
					##for clust_k in clust_k_list:
						##Cluster_Info_Dict[clust_i]._RemoveOutEdge(clust_k)
						##Cluster_Info_Dict[clust_k]._RemoveInEdge(clust_i)
		

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

