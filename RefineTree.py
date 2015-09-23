#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#--------------------------------------------------------
"""
this function is a shortcut to obtain the normalized expression 
used in the agglomerative clustering proposed in this code
as various methods are experimented, corresponding various forms of 
agglomerative clustering is tried
"""
def ObtainNormalizedVal(num, denom1, denom2):
	if ((denom1 + denom2) > 0):
		return (num * 1.0) / (denom1 + denom2)
	else:
		return 0

##---------------------------------------------
""" 
function to print the matrix content
N is the matrix dimension
"""
def PrintMatrixContent(N, TaxaList, inp_data, inp_str, textfile):
	fp = open(textfile, 'a')
	fp.write('\n printing contents of ' + str(inp_str) + ' ---- ')
	for i in range(N):
		fp.write('\n ' + str(i) + '--' + str(TaxaList[i]) + '--->>')
		for j in range(i+1):
			fp.write(' ' + str(inp_data[i][j]))
	fp.close()

#-------------------------------------------
""" 
this function takes inputs of two taxa clusters
for individual couplets (x,y) within these taxa clusters
it accumulates their average extra lineage information (already computed from the input trees and couplet based processing)
finally the average of this accumulated sum is returned as output
this average value is used for the proposed binary refinement
"""
def Compute_Avg_Clust_Dist_Score(clust_spec_list1, clust_spec_list2):
	clust_avg_dist = 0
	no_of_couplets = 0
	for x in clust_spec_list1:
		for y in clust_spec_list2:
			key1 = (x, y)
			key2 = (y, x)
			if (key1 in TaxaPair_Reln_Dict):
				target_key = key1
			elif (key2 in TaxaPair_Reln_Dict):
				target_key = key2
			else:
				continue

			# if there exists valid key then only this part will be processed
			clust_avg_dist = clust_avg_dist + TaxaPair_Reln_Dict[target_key]._GetAvgTreeXL()
			
			# increase the number of couplets
			no_of_couplets = no_of_couplets + 1      

	# error condition check
	if (no_of_couplets == 0):
		#print '--- couplet score : 0 --- '
		return 0, 1	# signifies that this couplet is not supported in the input trees
	return (clust_avg_dist * 1.0) / no_of_couplets, 0

#-------------------------------------------
""" 
this function resolves the multifurcation of a node within Curr_tree
the taxa clusters underlying the subtree rooted at the multifurcation node 
are given in the argument clust_species_list
"""
def ResolveMultifurcation(Curr_tree, clust_species_list, no_of_input_clusters, NJ_RULE_USED, Output_Text_File):

	# total number of clusters
	no_of_clust = no_of_input_clusters

	# allocate two 2D square matrices of no_of_clust dimension
	# first one stores the pairwise taxa cluster based score
	Clust_DistMat1 = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.float)

	# second one stores the normalized XL score employed in the proposed agglomerative clustering
	# computed using the used distance metric of all cluster pairs 
	Norm_Clust_DistMat1 = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.float)

	# compute the pairwise XL score of clusters
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n taxa list pair analysis:  i - ' + str(i) + '  j - ' + str(j) + \
					' 1st list : ' + str(clust_species_list[i]) + ' 2nd list - ' + str(clust_species_list[j]))
				fp.close()

			# this function computes the average XL score of the cluster pair
			clust_avg_dist, error_condition = Compute_Avg_Clust_Dist_Score(clust_species_list[i], clust_species_list[j], DIST_METRIC)    
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write(' clust_avg_dist: ' + str(clust_avg_dist) + '  error: ' + str(error_condition))
				fp.close()
			Clust_DistMat1[j][i] = Clust_DistMat1[i][j] = clust_avg_dist
			
	#--------------------------------------------------------
	# loop to execute the agglomerative clustering
	# for binary refinement
	while(no_of_clust > 2):      
		
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n iteration start --- number of clusters: ' + str(no_of_clust))
			fp.write('\n clust_species_list : ' + str(clust_species_list))
			fp.close()
			if 0:
				PrintMatrixContent(no_of_clust, clust_species_list, Clust_DistMat1, 'Clust_DistMat1', Output_Text_File)
		
		# allocate a new array for storing d(C,:) where C is a cluster and : denotes all other clusers
		# individual element of this list stores the distance sum for one particular cluster C 
		# as we allow 0 couplet score for even unsupported couplets, the sum is not affected
		sum_list = []
		for i in range(no_of_clust):
			t = 0
			for j in range(no_of_clust):
				t = t + Clust_DistMat1[i][j]
			sum_list.append(t)

		# normalize the cluster pair distance (with respect to the obtained XL count)
		# by dividing it with the sum of all cluster pair distance values
		for i in range(no_of_clust - 1):
			for j in range(i+1, no_of_clust):
				if (NJ_RULE_USED == AGGLO_CLUST):
					# this is agglomerative clustering
					Norm_Clust_DistMat1[i][j] = ObtainNormalizedVal(Clust_DistMat1[i][j], sum_list[i], sum_list[j])
					Norm_Clust_DistMat1[j][i] = Norm_Clust_DistMat1[i][j]
				else:
					# this is NJ based rule
					ri = sum_list[i] / (no_of_clust - 2)
					rj = sum_list[j] / (no_of_clust - 2)
					Norm_Clust_DistMat1[i][j] = (Clust_DistMat1[i][j] - ri - rj)
					Norm_Clust_DistMat1[j][i] = Norm_Clust_DistMat1[i][j]	  

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n printing contents of sum_list --- ' + str(sum_list))
			fp.close()
			PrintMatrixContent(no_of_clust, clust_species_list, Norm_Clust_DistMat1, 'Norm_Clust_DistMat1', Output_Text_File)
		#----------------------------------------------------------
		# now find the minimum of the distance matrix
		min_val = Norm_Clust_DistMat1[0][1]
		min_idx_i = 0
		min_idx_j = 1
		for i in range(no_of_clust - 1):
			for j in range(i+1, no_of_clust):
				if (i == j):
					continue
				if (Norm_Clust_DistMat1[i][j] < min_val):
					min_val = Norm_Clust_DistMat1[i][j]
					min_idx_i = i
					min_idx_j = j
				elif (Norm_Clust_DistMat1[i][j] == min_val):
					# here we prioritize the cluster pair having minimum number of species
					if (len(clust_species_list[i]) + len(clust_species_list[j])) < (len(clust_species_list[min_idx_i]) + len(clust_species_list[min_idx_j])):
						min_idx_i = i
						min_idx_j = j
			
		#----------------------------------------------------------
		# note down the taxa list in these two indices of the clust_species_list
		taxa_list = []
		for x in clust_species_list[min_idx_i]:
			taxa_list.append(x)
		for x in clust_species_list[min_idx_j]:
			taxa_list.append(x)

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n min_idx_i ' + str(min_idx_i) + ' min_idx_j : ' + str(min_idx_j))
			fp.write('\n min_idx_i species list ' + str(clust_species_list[min_idx_i]))
			fp.write('\n min_idx_j species list ' + str(clust_species_list[min_idx_j]))
			fp.write('\n complete species list ' + str(taxa_list))
			fp.close()
						
		#---------------------------------------------------------      
		# for individual clusters, we check if the cluster contains one or more species
						
		# case 1 - both the clusters have > 1 species
		# and the clusters are represented by an internal node which is the MRCA of the constituent species set
		if (len(clust_species_list[min_idx_i]) > 1) and (len(clust_species_list[min_idx_j]) > 1):
			first_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_i])
			second_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_j])
			all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
						
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n label of first_cluster_mrca_node: ' + str(Node_Label(first_cluster_mrca_node)))      
				fp.write('\n label of second_cluster_mrca_node: ' + str(Node_Label(second_cluster_mrca_node)))
				fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
				fp.close()
			
			# create new internal node 
			newnode = Node()  
			# its parent node will be the previous MRCA node of all the taxa in two clusters
			all_taxa_mrca_node.add_child(newnode)
			newnode.parent_node = all_taxa_mrca_node            
			all_taxa_mrca_node.remove_child(first_cluster_mrca_node)
			first_cluster_mrca_node.parent_node = None
			all_taxa_mrca_node.remove_child(second_cluster_mrca_node)
			second_cluster_mrca_node.parent_node = None
			# add these individual clusters' MRCA node as its children
			newnode.add_child(first_cluster_mrca_node)
			first_cluster_mrca_node.parent_node = newnode
			newnode.add_child(second_cluster_mrca_node)
			second_cluster_mrca_node.parent_node = newnode      
			# update splits of the resulting tree
			Curr_tree.update_splits(delete_outdegree_one=False)
			
		# case 2 and 3 - one cluster has at least 2 species, while other is a leaf
		elif (len(clust_species_list[min_idx_i]) == 1) and (len(clust_species_list[min_idx_j]) > 1):
			first_cluster_leaf_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_i][0])
			second_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_j])
			all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
			
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n first cluster is a leaf - its label: ' + str(Node_Label(first_cluster_leaf_node)))      
				fp.write('\n label of second_cluster_mrca_node: ' + str(Node_Label(second_cluster_mrca_node)))
				fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
				fp.close()
			
			# create new internal node 
			newnode = Node()        
			# its parent node will be the previous MRCA node of all the taxa in two clusters
			all_taxa_mrca_node.add_child(newnode)
			newnode.parent_node = all_taxa_mrca_node
			all_taxa_mrca_node.remove_child(first_cluster_leaf_node)
			first_cluster_leaf_node.parent_node = None
			all_taxa_mrca_node.remove_child(second_cluster_mrca_node)
			second_cluster_mrca_node.parent_node = None
			# add these individual clusters' MRCA node as its children
			newnode.add_child(first_cluster_leaf_node)
			first_cluster_leaf_node.parent_node = newnode
			newnode.add_child(second_cluster_mrca_node)
			second_cluster_mrca_node.parent_node = newnode
			# update splits of the resulting tree
			Curr_tree.update_splits(delete_outdegree_one=False)
			
		elif (len(clust_species_list[min_idx_i]) > 1) and (len(clust_species_list[min_idx_j]) == 1):
			first_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_i])
			second_cluster_leaf_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_j][0])
			all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
			
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n label of first_cluster_mrca_node: ' + str(Node_Label(first_cluster_mrca_node)))      
				fp.write('\n second cluster is a leaf - its label: ' + str(Node_Label(second_cluster_leaf_node)))
				fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
				fp.close()
			
			# create new internal node 
			newnode = Node()        
			# its parent node will be the previous MRCA node of all the taxa in two clusters
			all_taxa_mrca_node.add_child(newnode)
			newnode.parent_node = all_taxa_mrca_node
			all_taxa_mrca_node.remove_child(first_cluster_mrca_node)
			first_cluster_mrca_node.parent_node = None
			all_taxa_mrca_node.remove_child(second_cluster_leaf_node)
			second_cluster_leaf_node.parent_node = None
			# add these individual clusters' MRCA node as its children
			newnode.add_child(first_cluster_mrca_node)
			first_cluster_mrca_node.parent_node = newnode
			newnode.add_child(second_cluster_leaf_node)
			second_cluster_leaf_node.parent_node = newnode
			# update splits of the resulting tree
			Curr_tree.update_splits(delete_outdegree_one=False)
			
		# case 4 - when both child clusters are leaf nodes 
		else:
			first_cluster_leaf_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_i][0])
			second_cluster_leaf_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_j][0])      
			all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
			
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n first cluster is a leaf - its label: ' + str(Node_Label(first_cluster_leaf_node)))      
				fp.write('\n second cluster is a leaf - its label: ' + str(Node_Label(second_cluster_leaf_node)))
				fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
				fp.close()
			
			# create new internal node 
			newnode = Node()        
			# its parent node will be the previous MRCA node of all the taxa in two clusters
			all_taxa_mrca_node.add_child(newnode)
			newnode.parent_node = all_taxa_mrca_node
			all_taxa_mrca_node.remove_child(first_cluster_leaf_node)
			first_cluster_leaf_node.parent_node = None
			all_taxa_mrca_node.remove_child(second_cluster_leaf_node)
			second_cluster_leaf_node.parent_node = None
			# add these individual clusters' MRCA node as its children
			newnode.add_child(first_cluster_leaf_node)
			first_cluster_leaf_node.parent_node = newnode
			newnode.add_child(second_cluster_leaf_node)
			second_cluster_leaf_node.parent_node = newnode
			# update splits of the resulting tree
			Curr_tree.update_splits(delete_outdegree_one=False)
					
		#---------------------------------------------------------          
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n label of newnode: ' + str(Node_Label(newnode)))
			fp.write('\n label of all taxa mrca node (recomputed): ' + str(Node_Label(Curr_tree.mrca(taxon_labels=taxa_list))))                  
			#fp.write('\n before inserting row col, Clust_DistMat1 dimension: ' + str(Clust_DistMat1.size))
			fp.close()
		
		# adjust the Clust_DistMat1 by inserting one new row and column corresponding to the new cluster
		# and then deleting the information of earlier two clusters
		# first append one row
		Clust_DistMat1 = numpy.vstack((Clust_DistMat1, numpy.zeros((1, no_of_clust), dtype=numpy.float)))
		# then append one column
		Clust_DistMat1 = numpy.hstack((Clust_DistMat1, numpy.zeros((no_of_clust + 1, 1), dtype=numpy.float)))

		# then reshape the matrix with respect to the added dimension
		# note: the dimension is increased by (n+1) x (n+1) from the earlier n x n
		# where n is the number of clusters
		Clust_DistMat1 = numpy.reshape(Clust_DistMat1, ((no_of_clust + 1), (no_of_clust + 1)), order='C')
				
		# now fill the elements of the new added row and column of the Clust_DistMat1
		for k in range(no_of_clust):
			Clust_DistMat1[k][no_of_clust] = max(Clust_DistMat1[k][min_idx_i], Clust_DistMat1[k][min_idx_j], Clust_DistMat1[min_idx_i][min_idx_j])
			# symmetric property maintain
			Clust_DistMat1[no_of_clust][k] = Clust_DistMat1[k][no_of_clust]
			
		# now remove the rows and columns corresponding to min_idx_i and min_idx_j
		Clust_DistMat1 = numpy.delete(Clust_DistMat1, (min_idx_i), axis=0)	# delete the row
		Clust_DistMat1 = numpy.delete(Clust_DistMat1, (min_idx_i), axis=1)	# delete the column
		Clust_DistMat1 = numpy.delete(Clust_DistMat1, (min_idx_j - 1), axis=0)	# delete the row
		Clust_DistMat1 = numpy.delete(Clust_DistMat1, (min_idx_j - 1), axis=1)	# delete the column
						
		# clear Norm_Clust_DistMat1
		Norm_Clust_DistMat1.fill(0)
		Norm_Clust_DistMat1 = numpy.delete(Norm_Clust_DistMat1, (min_idx_i), axis=0)	# delete the row
		Norm_Clust_DistMat1 = numpy.delete(Norm_Clust_DistMat1, (min_idx_i), axis=1)	# delete the column
		
		# remove individual clusters' taxa information from the clust_species_list
		# and add taxa_list as a new element
		clust_species_list.pop(min_idx_i)
		clust_species_list.pop(min_idx_j - 1)
		clust_species_list.append(taxa_list)    
		
		# decrement the number of clusters considered
		no_of_clust = no_of_clust - 1

	return
      
#-------------------------------------------
# this function refines input supertree such that the supertree becomes binary
# this is required for proper benchmarking with existing binary tree construction methods on 
# ILS sorting
def Refine_Supertree_Binary_Form(Curr_tree, NJ_RULE_USED, Output_Text_File):
	# we traverse input tree internal nodes in postorder fashion
	# and list the child nodes of it
	# if the no of children > 2 then it is a case of multifurcation
	# for resolving
	for curr_node in Curr_tree.postorder_internal_node_iter():
		curr_node_children = curr_node.child_nodes()
		if (len(curr_node_children) > 2):
			# create a list which will contain the species list lying under 
			# individual child nodes of rhe current node
			clust_species_list = []
			# examine individual nodes of the current node's children list
			for x in curr_node_children:
				if (x.is_leaf() == True):
					subl = []
					subl.append(x.taxon.label)
				else:
					subl = GetTaxaUnderInternalNode(x)
				clust_species_list.append(subl)
			
			# call the resolving routine
			ResolveMultifurcation(Curr_tree, clust_species_list, len(curr_node_children), NJ_RULE_USED, Output_Text_File)
