#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#-------------------------------------------
"""
this function fills the distance matrix using normalized excess gene count
for a particular taxa cluster, it uses one representative taxon of that taxa cluster
"""
def Fill_DistMat_SingleEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE):
	"""
	check and explore each pair of taxa clusters
	"""
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			"""
			clust_species_list[i] and clust_species_list[j]
			contain two taxa list of one or more elements
			"""
			if (DIST_MAT_TYPE == 1) or (DIST_MAT_TYPE == 2):
				entry = FindAvgDistanceMeasure(clust_species_list[i], clust_species_list[j], DIST_MAT_TYPE, 1)
			elif (DIST_MAT_TYPE == 3):
				entry = FindAvgDistanceMeasure(clust_species_list[i], clust_species_list[j], DIST_MAT_TYPE, 1, 0)	#2)
			elif (DIST_MAT_TYPE == 4):
				entry = FindAvgDistanceMeasure(clust_species_list[i], clust_species_list[j], DIST_MAT_TYPE, 1, 0)	#2)
			else:
				entry = FindAvgDistanceMeasure(clust_species_list[i], clust_species_list[j], DIST_MAT_TYPE, 1, 0)
			"""
			copy the entry in the distance matrix
			"""
			DistMat[i][j] = entry
			DistMat[j][i] = DistMat[i][j]

	return
	
#-------------------------------------------
"""
this function fills the distance matrix using normalized excess gene count
for a particular taxa cluster, it uses aerage information of that taxa cluster
@param type_of_output: if 0, computes the average of XL measures
													1, returns the minimum of XL measures
													2, returns the maximum of XL measures
"""
def Fill_DistMat_AvgEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE):
	"""
	check and explore each pair of taxa clusters
	"""
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			"""
			clust_species_list[i] and clust_species_list[j]
			contain two taxa list of one or more elements
			"""
			if (DIST_MAT_TYPE == 1) or (DIST_MAT_TYPE == 2):
				entry = FindAvgDistanceMeasure(clust_species_list[i], clust_species_list[j], DIST_MAT_TYPE, 2)
			elif (DIST_MAT_TYPE == 3):
				entry = FindAvgDistanceMeasure(clust_species_list[i], clust_species_list[j], DIST_MAT_TYPE, 2, 0)	#2)
			elif (DIST_MAT_TYPE == 4):
				entry = FindAvgDistanceMeasure(clust_species_list[i], clust_species_list[j], DIST_MAT_TYPE, 2, 0)	#2)
			else:
				entry = FindAvgDistanceMeasure(clust_species_list[i], clust_species_list[j], DIST_MAT_TYPE, 2, 2)
			"""
			copy the entry in the distance matrix
			"""
			DistMat[i][j] = entry
			DistMat[j][i] = DistMat[i][j]
	
	return

#-------------------------------------------
"""
this function finds a single minimum from the input matrix
"""
def Find_Unique_Min_Valid(DistMat, Norm_DistMat, Valid_Mat, no_of_clust, clust_species_list):

	flag = False
	
	# traverse through the matrix elements
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			if (i == j):
				continue
			if (Valid_Mat[i][j] > 0):	# nonzero entry  means valid distance matrix entry
				if (flag == False):
					min_val = Norm_DistMat[i][j]
					min_idx_i = i
					min_idx_j = j
					flag = True
					
				elif (FlEq(Norm_DistMat[i][j], min_val) == True): 	#(Norm_DistMat[i][j] == min_val):
					if (DistMat[i][j] < DistMat[min_idx_i][min_idx_j]):
						min_idx_i = i
						min_idx_j = j
				
				elif (Norm_DistMat[i][j] < min_val):
					min_val = Norm_DistMat[i][j]
					min_idx_i = i
					min_idx_j = j
	
	
	if (flag == True):
		return min_idx_i, min_idx_j
	else:
		return -1, -1

#-------------------------------------------
"""
this function finds a single minimum from the input matrix
"""
def Find_Unique_Min(DistMat, Norm_DistMat, no_of_clust, clust_species_list):
	
	min_val = Norm_DistMat[0][1]
	min_idx_i = 0
	min_idx_j = 1
	
	# traverse through the matrix elements
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			if (i == j) or ((i == 0) and (j == 1)):
				continue
			elif (FlEq(Norm_DistMat[i][j], min_val) == True): 	#(Norm_DistMat[i][j] == min_val):
				if (DistMat[i][j] < DistMat[min_idx_i][min_idx_j]):
					min_idx_i = i
					min_idx_j = j
			elif (Norm_DistMat[i][j] < min_val):
				min_val = Norm_DistMat[i][j]
				min_idx_i = i
				min_idx_j = j
	
	return min_idx_i, min_idx_j
	
#--------------------------------------------------------
# this function is a shortcut to obtain the normalized expression 
# used in the agglomerative clustering proposed in this code
# as various methods are experimented, corresponding various forms of 
# agglomerative clustering is tried
#--------------------------------------------------------
def ObtainNormalizedVal(num, denom1, denom2):
  if ((denom1 + denom2) > 0):
    return (num * 1.0) / (denom1 + denom2)
  else:
    return 0

#---------------------------------------------
""" 
function to print the matrix content
@parameters: 
	N = matrix dimension
	TaxaList = input set of taxa (clusters)
	inp_data = matrix data
	inp_str = a string, depicting the matrix name
	textfile = output file 
	format_mat  can be one of the three values:
		1: Lower triangular format 
		2: Upper triangular format
		0: full matrix print
"""
def PrintMatrixContent(N, TaxaList, inp_data, inp_str, textfile, format_mat=1):
	fp = open(textfile, 'a')
	
	fp.write('\n\n\n printing contents of ' + str(inp_str) + ' ---- ')
	if (format_mat == 1):	# lower triangular format printing
		for i in range(N):
			#fp.write('\n ' + str(i) + '--' + str(TaxaList[i]) + '--->>')
			fp.write('\n ' + str(i) + '------>>')
			if 1:	#(i > 0):
				for j in range(i+1):	#i+1
					fp.write(' ' + str(inp_data[i][j]))
	elif (format_mat == 0):	# full matrix printing
		for i in range(N):
			#fp.write('\n ' + str(i) + '--' + str(TaxaList[i]) + '--->>')
			fp.write('\n ' + str(i) + '------>>')
			if 1:
				for j in range(N):	#i+1
					fp.write(' ' + str(inp_data[i][j]))
	else:	# upper triangular format printing
		for i in range(N):
			#fp.write('\n ' + str(i) + '--' + str(TaxaList[i]) + '--->>')
			fp.write('\n ' + str(i) + '------>>')
			if 1:	#(i > 0):
				for j in range(i, N):
					fp.write(' ' + str(inp_data[i][j]))
					
	fp.close()

#-------------------------------------------
"""
this function processes input distance matrix in every iteration
and finds the pair of indices satisfying minimum distance criterion 
used in NJ based algorithm
"""
def Get_NJ_Based_Min_Pair_Idx(DistMat, Norm_DistMat, no_of_clust, clust_species_list, NJ_RULE_USED, Output_Text_File):
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n\n\n\n **** iteration start --- number of clusters: ' + str(no_of_clust))
		fp.write('\n\n *** clust_species_list : **** ')
		for i in range(len(clust_species_list)):
			fp.write('\n Index: ' + str(i) + ' Taxa list: ' + str(clust_species_list[i]))
		fp.close()
		PrintMatrixContent(no_of_clust, clust_species_list, DistMat, 'DistMat', Output_Text_File)
	
	"""
	allocate one new square matrix which will contain the 
	normalized matrix elements (w.r.t the sum of sum of rows and columns)
	"""
	sum_list = []
	for i in range(no_of_clust):
		t = 0
		for j in range(no_of_clust):
			t = t + DistMat[i][j]
		sum_list.append(t)

	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			if (NJ_RULE_USED == AGGLO_CLUST):
				# we normalize the extra lineage based score
				# by the sum of extra lineages for all other taxa from the taxa indices i and j
				# modified - sourya
				#Norm_DistMat[i][j] = ObtainNormalizedVal(DistMat[i][j], sum_list[i], sum_list[j])
				# add - sourya
				Norm_DistMat[i][j] = ObtainNormalizedVal(DistMat[i][j], sum_list[i] - DistMat[i][j], sum_list[j] - DistMat[i][j])
				# end add - sourya
				Norm_DistMat[j][i] = Norm_DistMat[i][j]
			else:
				#ri = (sum_list[i] * 1.0) / (no_of_clust - 2)
				#rj = (sum_list[j] * 1.0) / (no_of_clust - 2)
				#Norm_DistMat[i][j] = (DistMat[i][j] - ri - rj)
				Norm_DistMat[i][j] = (no_of_clust - 2) * DistMat[i][j] - sum_list[i] - sum_list[j]
				Norm_DistMat[j][i] = Norm_DistMat[i][j]
				
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n printing contents of sum_list --- ' + str(sum_list))
		fp.close()
		PrintMatrixContent(no_of_clust, clust_species_list, Norm_DistMat, 'Norm_DistMat', Output_Text_File)
		
	"""
	now we have to find the minimum among these elements 
	present in the matrix Norm_DistMat
	"""
	min_idx_i, min_idx_j = Find_Unique_Min(DistMat, Norm_DistMat, no_of_clust, clust_species_list)

	return min_idx_i, min_idx_j

#-------------------------------------------
"""
this function processes input distance matrix in every iteration
and finds the pair of indices satisfying minimum distance criterion 
used in NJ based algorithm
"""
def Get_NJ_Based_Min_Pair_Idx_Valid(DistMat, Norm_DistMat, Valid_Mat, no_of_clust, clust_species_list, NJ_RULE_USED, Output_Text_File):
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n\n\n\n **** iteration start --- number of clusters: ' + str(no_of_clust))
		fp.write('\n\n *** clust_species_list : **** ')
		for i in range(len(clust_species_list)):
			fp.write('\n Index: ' + str(i) + ' Taxa list: ' + str(clust_species_list[i]))
		fp.close()
		PrintMatrixContent(no_of_clust, clust_species_list, DistMat, 'DistMat', Output_Text_File)
	
	"""
	allocate one new square matrix which will contain the 
	normalized matrix elements (w.r.t the sum of sum of rows and columns)
	"""
	sum_list = []
	for i in range(no_of_clust):
		t = 0
		for j in range(no_of_clust):
			t = t + DistMat[i][j]
		sum_list.append(t)

	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			if (NJ_RULE_USED == AGGLO_CLUST):
				# we normalize the extra lineage based score
				# by the sum of extra lineages for all other taxa from the taxa indices i and j
				# modified - sourya
				#Norm_DistMat[i][j] = ObtainNormalizedVal(DistMat[i][j], sum_list[i], sum_list[j])
				# add - sourya
				Norm_DistMat[i][j] = ObtainNormalizedVal(DistMat[i][j], sum_list[i] - DistMat[i][j], sum_list[j] - DistMat[i][j])
				# end add - sourya
				Norm_DistMat[j][i] = Norm_DistMat[i][j]
			else:
				#ri = (sum_list[i] * 1.0) / (no_of_clust - 2)
				#rj = (sum_list[j] * 1.0) / (no_of_clust - 2)
				#Norm_DistMat[i][j] = (DistMat[i][j] - ri - rj)
				Norm_DistMat[i][j] = (no_of_clust - 2) * DistMat[i][j] - sum_list[i] - sum_list[j]
				Norm_DistMat[j][i] = Norm_DistMat[i][j]
				
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n printing contents of sum_list --- ' + str(sum_list))
		fp.close()
		PrintMatrixContent(no_of_clust, clust_species_list, Norm_DistMat, 'Norm_DistMat', Output_Text_File)
		
	"""
	now we have to find the minimum among these elements 
	present in the matrix Norm_DistMat
	"""
	min_idx_i, min_idx_j = Find_Unique_Min_Valid(DistMat, Norm_DistMat, Valid_Mat, no_of_clust, clust_species_list)

	return min_idx_i, min_idx_j

#-------------------------------------------
"""
checks whether a taxa cluster specified by the input index is a leaf
"""
def IsLeafCluster(clust_species_list, idx):
	if (len(clust_species_list[idx]) == 1):
		return True
	return False

#-------------------------------------------
"""
this function has following parameters:
1) first_cluster_mrca_node: root of 1st subtree 
2) second_cluster_mrca_node: root of 2nd subtree 
3) all_taxa_mrca_node: root of all these trees
4) Curr_tree: Tree containing all these subtrees

It creates one new internal node as a child of all_taxa_mrca_node
and places above mentioned subtrees as its children
"""
def MergeSubtrees(Curr_tree, first_cluster_mrca_node, second_cluster_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File):
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n label of first_cluster_mrca_node: ' + str(Node_Label(first_cluster_mrca_node)))      
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

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n label of newnode: ' + str(Node_Label(newnode)))
		fp.write('\n label of all taxa mrca node (recomputed): ' + str(Node_Label(Curr_tree.mrca(taxon_labels=taxa_list))))  
		fp.close()

	return Curr_tree

#-------------------------------------------
"""
this function merges a pair of clusters whose indices are pointed by the min_idx_i and min_idx_j entries
this is part of the proposed agglomerative clustering
taxa_list is the union of these two clusters (species contents)
"""
def Merge_Cluster_Pair(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File):
	isleaf_clust1 = IsLeafCluster(clust_species_list, min_idx_i)
	isleaf_clust2 = IsLeafCluster(clust_species_list, min_idx_j)
	
	if (isleaf_clust1):
		first_cluster_mrca_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_i][0])
	else:
		first_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_i])
	
	if (isleaf_clust2):
		second_cluster_mrca_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_j][0])
	else:
		second_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_j])
	
	all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	
	Curr_tree = MergeSubtrees(Curr_tree, first_cluster_mrca_node, second_cluster_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
	
	return Curr_tree

#-------------------------------------------
"""
this function processes one internal node (basically the children list)
to resolve multifurcation
"""
def ResolveMultifurcation_old(Curr_tree, clust_species_list, no_of_input_clusters, Output_Text_File, NJ_RULE_USED, DIST_MAT_TYPE):
	# total number of clusters
	no_of_clust = no_of_input_clusters

	# allocate a 2D square matrix of no_of_clust dimension
	DistMat = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.float)
	Norm_DistMat = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.float)

	# here we compute the ILS score of the current cluster pair
	# with respect to input gene tree list
	if (DEBUG_LEVEL > 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n Examining ILS score for individual cluster pairs ')
		fp.close()      

	#---------------------------------------
	## using single taxon as a representative of the taxa cluster
	#Fill_DistMat_SingleEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE)

	# using average information from a taxa cluster
	Fill_DistMat_AvgEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE)
	#---------------------------------------

	# loop to execute the agglomerative clustering
	while(no_of_clust > 2):
		min_idx_i, min_idx_j = Get_NJ_Based_Min_Pair_Idx(DistMat, Norm_DistMat, no_of_clust, clust_species_list, NJ_RULE_USED, Output_Text_File)

		# note down the taxa list in these two indices of the clust_species_list
		taxa_list = []
		for x in clust_species_list[min_idx_i]:
			taxa_list.append(x)
		for x in clust_species_list[min_idx_j]:
			taxa_list.append(x)

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n\n\n min_idx_i ' + str(min_idx_i) + ' min_idx_j : ' + str(min_idx_j))
			fp.write('\n min_idx_i species list ' + str(clust_species_list[min_idx_i]))
			fp.write('\n min_idx_j species list ' + str(clust_species_list[min_idx_j]))
			fp.write('\n complete species list ' + str(taxa_list))
			fp.close()
			
		"""
		now we merge the pair of clusters pointed by these indices
		"""
		Curr_tree = Merge_Cluster_Pair(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File)
		#---------------------------------------------------------      
		# remove individual clusters' taxa information from the clust_species_list
		# and add taxa_list as a new element
		clust_species_list.pop(min_idx_i)
		clust_species_list.pop(min_idx_j - 1)
		
		# comment - sourya
		clust_species_list.append(taxa_list)    
		
		## add - sourya
		#taxa_list_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
		#preorder_taxa_list = []
		#for n in taxa_list_mrca_node.preorder_iter():
			#if (n.is_leaf() == True):
				#if n.taxon.label in taxa_list:
					#preorder_taxa_list.append(n.taxon.label)
		#clust_species_list.append(preorder_taxa_list)   
		## end add - sourya
		#---------------------------------------------------------      
		"""
		adjust the DistMat by inserting one new row and column corresponding to the new cluster
		and then deleting the information of earlier two clusters
		"""
		# first append one row
		DistMat = numpy.vstack((DistMat, numpy.zeros((1, no_of_clust), dtype=numpy.float)))
		# then append one column
		DistMat = numpy.hstack((DistMat, numpy.zeros((no_of_clust + 1, 1), dtype=numpy.float)))
		# now apply reshape operation to get proper square matrix dimension
		DistMat = numpy.reshape(DistMat, ((no_of_clust + 1), (no_of_clust + 1)), order='C')
		
		# now fill the elements of the new added row and column
		for k in range(no_of_clust):
			if (NJ_RULE_USED == AGGLO_CLUST):
				DistMat[k][no_of_clust] = (DistMat[k][min_idx_i] + DistMat[k][min_idx_j]) / 2.0
			else:
				DistMat[k][no_of_clust] = (DistMat[k][min_idx_i] + DistMat[k][min_idx_j] - DistMat[min_idx_i][min_idx_j]) / 2.0
			# symmetric property
			DistMat[no_of_clust][k] = DistMat[k][no_of_clust]
		
		# now remove the rows and columns corresponding to min_idx_i and min_idx_j
		DistMat = numpy.delete(DistMat, (min_idx_i), axis=0)	# delete the row
		DistMat = numpy.delete(DistMat, (min_idx_i), axis=1)	# delete the column
		DistMat = numpy.delete(DistMat, (min_idx_j - 1), axis=0)	# delete the row
		DistMat = numpy.delete(DistMat, (min_idx_j - 1), axis=1)	# delete the column

		# clear Norm_DistMat
		Norm_DistMat = numpy.delete(Norm_DistMat, (min_idx_i), axis=0)	# delete the row
		Norm_DistMat = numpy.delete(Norm_DistMat, (min_idx_i), axis=1)	# delete the column
		Norm_DistMat.fill(0)
		
		# decrement the number of clusters considered
		no_of_clust = no_of_clust - 1
	
	# add - sourya
	"""
	delete the distance matrices
	"""
	del DistMat
	del Norm_DistMat
	# end add - sourya
	
	return

#-------------------------------------------
"""
here we process all the input taxa of all the taxa clusters, and club their supporting tree information
"""
def Form_Complete_Supporting_TreeDict(SpecList, Complete_Supporting_Tree_Dict, no_of_clust, outfile):
	for i in range(no_of_clust):
		for j in range(len(SpecList[i])):
			curr_taxa = SpecList[i][j]
			curr_taxa_idx = COMPLETE_INPUT_TAXA_LIST.index(curr_taxa)
			curr_taxa_support_tree = Taxa_Info_Dict[curr_taxa_idx]._GetSupportTreeList()
			if (DEBUG_LEVEL > 2):
				fp = open(outfile, 'a')
				fp.write('\n Form_Complete_Supporting_TreeDict -- taxa cluster index: ' + str(i) + \
					' curr_taxa_idx: ' + str(curr_taxa_idx) + ' curr_taxa: ' + str(curr_taxa) + ' supporting tree set: ' + str(curr_taxa_support_tree)) 
				fp.close()
			"""
			for individual trees within the set "curr_taxa_support_tree"
			initiate the entry of "Complete_Supporting_Tree_Dict"
			"""
			for t in curr_taxa_support_tree:
				if t not in Complete_Supporting_Tree_Dict:
					"""
					initiate the dictionary entry
					the value is an empty list array of "no_of_clust" dimension
					"""
					Complete_Supporting_Tree_Dict.setdefault(t, [[] for k in range(no_of_clust)])
				"""
				now add the "j" in the "i'th" list of the dictionary entry
				"""
				Complete_Supporting_Tree_Dict[t][i].append(j)

	return

#-------------------------------------------
"""
this function updates the distance matrices:
Based on branch count for individual couplets

existing_taxa_cluster_list: contains the list of taxa clusters those are supported by 
the current tree ("tree_idx")
"""
def Update_DistanceMat_Couplet_Measures(Complete_Supporting_Tree_Dict, existing_taxa_cluster_list, \
	tree_idx, Coal_Rank_DistMat, XL_DistMat, Valid_Mat, SpecList, curr_Inp_tree, outfile):

	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')

	"""
	this is the array of taxa labels corresponding to the existing_taxa_cluster_list
	"""
	Taxa_Label_Existing_Taxa_Cluster_List = [[] for k in range(len(existing_taxa_cluster_list))]
	for taxa_clust_idx in existing_taxa_cluster_list:
		for j in range(len(Complete_Supporting_Tree_Dict[tree_idx][taxa_clust_idx])):
			taxa_idx = Complete_Supporting_Tree_Dict[tree_idx][taxa_clust_idx][j]
			Taxa_Label_Existing_Taxa_Cluster_List[existing_taxa_cluster_list.index(taxa_clust_idx)].append(SpecList[taxa_clust_idx][taxa_idx])
	
	if (DEBUG_LEVEL > 2):
		fp.write('\n\n ===>>> Analyzing input treelist index: ' + str(tree_idx))
		for i in range(len(existing_taxa_cluster_list)):
			fp.write('\n\n existing_taxa_cluster_list index: ' + str(i) + ' val: ' + str(existing_taxa_cluster_list[i]) + \
				'  Taxa_Label_Existing_Taxa_Cluster_List: ' + str(Taxa_Label_Existing_Taxa_Cluster_List[i]))
	
	"""
	this is the current tree for analysis
	complete set of taxa (to be used for restriction operation)
	"""
	curr_tree_complete_taxa_list = []
	for i in range(len(existing_taxa_cluster_list)):
		curr_tree_complete_taxa_list.extend(Taxa_Label_Existing_Taxa_Cluster_List[i])

	if (DEBUG_LEVEL > 2):
		fp.write('\n\n ===>>> The initial tree: ' + str(curr_Inp_tree))
		fp.write('\n\n ===>>> curr_tree_complete_taxa_list: ' + str(curr_tree_complete_taxa_list))
	
	"""
	first prune the input tree using only the taxa set belonging to "curr_tree_complete_taxa_list"
	"""
	curr_Inp_tree.retain_taxa_with_labels(curr_tree_complete_taxa_list)	#, update_splits=False)
	#curr_Inp_tree.update_splits(delete_outdegree_one=False)
	if (DEBUG_LEVEL >= 2):
		fp.write('\n\n ===>>> Current tree with only "curr_tree_complete_taxa_list" retained: ' + str(curr_Inp_tree))
	
	"""
	pre-compute the LCA nodes for all different taxa clusters, and store them in the below mentioned list
	"""
	Array_of_LCA_nodes = []
	for i in range(len(existing_taxa_cluster_list)):
		"""
		LCA of the taxa cluster whose species list is maintained in "Taxa_Label_Existing_Taxa_Cluster_List[i]"
		if a single taxon is present at Taxa_Label_Existing_Taxa_Cluster_List[i]
		then LCA node returns the parent (internal) node of the corresponding taxon
		"""
		if (len(Taxa_Label_Existing_Taxa_Cluster_List[i]) > 1):
			"""
			here the LCA node is an internal node
			"""
			curr_LCA_node = curr_Inp_tree.mrca(taxon_labels=Taxa_Label_Existing_Taxa_Cluster_List[i])
		else:
			"""
			here the LCA node is a leaf - so we use its parent (internal) node as the target LCA node
			"""
			curr_LCA_node = (curr_Inp_tree.find_node_with_taxon_label(Taxa_Label_Existing_Taxa_Cluster_List[i][0])).parent_node
		"""
		add the LCA node in the Array_of_LCA_nodes
		"""
		Array_of_LCA_nodes.append(curr_LCA_node)
		if (DEBUG_LEVEL > 2):
			fp.write('\n\n ===>>> LCA node of taxa cluster index: ' + str(i) + ' taxa cluster: ' + str(existing_taxa_cluster_list[i]) + ' is: ' + \
				str(curr_LCA_node) + '  its level: ' + str(curr_LCA_node.level()))

	"""
	now modify the tree so that the taxa cluster representative taxon are placed as children of the 
	pre computed LCA nodes
	"""
	for i in range(len(existing_taxa_cluster_list)):
		"""
		create a new node having the taxon label = existing_taxa_cluster_list[i]
		"""
		curr_LCA_node = Array_of_LCA_nodes[i]
		new_taxon_label = 'A' + str(existing_taxa_cluster_list[i])
		curr_LCA_node.new_child(taxon=Taxon(label=new_taxon_label)) 
	
	"""
	first print the new tree
	"""
	if (DEBUG_LEVEL > 2):
		fp.write('\n\n ===>>> After appending taxa labels -- modified tree: ' + str(curr_Inp_tree))
	
	"""
	then prune all the earlier taxa labels so that only the new representatives remain
	"""
	curr_Inp_tree.prune_taxa_with_labels(curr_tree_complete_taxa_list)	#, update_splits=False)
	#curr_Inp_tree.update_splits(delete_outdegree_one=False)
	curr_tree_taxa_label_list = curr_Inp_tree.infer_taxa().labels()
	
	if (DEBUG_LEVEL >= 2):
		fp.write('\n\n ===>>> Restricted tree: ' + str(curr_Inp_tree))
		fp.write('\n\n ===>>> Its labels: ' + str(curr_tree_taxa_label_list))
	
	#-----------------------------------------------------
	"""
	now process individual couplets of the tree, and update various input distance matrices 
	accordingly
	"""
	for i in range(len(curr_tree_taxa_label_list) - 1):
		for j in range(i+1, len(curr_tree_taxa_label_list)):
			"""
			compute the LCA between the taxa cluster placed at the indices i and j of the "curr_tree_taxa_label_list"
			"""
			node_i = curr_Inp_tree.find_node_with_taxon_label(curr_tree_taxa_label_list[i])
			node_j = curr_Inp_tree.find_node_with_taxon_label(curr_tree_taxa_label_list[j])
			node_i_level = node_i.level()
			node_j_level = node_j.level()

			node_i_label = node_i.taxon.label
			node_j_label = node_j.taxon.label

			"""
			important - sourya
			somehow, the LCA operation was not working
			so, we have used the custom LCA function
			"""
			curr_LCA_node = Find_MRCA(curr_Inp_tree, [node_i_label, node_j_label])
			#curr_LCA_node = curr_Inp_tree.mrca(taxon_labels=[node_i_label, node_j_label])
			
			curr_LCA_level = curr_LCA_node.level()

			coal_rank_val = (len(curr_tree_taxa_label_list) - curr_LCA_level) * 1.0 / len(curr_tree_taxa_label_list)
			XL_val = (len(curr_LCA_node.leaf_nodes()) - 2) * 1.0 / len(curr_tree_taxa_label_list)
			
			"""
			distance matrix indices and debug
			"""
			distmat_idx_i = existing_taxa_cluster_list[i]
			distmat_idx_j = existing_taxa_cluster_list[j]
			if (DEBUG_LEVEL > 2):
				fp.write('\n\n ===>>> Couplet indices: (' + str(i) + ',' + str(j) + \
					') DistMat indices: (' + str(distmat_idx_i) + ',' + str(distmat_idx_j) + ') XL value: ' + str(XL_val) + \
						'  coal_rank_val: ' + str(coal_rank_val))
			
			Coal_Rank_DistMat[distmat_idx_i][distmat_idx_j] = max(Coal_Rank_DistMat[distmat_idx_i][distmat_idx_j], coal_rank_val)
			Coal_Rank_DistMat[distmat_idx_j][distmat_idx_i] = Coal_Rank_DistMat[distmat_idx_i][distmat_idx_j]
			
			XL_DistMat[distmat_idx_i][distmat_idx_j] = max(XL_DistMat[distmat_idx_i][distmat_idx_j], XL_val)
			XL_DistMat[distmat_idx_j][distmat_idx_i] = XL_DistMat[distmat_idx_i][distmat_idx_j]
			
			"""
			also set the valid matrix
			"""
			Valid_Mat[distmat_idx_i][distmat_idx_j] = 1
			Valid_Mat[distmat_idx_j][distmat_idx_i] = 1
			

	if (DEBUG_LEVEL >= 2):
		fp.close()
	
	"""
	delete the variables
	"""
	Array_of_LCA_nodes = []
	
	return

#-------------------------------------------
"""
new added function for refinement of the supertree in binary form
"""
def ResolveMultifurcation_Couplet_Measures(LCA_node, Input_Treelist, Curr_Suptree, SpecList, nclust, NJ_RULE_USED, outfile):
	"""
	initialize the number of clusters - it is a dynamic variable
	"""
	no_of_clust = nclust
	
	"""
	first create the distance matrix which will contain the pairwise distance between the pair of clusters
	Note: Here we create three distance matrices (for experimentation):
	1) Containing the coalescence rank between individual couplets
	2) Containing the branch count between individual couplets
	3) Containing the XL measures between individual couplets
	"""
	Coal_Rank_DistMat = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.float)
	XL_DistMat = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.float)
	Valid_Mat = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.int)
	
	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')
		fp.write('\n\n\n ******************* \n --> within function ResolveMultifurcation_Latest --- \n **************** \n Clust Species list: ') 
		for i in range(len(SpecList)): 
			fp.write('\n Index : ' + str(i) + ' --->> Taxa list: ' + str(SpecList[i]))
		fp.close()

	"""
	this dictionary stores the support tree information for all the taxa for all the taxa clusters
	used in this refinement stage
	"""
	Complete_Supporting_Tree_Dict = dict()
	Form_Complete_Supporting_TreeDict(SpecList, Complete_Supporting_Tree_Dict, no_of_clust, outfile)
	
	"""
	print the tree dictionary
	"""
	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')
		for t in Complete_Supporting_Tree_Dict:
			fp.write('\n\n Printing tree dictionary elements for the input tree index: ' + str(t))
			for i in range(no_of_clust):
				fp.write('\n Taxa cluster index: ' + str(i) + '  Underlying taxa index list: ' + str(Complete_Supporting_Tree_Dict[t][i]))
		fp.close()
	
	"""
	if a tree t in this dictionary contains more than two non empty lists (corresponding to more than two different taxa clusters)
	then that tree needs to be included in the global supporting tree set T
	"""
	for t in Complete_Supporting_Tree_Dict:
		"""
		no of different taxa clusters supported by this tree
		"""
		existing_taxa_cluster_list = []
		
		for cl_idx in range(no_of_clust):
			if (len(Complete_Supporting_Tree_Dict[t][cl_idx]) > 0):
				existing_taxa_cluster_list.append(cl_idx)
		
		if (len(existing_taxa_cluster_list) > 2):
			if (DEBUG_LEVEL >= 2):
				fp = open(outfile, 'a')
				fp.write('\n\n\n ****** Supporting tree: ' + str(t) + '  needs to be checked ') 
				fp.close()
			"""
			store the current source tree in a separate structure, for processing
			"""
			tree_idx = int(t)
			curr_Inp_tree = dendropy.Tree(Input_Treelist[tree_idx])
			"""
			update the global distance matrix entries for this supporting tree
			"""
			Update_DistanceMat_Couplet_Measures(Complete_Supporting_Tree_Dict, existing_taxa_cluster_list, \
				tree_idx, Coal_Rank_DistMat, XL_DistMat, Valid_Mat, SpecList, curr_Inp_tree, outfile)

	"""
	allocate the distance matrices which contains the branch count distance and the relative 
	distance entries for individual couplets
	"""
	DistMat = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.float)
	Norm_DistMat = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.float)

	"""
	now assign individual elements of the DistMat
	"""
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			"""
			individual distance matrix entry is the sum of distance 
			from individual leaf nodes to their LCA node
			"""
			#DistMat[i][j] = XL_DistMat[i][j]
			DistMat[i][j] = Coal_Rank_DistMat[i][j]	#sourya
			DistMat[j][i] = DistMat[i][j]
	
	#-------------------------------------------------------
	"""
	now we loop through the distance matrices to refine the supertree using NJ based 
	agglomeration
	"""
	#-------------------------------------------------------
	while(no_of_clust > 2):
		"""
		extract the current minimum for NJ based agglomeration
		"""
		#min_idx_i, min_idx_j = Get_NJ_Based_Min_Pair_Idx_Valid(DistMat, Norm_DistMat, Valid_Mat, no_of_clust, SpecList, AGGLO_CLUST, outfile)
		min_idx_i, min_idx_j = Get_NJ_Based_Min_Pair_Idx_Valid(DistMat, Norm_DistMat, Valid_Mat, no_of_clust, SpecList, NJ_RULE_USED, outfile)	#sourya
		
		"""
		if no entry is valid
		"""
		if (min_idx_i == -1) and (min_idx_j == -1):
			break

		# note down the taxa list in these two indices of the SpecList
		taxa_list = []
		for x in SpecList[min_idx_i]:
			taxa_list.append(x)
		for x in SpecList[min_idx_j]:
			taxa_list.append(x)

		if (DEBUG_LEVEL >= 2):
			fp = open(outfile, 'a')
			fp.write('\n min_idx_i ' + str(min_idx_i) + ' min_idx_j : ' + str(min_idx_j))
			fp.write('\n min_idx_i species list ' + str(SpecList[min_idx_i]))
			fp.write('\n min_idx_j species list ' + str(SpecList[min_idx_j]))
			fp.write('\n complete species list ' + str(taxa_list))
			fp.close()
			
		"""
		now we merge the pair of clusters pointed by these indices
		"""
		Curr_Suptree = Merge_Cluster_Pair(Curr_Suptree, SpecList, min_idx_i, min_idx_j, taxa_list, outfile)
		#---------------------------------------------------------      
		"""
		remove individual clusters' taxa information from the SpecList
		and add taxa_list as a new element
		"""
		SpecList.pop(min_idx_i)
		SpecList.pop(min_idx_j - 1)
		
		SpecList.append(taxa_list)    
		#---------------------------------------------------------      
		"""
		adjust the DistMat by inserting one new row and column corresponding to the new cluster
		and then deleting the information of earlier two clusters
		"""
		# first append one row
		DistMat = numpy.vstack((DistMat, numpy.zeros((1, no_of_clust), dtype=numpy.float)))
		# then append one column
		DistMat = numpy.hstack((DistMat, numpy.zeros((no_of_clust + 1, 1), dtype=numpy.float)))
		# now apply reshape operation to get proper square matrix dimension
		DistMat = numpy.reshape(DistMat, ((no_of_clust + 1), (no_of_clust + 1)), order='C')
		
		"""
		now fill the elements of the new added row and column
		"""
		for k in range(no_of_clust):
			if (NJ_RULE_USED == TRADITIONAL_NJ):
				"""
				standard NJ rule based updation, along with the symmetric property
				"""
				DistMat[k][no_of_clust] = (DistMat[k][min_idx_i] + DistMat[k][min_idx_j] - DistMat[min_idx_i][min_idx_j]) / 2.0
			else:
				DistMat[k][no_of_clust] = (DistMat[k][min_idx_i] + DistMat[k][min_idx_j]) / 2.0
			"""
			maintain the symmetric property
			"""
			DistMat[no_of_clust][k] = DistMat[k][no_of_clust]

		"""
		now remove the rows and columns corresponding to min_idx_i and min_idx_j
		"""
		DistMat = numpy.delete(DistMat, (min_idx_i), axis=0)	# delete the row
		DistMat = numpy.delete(DistMat, (min_idx_i), axis=1)	# delete the column
		DistMat = numpy.delete(DistMat, (min_idx_j - 1), axis=0)	# delete the row
		DistMat = numpy.delete(DistMat, (min_idx_j - 1), axis=1)	# delete the column

		"""
		clear Norm_DistMat
		"""
		Norm_DistMat = numpy.delete(Norm_DistMat, (min_idx_i), axis=0)	# delete the row
		Norm_DistMat = numpy.delete(Norm_DistMat, (min_idx_i), axis=1)	# delete the column
		Norm_DistMat.fill(0)
		
		"""
		decrement the number of clusters considered
		"""
		no_of_clust = no_of_clust - 1
	
	"""
	delete the distance matrices
	"""
	del XL_DistMat
	del Coal_Rank_DistMat
	#del Branch_Count_DistMat
	#del Reln_Mat
	del DistMat
	del Norm_DistMat
	
	Complete_Supporting_Tree_Dict.clear()

	return Curr_Suptree

#-------------------------------------------
"""
this function updates the distance matrices:
Based on branch count for individual couplets

existing_taxa_cluster_list: contains the list of taxa clusters those are supported by 
the current tree ("tree_idx")
"""
def Update_DistanceMat_New(Restricted_Treelist, Complete_Supporting_Tree_Dict, existing_taxa_cluster_list, \
	tree_idx, SpecList, curr_Inp_tree, outfile):

	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')

	"""
	this is the array of taxa labels corresponding to the existing_taxa_cluster_list
	"""
	Taxa_Label_Existing_Taxa_Cluster_List = [[] for k in range(len(existing_taxa_cluster_list))]
	for taxa_clust_idx in existing_taxa_cluster_list:
		for j in range(len(Complete_Supporting_Tree_Dict[tree_idx][taxa_clust_idx])):
			taxa_idx = Complete_Supporting_Tree_Dict[tree_idx][taxa_clust_idx][j]
			Taxa_Label_Existing_Taxa_Cluster_List[existing_taxa_cluster_list.index(taxa_clust_idx)].append(SpecList[taxa_clust_idx][taxa_idx])
	
	if (DEBUG_LEVEL > 2):
		fp.write('\n\n ===>>> Analyzing input treelist index: ' + str(tree_idx))
		for i in range(len(existing_taxa_cluster_list)):
			fp.write('\n\n existing_taxa_cluster_list index: ' + str(i) + ' val: ' + str(existing_taxa_cluster_list[i]) + \
				'  Taxa_Label_Existing_Taxa_Cluster_List: ' + str(Taxa_Label_Existing_Taxa_Cluster_List[i]))
	
	"""
	this is the current tree for analysis
	complete set of taxa (to be used for restriction operation)
	"""
	curr_tree_complete_taxa_list = []
	for i in range(len(existing_taxa_cluster_list)):
		curr_tree_complete_taxa_list.extend(Taxa_Label_Existing_Taxa_Cluster_List[i])

	if (DEBUG_LEVEL > 2):
		fp.write('\n\n ===>>> The initial tree: ' + str(curr_Inp_tree))
		fp.write('\n\n ===>>> curr_tree_complete_taxa_list: ' + str(curr_tree_complete_taxa_list))
	
	"""
	first prune the input tree using only the taxa set belonging to "curr_tree_complete_taxa_list"
	"""
	curr_Inp_tree.retain_taxa_with_labels(curr_tree_complete_taxa_list)	#, update_splits=False)
	#curr_Inp_tree.update_splits(delete_outdegree_one=False)
	if (DEBUG_LEVEL >= 2):
		fp.write('\n\n ===>>> Current tree with only "curr_tree_complete_taxa_list" retained: ' + str(curr_Inp_tree))
	
	"""
	here we first create few taxon, labeled by the taxa clusters in "existing_taxa_cluster_list"
	say the cluster is 1
	we insert the taxon 1 as a child of the LCA node of the taxa labels belonging to the taxa cluster 1
	similarly for other taxa clusters
	The objective is to reduce the current tree to have a single representative per distinct taxa cluster
	"""

	"""
	pre-compute the LCA nodes for all different taxa clusters, and store them in the below mentioned list
	"""
	Array_of_LCA_nodes = []
	for i in range(len(existing_taxa_cluster_list)):
		"""
		LCA of the taxa cluster whose species list is maintained in "Taxa_Label_Existing_Taxa_Cluster_List[i]"
		if a single taxon is present at Taxa_Label_Existing_Taxa_Cluster_List[i]
		then LCA node returns the parent (internal) node of the corresponding taxon
		"""
		if (len(Taxa_Label_Existing_Taxa_Cluster_List[i]) > 1):
			"""
			here the LCA node is an internal node
			"""
			curr_LCA_node = curr_Inp_tree.mrca(taxon_labels=Taxa_Label_Existing_Taxa_Cluster_List[i])
		else:
			"""
			here the LCA node is a leaf - so we use its parent (internal) node as the target LCA node
			"""
			curr_LCA_node = (curr_Inp_tree.find_node_with_taxon_label(Taxa_Label_Existing_Taxa_Cluster_List[i][0])).parent_node
		"""
		add the LCA node in the Array_of_LCA_nodes
		"""
		Array_of_LCA_nodes.append(curr_LCA_node)
		if (DEBUG_LEVEL > 2):
			fp.write('\n\n ===>>> LCA node of taxa cluster index: ' + str(i) + ' taxa cluster: ' + str(existing_taxa_cluster_list[i]) + ' is: ' + \
				str(curr_LCA_node) + '  its level: ' + str(curr_LCA_node.level()))

	"""
	now modify the tree so that the taxa cluster representative taxon are placed as children of the 
	pre computed LCA nodes
	"""
	for i in range(len(existing_taxa_cluster_list)):
		"""
		create a new node having the taxon label = existing_taxa_cluster_list[i]
		"""
		curr_LCA_node = Array_of_LCA_nodes[i]
		new_taxon_label = 'A' + str(existing_taxa_cluster_list[i])
		curr_LCA_node.new_child(taxon=Taxon(label=new_taxon_label)) 
	
	"""
	first print the new tree
	"""
	if (DEBUG_LEVEL > 2):
		fp.write('\n\n ===>>> After appending taxa labels -- modified tree: ' + str(curr_Inp_tree))
	
	"""
	then prune all the earlier taxa labels so that only the new representatives remain
	"""
	curr_Inp_tree.prune_taxa_with_labels(curr_tree_complete_taxa_list)	#, update_splits=False)
	#curr_Inp_tree.update_splits(delete_outdegree_one=False)
	curr_tree_taxa_label_list = curr_Inp_tree.infer_taxa().labels()
	
	if (DEBUG_LEVEL >= 2):
		fp.write('\n\n ===>>> Restricted tree: ' + str(curr_Inp_tree))
		fp.write('\n\n ===>>> Its labels: ' + str(curr_tree_taxa_label_list))
	
	"""
	append the restricted input tree (with respect to the set of taxa clusters)
	to the treelist
	"""
	Restricted_Treelist.append(curr_Inp_tree)
	
	if (DEBUG_LEVEL >= 2):
		fp.close()
	
	"""
	delete the variables
	"""
	Array_of_LCA_nodes = []
	
	return

#-------------------------------------------
"""
new added function for refinement of the supertree in binary form
"""
def ResolveMultifurcation_Latest(LCA_node, Input_Treelist, Curr_Suptree, SpecList, nclust, outfile, nnode):
	"""
	initialize the number of clusters - it is a dynamic variable
	"""
	no_of_clust = nclust
	
	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')
		fp.write('\n\n\n ******************* \n --> within function ResolveMultifurcation_Latest --- \n **************** \n Clust Species list: ') 
		for i in range(len(SpecList)): 
			fp.write('\n Index : ' + str(i) + ' --->> Taxa list: ' + str(SpecList[i]))
		fp.close()

	"""
	this dictionary stores the support tree information for all the taxa for all the taxa clusters
	used in this refinement stage
	"""
	Complete_Supporting_Tree_Dict = dict()
	Form_Complete_Supporting_TreeDict(SpecList, Complete_Supporting_Tree_Dict, no_of_clust, outfile)
	
	"""
	print the tree dictionary
	"""
	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')
		for t in Complete_Supporting_Tree_Dict:
			fp.write('\n\n Printing tree dictionary elements for the input tree index: ' + str(t))
			for i in range(no_of_clust):
				fp.write('\n Taxa cluster index: ' + str(i) + '  Underlying taxa index list: ' + str(Complete_Supporting_Tree_Dict[t][i]))
		fp.close()
	
	"""
	this is a treelist containing the input trees restricted to the taxa subsets 
	analyzed here
	"""
	Restricted_Treelist = TreeList()
	
	Taxa_Clusters_Covered_Total = []
	
	"""
	if a tree t in this dictionary contains more than two non empty lists (corresponding to more than two different taxa clusters)
	then that tree needs to be included in the global supporting tree set T
	"""
	for t in Complete_Supporting_Tree_Dict:
		"""
		no of different taxa clusters supported by this tree
		"""
		existing_taxa_cluster_list = []
		
		for cl_idx in range(no_of_clust):
			if (len(Complete_Supporting_Tree_Dict[t][cl_idx]) > 0):
				existing_taxa_cluster_list.append(cl_idx)
		
		if (len(existing_taxa_cluster_list) > 2):
			for x in existing_taxa_cluster_list:
				if x not in Taxa_Clusters_Covered_Total:
					Taxa_Clusters_Covered_Total.append(x)
			
			if (DEBUG_LEVEL >= 2):
				fp = open(outfile, 'a')
				fp.write('\n\n\n ****** Supporting tree: ' + str(t) + '  needs to be checked ') 
				fp.close()
			"""
			store the current source tree in a separate structure, for processing
			"""
			tree_idx = int(t)
			curr_Inp_tree = dendropy.Tree(Input_Treelist[tree_idx])
			"""
			update the global distance matrix entries for this supporting tree
			"""
			Update_DistanceMat_New(Restricted_Treelist, Complete_Supporting_Tree_Dict, existing_taxa_cluster_list, \
				tree_idx, SpecList, curr_Inp_tree, outfile)

	#------------------------------------------------------------------------
	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')

	Taxa_Clusters_Covered_Total.sort()
	fp.write('\n\n\n ************ Taxa_Clusters_Covered_Total: ' + str(Taxa_Clusters_Covered_Total) + '**************')

	temp_outfile = 'temp.txt'

	"""
	now write the restricted treelist in a file
	"""
	treefile = 'input_treelist_' + str(nnode) + '_nexus.tre'
	if (len(Restricted_Treelist) > 0):
		if (DEBUG_LEVEL >= 2):
			fp.write('\n Writing the restricted treelist to the file: ' + str(treefile)) 
		Restricted_Treelist.write_to_path(treefile, 'nexus', suppress_rooting=True, simple=True)
		"""
		now process the nexus file and derive a supertree from the nexus treelist file
		"""
		SupTree_command = 'java -Xmx1g -jar Triplet_Supertree_Lin.jar ' + str(treefile)
		os.system(SupTree_command)
	
		"""
		the supertree will be stored in the following file name
		read the supertree
		the supertree is saved in nexus format
		"""
		supertree_filename = 'super_' + treefile
		SupTree_Bin = dendropy.Tree.get_from_path(supertree_filename, schema='nexus')
		if (DEBUG_LEVEL >= 2):
			fp.write('\n\n SupTree_Bin (supertree with respect to taxa cluster indices): ' + str(SupTree_Bin)) 
	
		"""
		extract the root node of the SupTree_Bin
		"""
		for n in SupTree_Bin.preorder_internal_node_iter():
			Suptree_Bin_Root_node = n
			break
	
		"""
		set of taxa labels belonging to the supertree obtained by triplet based approach
		"""
		SupTree_Bin_taxa_Labels = SupTree_Bin.infer_taxa().labels()
		if (DEBUG_LEVEL >= 2):
			fp.write('\n\n Label of Suptree_Bin_Root_node: ' + str(Node_Label(Suptree_Bin_Root_node)))
	
		"""
		insert the "Suptree_Bin_Root_node" as a child to the "LCA_node" of the input tree
		"""
		LCA_node.add_child(Suptree_Bin_Root_node)
		Suptree_Bin_Root_node.parent_node = LCA_node
		if (DEBUG_LEVEL >= 2):
			fp.write('\n\n After appending the Suptree_Bin_Root_node as a child to the LCA node, the input supertree becomes: ' + str(Curr_Suptree))
	
		"""
		now for individual taxa clusters within SpecList (denoted as SpecList[i]) where i ranges from 0 to (nclust - 1)
		1) Generate a tree copy from Curr_Suptree. Say it is termed as "Curr_Suptree_Copy"
		2) Restrict "Curr_Suptree_Copy" to SpecList[i] 
		3) Prune from the "Curr_Suptree", the taxa list SpecList[i]
		4) In "Curr_Suptree", search for the node n corresponding to the taxon i
		5) Place Curr_Suptree_Copy (its root) as a child to the parent node of n
		"""
		for i in Taxa_Clusters_Covered_Total:	#range(no_of_clust):
			Curr_Suptree_Copy = Tree(Curr_Suptree)
			Curr_Suptree_Copy.write_to_path(temp_outfile, 'newick')
			Curr_Suptree_Copy = dendropy.Tree.get_from_path(temp_outfile, schema='newick')
			Curr_Suptree_Copy.retain_taxa_with_labels(SpecList[i])
			if (DEBUG_LEVEL >= 2):
				fp.write('\n\n Taxa cluster index: ' + str(i) + '  Corresponding taxa list: ' + str(SpecList[i]))
				fp.write('\n\n Curr_Suptree_Copy after retaining only this taxa subset: ' + str(Curr_Suptree_Copy))

			Curr_Suptree.write_to_path(temp_outfile, 'newick')
			Curr_Suptree = dendropy.Tree.get_from_path(temp_outfile, schema='newick')
			Curr_Suptree.prune_taxa_with_labels(SpecList[i])
			if (DEBUG_LEVEL >= 2):
				fp.write('\n\n After Pruning the taxa set, Curr_Suptree: ' + str(Curr_Suptree))
			"""
			find the label l in Curr_Suptree
			where the label l is 'A' + str(i)
			"""
			l = 'A' + str(i)
			n = Curr_Suptree.find_node_with_taxon_label(l)
			if (DEBUG_LEVEL >= 2):
				fp.write('\n\n The node label corresponding to the taxon label ' + str(l) + ' in the Curr_Suptree: ' + str(Node_Label(n)))
			
			if (len(SpecList[i]) > 1):
				for n1 in Curr_Suptree_Copy.preorder_internal_node_iter():
					Curr_Suptree_Copy_Root_node = n1
					break
				n.parent_node.add_child(Curr_Suptree_Copy_Root_node)
			else:
				Curr_Suptree_Copy_node = Curr_Suptree_Copy.find_node_with_taxon_label(SpecList[i][0])
				n.parent_node.add_child(Curr_Suptree_Copy_node)

			if (DEBUG_LEVEL >= 2):
				fp.write('\n\n After Re-inserting the taxa list: Curr_Suptree: ' + str(Curr_Suptree))
			
			Curr_Suptree.write_to_path(temp_outfile, 'newick')
			Curr_Suptree = dendropy.Tree.get_from_path(temp_outfile, schema='newick')
			Curr_Suptree.prune_taxa_with_labels([l])
			if (DEBUG_LEVEL >= 2):
				fp.write('\n\n After Pruning the taxon :' + str(l) + '  Curr_Suptree: ' + str(Curr_Suptree))

	if (DEBUG_LEVEL >= 2):
		fp.close()

	"""
	remove the file 'temp.txt'
	"""
	if os.path.exists(temp_outfile):
		system_cmd = 'rm ' + str(temp_outfile)
		os.system(system_cmd)


	Complete_Supporting_Tree_Dict.clear()
	Restricted_Treelist = []
	Curr_Suptree_Copy = []
	SupTree_Bin = []

	return Curr_Suptree

#-------------------------------------------
"""
this function refines input supertree such that the supertree becomes binary
this is required for proper benchmarking with existing binary tree construction methods on 
ILS sorting
"""
def Refine_Supertree_Binary_Form(Input_Treelist, Curr_Suptree, Output_Text_File, NJ_RULE_USED, DIST_MAT_TYPE):

	"""
	we traverse input tree internal nodes in postorder fashion
	and list the child nodes of it
	if the no of children > 2 then it is a case of multifurcation
	for resolving
	"""
	
	#------------------------------------------------------
	# important  code - sourya
	#------------------------------------------------------
	
	"""
	old version of the binary refinement
	"""
	for curr_node in Curr_Suptree.postorder_internal_node_iter():
		curr_node_children = curr_node.child_nodes()
		if (len(curr_node_children) > 2):
			"""
			create a list which will contain the species list lying under 
			individual child nodes of rhe current node
			"""
			clust_species_list = []
			for x in curr_node_children:
				subl = []
				for n in x.preorder_iter():
					if (n.is_leaf() == True):
						subl.append(n.taxon.label)
				clust_species_list.append(subl)
	
			# call the resolving routine
			ResolveMultifurcation_old(Curr_Suptree, clust_species_list, len(curr_node_children), Output_Text_File, NJ_RULE_USED, DIST_MAT_TYPE)
			#Curr_Suptree = ResolveMultifurcation_Couplet_Measures(curr_node, Input_Treelist, Curr_Suptree, \
				#clust_species_list, len(curr_node_children), NJ_RULE_USED, Output_Text_File)
	
	#------------------------------------------------------
	# important  code - sourya
	#------------------------------------------------------
	
	##------------------------------------------------
	## contains all the taxa list of all multifurcating nodes
	#Global_Clust_Species_List = []
	
	#for curr_node in Curr_Suptree.postorder_internal_node_iter():
		#curr_node_children = curr_node.child_nodes()
		#if (len(curr_node_children) > 2):
			#"""
			#create a list which will contain the species list lying under 
			#individual child nodes of rhe current node
			#"""
			#clust_species_list = []
			#for x in curr_node_children:
				#subl = []
				#for n in x.preorder_iter():
					#if (n.is_leaf() == True):
						#subl.append(n.taxon.label)
				#clust_species_list.append(subl)
			#"""
			#append the clust_species_list in the list "Global_Clust_Species_List"
			#"""
			#Global_Clust_Species_List.append(clust_species_list)
	
	#"""
	#now navigate through individual elements of Global_Clust_Species_List
	#"""
	#for nnode in range(len(Global_Clust_Species_List)):
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n\n *** Examining the index ' + str(nnode) + '  of Global_Clust_Species_List')
			#fp.write('\n The set of taxa for analysis: ' + str(Global_Clust_Species_List[nnode]))
			#fp.write('\n Degree of multifurcation: ' + str(len(Global_Clust_Species_List[nnode])))
			#fp.close()
		
		#clust_species_list = Global_Clust_Species_List[nnode]
		#taxa_list = []
		#for i in range(len(clust_species_list)):
			#taxa_list.extend(clust_species_list[i])

		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n taxa_list: ' + str(taxa_list))
			#fp.close()
		
		#LCA_node = Curr_Suptree.mrca(taxon_labels=taxa_list)

		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n Label of LCA_node: ' + str(Node_Label(LCA_node)))
			#fp.close()
			
		#Curr_Suptree = ResolveMultifurcation_Latest(LCA_node, Input_Treelist, Curr_Suptree, \
			#Global_Clust_Species_List[nnode], len(Global_Clust_Species_List[nnode]), Output_Text_File, nnode)
		
	##------------------------------------------------









	
	#no_of_multifurcating_nodes = 0
	
	#flag = True
	
	#Existing_Taxa_List_Checked = []
	
	#while (flag == True):
		
		#flag = False
		
		#for curr_node in Curr_Suptree.postorder_internal_node_iter():
			#curr_node_children = curr_node.child_nodes()
			#if (len(curr_node_children) > 2):
				#"""
				#create a list which will contain the species list lying under 
				#individual child nodes of rhe current node
				#"""
				#clust_species_list = []
				#for x in curr_node_children:
					#subl = []
					#for n in x.preorder_iter():
						#if (n.is_leaf() == True):
							#subl.append(n.taxon.label)
					#clust_species_list.append(subl)
				
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n\n *** Before soring the clust_species_list: ' + str(clust_species_list))
				
				#clust_species_list.sort()

				#if (DEBUG_LEVEL >= 2):
					#fp.write('\n\n *** After soring the clust_species_list: ' + str(clust_species_list))
					#fp.close()

				#if clust_species_list not in Existing_Taxa_List_Checked:
					#if (DEBUG_LEVEL >= 2):
						#fp = open(Output_Text_File, 'a')
						#fp.write('\n This taxa list was not tested before - check for binary refinement ')
						#fp.close()
					
					#Existing_Taxa_List_Checked.append(clust_species_list)
					#flag = True
					#no_of_multifurcating_nodes = no_of_multifurcating_nodes + 1
					#break
				
			
		#if (flag == True):
			#Curr_Suptree = ResolveMultifurcation_Latest(curr_node, Input_Treelist, Curr_Suptree, clust_species_list, len(curr_node_children), \
				#NJ_RULE_USED, Output_Text_File, no_of_multifurcating_nodes)

	#------------------------------------------------
		
		
	"""
	now delete all the files within the current directory which ends with '_nexus.tre' and 'temp.txt'
	"""
	for f in os.listdir("."):
		if f.endswith("_nexus.tre"):
			sys_cmd = 'rm ' + str(f)
			os.system(sys_cmd)
			

	return Curr_Suptree

