#!/usr/bin/env python

import Header
from Header import *
    
##-----------------------------------------------------    
#""" 
#following code removes extra paranthesis (thus producing insignificant edges) 
#from the supertree expression contained in the string Final_Supertree_Str
#this is important since the original string (and thus the tree) may contain 
#internal nodes of outdegree 1, thus can be removed 
#employing the update_splits routine in the dendropy often causes trouble
#it is better to manually track the number of taxa and the occurrences of 
#enclosing brackets within the string expression
#"""
#def Remove_Extra_Paranthesis(Final_Supertree_Str):
	## first we convert the input string expression in a list format
	## so as to traverse it in position specific manner
	#L = list(Final_Supertree_Str)
	## this is a stack (list structure) which contains the positions of the first brackets encountered
	## within the input string expression
	## this is required to keep track of the correspondence between the opening and the closing brackets
	#SL = []
	## this is a dictionary storing the positions of respective opening and closing brackets
	## the position of opening bracket is the dictionary key
	## the value for one key of the dictionary contains the position of the closing bracket
	#first_bracket_dict = dict()

	## scan through the tree expression string
	#for i in range(len(L)):
		## append the position of an opening bracket in the stack
		#if (L[i] == '('):
			#SL.append(i)
		## for a closing bracket, retrieve the corresponding opening bracket position
		## and enter those positions in the dictionary
		#elif (L[i] == ')'):
			#first_bracket_idx = SL.pop()
			#first_bracket_dict.setdefault(first_bracket_idx, i)

	#if 0:	#(DEBUG_LEVEL > 2):
		#print 'L : ', L
		#print 'first_bracket_dict: ', first_bracket_dict

	#for i in range(len(L) - 1):
		#if (L[i] == '(') and (L[i+1] == '('):
			#sb1 = first_bracket_dict[i]
			#sb2 = first_bracket_dict[i+1]
			#if (sb1 - sb2 == 1):
				#L.pop(i)
				#L.insert(i, '$')
				#L.pop(sb1)
				#L.insert(sb1, '$')

		## add - sourya
		## if enclosing brackets contain only one species 
		## then those brackets need to be removed
		#if (L[i] == '('):
			#sb = first_bracket_dict[i]
			#f = 0
			#for j in range(i,sb):
				#if (L[j] == ','):
					#f = 1
					#break
			#if (f == 0):
				#L.pop(i)
				#L.insert(i, '$')
				#L.pop(sb)
				#L.insert(sb, '$')
		## end add - sourya

	#if 0:	#(DEBUG_LEVEL > 2):
		#print 'L : ', L

	#while (1):
		#if '$' in L:
			#L.remove('$')
		#else:
			#break

	#if 0:	#(DEBUG_LEVEL > 2):
		#print 'L : ', L

	## construct the string containing final supertree
	#outstr = ''.join(L)
	#return outstr
  
#-----------------------------------------------------
# this function prints the tree in Newick format
def PrintNewick(root_clust_node_idx):
	if 0:
		print 'in function printnewick:   root_clust_node_idx: ', root_clust_node_idx
		print 'taxa set: ', Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()  
		print 'out clust list: ', Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1)

	Tree_Str_List = ''
	# process the node provided it has not been explored yet
	if (Cluster_Info_Dict[root_clust_node_idx]._GetExploredStatus() == 0):  
		# set the explored status of the current node to true
		Cluster_Info_Dict[root_clust_node_idx]._SetExploredStatus()
		# get the out edge list of the current node which are not explored yet 
		outnodes = []
		for l in Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1):
			if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				outnodes.append(l)
		# comment - sourya
		if (len(outnodes) == 0):
		# add - sourya
		#if (len(outnodes) <= 1):
			spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
			if (len(spec_list) > 1):
				Tree_Str_List = Tree_Str_List + '('
			Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
			if (len(spec_list) > 1):
				Tree_Str_List = Tree_Str_List + ')'
		else:
			Tree_Str_List = Tree_Str_List + '('
			Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList())
			Tree_Str_List = Tree_Str_List + ','    
			Tree_Str_List = Tree_Str_List + '('
			for i in range(len(outnodes)):
				if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):  
					Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i])
					if (i < (len(outnodes) - 1)):
						# we check whether any subsequent node belonging to the outnodes list
						# is left for traverse
						j = i + 1
						while (j < len(outnodes)):
							if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
								break
							j = j + 1	      
						# in this case, we append one comma
						if (j < len(outnodes)):
							Tree_Str_List = Tree_Str_List + ','
			
			Tree_Str_List = Tree_Str_List + ')'
			Tree_Str_List = Tree_Str_List + ')'
		
	return Tree_Str_List    

#--------------------------------------------------------
"""
this function defines relationship between a pair of nodes in a tree
the relationship is either ancestor / descendant, or siblings, or no relationship 
@parameters: 
	wt_taxa_subset: If True, then the intersection between 
									the und_tax_list and curr_tree_taxa is accounted
	xl_val: excess gene count (normalized) between this couplet
	lca_level: level of the LCA node of this couplet
	curr_tree_taxa: set of taxa belonging to the current gene tree (indices of taxon)
	tr_idx: Index of the input tree
"""
def DefineLeafPairReln(xl_val, lca_level, lca_rank, node1, node2, reln_type, curr_tree_taxa, wt_taxa_subset, tr_idx):

	"""
	compute the levels of individual nodes
	"""
	node1_level = node1.level()
	node2_level = node2.level()

	"""
	using normalized internode count value
	"""
	sum_of_branch_count = (((node1_level - lca_level) + (node2_level - lca_level) - 1) * 1.0) / len(curr_tree_taxa)

	"""
	index of the node1 taxon with respect to COMPLETE_INPUT_TAXA_LIST
	"""
	node1_idx = COMPLETE_INPUT_TAXA_LIST.index(node1.taxon.label)
	"""
	index of the node1 taxon with respect to COMPLETE_INPUT_TAXA_LIST
	"""
	node2_idx = COMPLETE_INPUT_TAXA_LIST.index(node2.taxon.label)
	"""
	a couplet is referred by the indices of the taxa, sorted in ascending order
	"""
	if (node1_idx < node2_idx):
		target_key = (node1_idx, node2_idx)
		target_reln_type = reln_type
	else:
		target_key = (node2_idx, node1_idx)
		target_reln_type = Complementary_Reln(reln_type)
		
	"""
	check if the key exists - else first create the key
	"""
	if target_key not in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict.setdefault(target_key, Reln_TaxaPair())

	"""
	THE VARIABLE "intersect_ratio" is used for the weighted frequency value
	"""
	if (wt_taxa_subset == True):
		und_tax_list = TaxaPair_Reln_Dict[target_key]._GetUnderlyingTaxonList()
		intersect_taxa_set = set(und_tax_list) & set(curr_tree_taxa)
		intersect_ratio = (len(intersect_taxa_set) * 1.0) / len(und_tax_list)
	else:
		intersect_ratio = (len(curr_tree_taxa) * 1.0) / len(COMPLETE_INPUT_TAXA_LIST)	#1	#lca_rank

	TaxaPair_Reln_Dict[target_key]._AddSupportingTree()
	TaxaPair_Reln_Dict[target_key]._AddXLVal(xl_val / intersect_ratio)
	TaxaPair_Reln_Dict[target_key]._AddEdgeCount(target_reln_type, tr_idx, intersect_ratio)
	TaxaPair_Reln_Dict[target_key]._AddLevel(sum_of_branch_count)
	return

#--------------------------------------------------------
"""
this function derives couplet relations belonging to one tree
that is provided as an input argument to this function
@parameters:  
	WEIGHT_TAXA_SUBSET: If True, the relation takes care of the 
											set of taxa underlying the LCA node for this couplet
	tr_idx: Index of the input tree (starting from 0 to ntrees - 1)
"""
def DeriveCoupletRelations(Curr_tree, WEIGHT_TAXA_SUBSET, tr_idx):
  
	"""
	taxa set of the current tree, and also the count of taxa
	"""
	curr_tree_taxa = [COMPLETE_INPUT_TAXA_LIST.index(x) for x in Curr_tree.infer_taxa().labels()]
	no_of_taxa = len(curr_tree_taxa)  
	
	"""
	traverse the internal nodes of the tree in postorder fashion
	"""
	for curr_node in Curr_tree.postorder_internal_node_iter():
		"""
		compute the rank associated with this node
		"""
		curr_node_level = curr_node.level()
		curr_node_rank = no_of_taxa - curr_node_level
		
		"""
		compute the excess gene count value associated with this node    
		"""
		if (WEIGHT_TAXA_SUBSET == True):
			xl_val = (len(curr_node.leaf_nodes()) - 2)
		else:
			"""
			normalized value of excess gene count 
			with respect to the number of taxa of the current tree
			"""
			xl_val = ((len(curr_node.leaf_nodes()) - 2) * 1.0) / no_of_taxa
		
		"""
		list the leaf and internal children of the current node
		"""
		curr_node_child_leaf_nodes = []
		curr_node_child_internal_nodes = []
		for x in curr_node.child_nodes():
			if (x.is_leaf() == True):
				curr_node_child_leaf_nodes.append(x)
			else:
				curr_node_child_internal_nodes.append(x)
		
		"""
		pair of leaf nodes will be related by sibling relations
		"""
		if (len(curr_node_child_leaf_nodes) > 1):
			for i in range(len(curr_node_child_leaf_nodes) - 1):
				for j in range(i+1, len(curr_node_child_leaf_nodes)):
					DefineLeafPairReln(xl_val, curr_node_level, curr_node_rank, \
						curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j], \
						RELATION_R3, curr_tree_taxa, WEIGHT_TAXA_SUBSET, tr_idx)
		
		"""
		one leaf node (direct descendant) and another leaf node (under one internal node)
		will be related by ancestor / descendant relations
		"""
		if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
			for p in curr_node_child_leaf_nodes:
				for q in curr_node_child_internal_nodes:
					for r in q.leaf_nodes():
						DefineLeafPairReln(xl_val, curr_node_level, curr_node_rank, p, r, RELATION_R1, \
							curr_tree_taxa, WEIGHT_TAXA_SUBSET, tr_idx)
		
		"""
		finally a pair of leaf nodes which are descendant 
		of internal nodes will be related by RELATION_R4 relation
		"""
		if (len(curr_node_child_internal_nodes) > 1):
			for i in range(len(curr_node_child_internal_nodes) - 1):
				for p in curr_node_child_internal_nodes[i].leaf_nodes():
					for j in range(i+1, len(curr_node_child_internal_nodes)):
						for q in curr_node_child_internal_nodes[j].leaf_nodes():
							DefineLeafPairReln(xl_val, curr_node_level, curr_node_rank, p, q, RELATION_R4, \
								curr_tree_taxa, WEIGHT_TAXA_SUBSET, tr_idx)

#--------------------------------------------------------
"""
auxiliary function to append underlying taxa information for a couplet
"""
def AppendUnderlyingTaxa(node1, node2, taxa_under_curr_node):
	# index of the node1 taxon with respect to COMPLETE_INPUT_TAXA_LIST
	node1_idx = COMPLETE_INPUT_TAXA_LIST.index(node1.taxon.label)
	# index of the node1 taxon with respect to COMPLETE_INPUT_TAXA_LIST
	node2_idx = COMPLETE_INPUT_TAXA_LIST.index(node2.taxon.label)
	"""
	a couplet is referred by the indices of the taxa, sorted in ascending order
	"""
	if (node1_idx < node2_idx):
		target_key = (node1_idx, node2_idx)
	else:
		target_key = (node2_idx, node1_idx)
	"""
	check if the key exists - else first create the key
	"""
	if target_key not in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict.setdefault(target_key, Reln_TaxaPair())
	"""
	now add the underlying taxa set 
	"""
	TaxaPair_Reln_Dict[target_key]._AppendUnderlyingTaxonList(taxa_under_curr_node)
	return

#--------------------------------------------------------
"""
For a particular couplet, this function checks all the taxa 
belonging under its LCA node 
for all of the input trees supporting this couplet
"""
def FindCoupletUnderlyingTaxon(Curr_tree):

	"""
	traverse the internal nodes of the tree in postorder fashion
	check all the couplets whose LCA node is the current internal node 
	"""
	for curr_node in Curr_tree.postorder_internal_node_iter():
		"""
		taxa set belonging under this internal node
		"""
		taxa_under_curr_node = GetTaxaUnderInternalNode(curr_node)

		"""
		list all the leaf and internal nodes under curr_node 
		"""
		curr_node_child_leaf_nodes = []
		curr_node_child_internal_nodes = []
		for x in curr_node.child_nodes():
			if (x.is_leaf() == True):
				curr_node_child_leaf_nodes.append(x)
			else:
				curr_node_child_internal_nodes.append(x)
		
		"""
		pair of leaf nodes under curr_node
		"""
		if (len(curr_node_child_leaf_nodes) > 1):
			for i in range(len(curr_node_child_leaf_nodes) - 1):
				node1 = curr_node_child_leaf_nodes[i]
				for j in range(i+1, len(curr_node_child_leaf_nodes)):
					node2 = curr_node_child_leaf_nodes[j]
					AppendUnderlyingTaxa(node1, node2, taxa_under_curr_node)

		"""
		one leaf node (direct descendant) and another leaf node (under one internal node)
		"""
		if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
			for p in curr_node_child_leaf_nodes:
				for q in curr_node_child_internal_nodes:
					for r in q.leaf_nodes():
						AppendUnderlyingTaxa(p, r, taxa_under_curr_node)
						
		"""
		a pair of leaf nodes which are descendant of internal nodes 
		"""
		if (len(curr_node_child_internal_nodes) > 1):
			for i in range(len(curr_node_child_internal_nodes) - 1):
				for p in curr_node_child_internal_nodes[i].leaf_nodes():
					for j in range(i+1, len(curr_node_child_internal_nodes)):
						for q in curr_node_child_internal_nodes[j].leaf_nodes():
							AppendUnderlyingTaxa(p, q, taxa_under_curr_node)
	return

#--------------------------------------------------------
""" 
this function prints the elements of the queue (which stores the couplet scores 
for individual relations 
"""
def PrintQueueInfo(inp_queue, Output_Text_File):
	fp = open(Output_Text_File, 'a')
	for elem in inp_queue:
		fp.write('\n' + str(elem))
	fp.close()

##-----------------------------------------------------
"""
this function reads the input tree list file
@parameters: 
	ROOTED_TREE - whether the treelist to be read as rooted format
	PRESERVE_UNDERSCORE: whether underscores of the taxa name will be preserved or not
	INPUT_FILE_FORMAT: data is read from the file according to NEWICK or NEXUS format
	INPUT_FILENAME: file containing the input treelist
"""
def Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME):
	Inp_TreeList = dendropy.TreeList.get_from_path(INPUT_FILENAME, schema=INPUT_FILE_FORMAT, \
							preserve_underscores=PRESERVE_UNDERSCORE, \
							default_as_rooted=ROOTED_TREE)

	return Inp_TreeList

#--------------------------------------------------
# this function returns the label of an internal or a leaf node 
# in terms of newick representation
def Node_Label(inp_node):
	return str(inp_node.as_newick_string(suppress_edge_lengths=True))

#-----------------------------------------------------
"""
this function returns the list of taxa underlying the given internal node
@param: curr_node: Input node under which the taxa set will be explored
"""
def GetTaxaUnderInternalNode(curr_node):
	taxa_list_from_curr_internal_node = []
	for n in curr_node.leaf_nodes():
		taxa_list_from_curr_internal_node.append(n.taxon.label)
	return taxa_list_from_curr_internal_node

#----------------------------------------
def Complementary_Reln(inp_reln):
	if (inp_reln == RELATION_R3) or (inp_reln == RELATION_R4):
		return inp_reln
	elif (inp_reln == RELATION_R1):
		return RELATION_R2
	else:
		return RELATION_R1

#-----------------------------------------------------------------
"""
this function returns the list of taxa underlying the given internal node
in preorder traversal
@param: inp_node: Input node under which the taxa set will be explored
				taxa_list: Output taxa list in preorder traversal order
				inp_set_of_taxa: A superset of taxon; the 'taxa_list' should be a subset of it
"""
def GetPreorderTaxaList(inp_node, taxa_list, inp_set_of_taxa):
	for n in inp_node.preorder_iter():
		if (n.is_leaf() == True):
			if n.taxon.label in inp_set_of_taxa:
				taxa_list.append(n.taxon.label)
	
	return taxa_list

#------------------------------------------------
"""
this function computes average XL information between a pair of taxa clusters
@param: taxa_clust1: first taxa list
				taxa_clust2: second taxa list
				DIST_MAT_TYPE: Type of distance employed
				single_elem: can contain one of possible three values
				0: only one element of taxa_clust1 and one element of taxa_clust2 will be compared
				1: cluster containing taxa_clust1[0] and cluster containing taxa_clust2[0] will be compared
				2: All pairs of elements of taxa_clust1 and taxa_clust2 will be compared
"""
def FindAvgXL(taxa_clust1, taxa_clust2, DIST_MAT_TYPE, single_elem=2, type_of_output=0):
	"""
	if single_elem = 0
	we compare taxa_clust1[0] and taxa_clust2[0], in terms of the preorder level
	
	if single_elem = 1
	we check the first preorder level taxon of both lists taxa_clust1 and taxa_clust2
	suppose the taxon names are taxa1 and taxa2
	but instead of comparing taxa1 and taxa2 only
	we compare the original taxa clusters (may have cardinality > 1) containing taxa1 and taxa2
	
	if single_elem = 2
	we compare pairwise all the elements belonging to taxa_clust1 and taxa_clust2
	"""
	if (single_elem == 1):
		taxa1 = taxa_clust1[0]
		taxa2 = taxa_clust2[0]
		clust1 = Taxa_Info_Dict[taxa1]._Get_Taxa_Part_Clust_Idx()
		clust2 = Taxa_Info_Dict[taxa2]._Get_Taxa_Part_Clust_Idx()
		taxa_list1 = Cluster_Info_Dict[clust1]._GetSpeciesList()
		taxa_list2 = Cluster_Info_Dict[clust2]._GetSpeciesList()
	elif (single_elem == 2):
		taxa_list1 = taxa_clust1
		taxa_list2 = taxa_clust2 
	else:
		taxa_list1 = []
		taxa_list1.append(taxa_clust1[0])
		taxa_list2 = []
		taxa_list2.append(taxa_clust2[0])
		
	curr_taxa_pair_list = []
	for x1 in taxa_list1:
		for x2 in taxa_list2:  
			key1 = (x1, x2)
			key2 = (x2, x1)
			#print 'key1: ', key1, ' key2: ', key2
			if key1 in TaxaPair_Reln_Dict:
				val = TaxaPair_Reln_Dict[key1]._GetNormalizedXLSumGeneTrees(DIST_MAT_TYPE)
				curr_taxa_pair_list.append(val)
			elif key2 in TaxaPair_Reln_Dict:
				val = TaxaPair_Reln_Dict[key2]._GetNormalizedXLSumGeneTrees(DIST_MAT_TYPE)
				curr_taxa_pair_list.append(val)
	
	# average of this pairwise list is used as the XL approximation
	if (len(curr_taxa_pair_list) > 0):
		if (type_of_output == 0):
			return (sum(curr_taxa_pair_list) * 1.0) / len(curr_taxa_pair_list)
		else:
			return max(curr_taxa_pair_list)
			#return min(curr_taxa_pair_list)
	else:
		return 0

#------------------------------------------------
"""
this function computes average XL information between a pair of taxa clusters
@param: taxa_clust1: first taxa list
				taxa_clust2: second taxa list
				DIST_MAT_TYPE: Type of distance employed
				single_elem: can contain one of possible three values
				0: only one element of taxa_clust1 and one element of taxa_clust2 will be compared
				1: cluster containing taxa_clust1[0] and cluster containing taxa_clust2[0] will be compared
				2: All pairs of elements of taxa_clust1 and taxa_clust2 will be compared
"""
def FindAvgInternodeCount(taxa_clust1, taxa_clust2, single_elem=2, type_of_output=1):
	"""
	if single_elem = 0
	we compare taxa_clust1[0] and taxa_clust2[0], in terms of the preorder level
	
	if single_elem = 1
	we check the first preorder level taxon of both lists taxa_clust1 and taxa_clust2
	suppose the taxon names are taxa1 and taxa2
	but instead of comparing taxa1 and taxa2 only
	we compare the original taxa clusters (may have cardinality > 1) containing taxa1 and taxa2
	
	if single_elem = 2
	we compare pairwise all the elements belonging to taxa_clust1 and taxa_clust2
	"""
	if (single_elem == 1):
		taxa1 = taxa_clust1[0]
		taxa2 = taxa_clust2[0]
		clust1 = Taxa_Info_Dict[taxa1]._Get_Taxa_Part_Clust_Idx()
		clust2 = Taxa_Info_Dict[taxa2]._Get_Taxa_Part_Clust_Idx()
		taxa_list1 = Cluster_Info_Dict[clust1]._GetSpeciesList()
		taxa_list2 = Cluster_Info_Dict[clust2]._GetSpeciesList()
	elif (single_elem == 2):
		taxa_list1 = taxa_clust1
		taxa_list2 = taxa_clust2 
	else:
		taxa_list1 = []
		taxa_list1.append(taxa_clust1[0])
		taxa_list2 = []
		taxa_list2.append(taxa_clust2[0])
		
	curr_taxa_pair_list = []
	for x1 in taxa_list1:
		for x2 in taxa_list2:  
			key1 = (x1, x2)
			key2 = (x2, x1)
			#print 'key1: ', key1, ' key2: ', key2
			if key1 in TaxaPair_Reln_Dict:
				val = TaxaPair_Reln_Dict[key1]._GetAvgSumLevel()
				curr_taxa_pair_list.append(val)
			elif key2 in TaxaPair_Reln_Dict:
				val = TaxaPair_Reln_Dict[key2]._GetAvgSumLevel()
				curr_taxa_pair_list.append(val)
	
	# average of this pairwise list is used as the XL approximation
	if (len(curr_taxa_pair_list) > 0):
		if (type_of_output == 0):
			return (sum(curr_taxa_pair_list) * 1.0) / len(curr_taxa_pair_list)
		else:
			return max(curr_taxa_pair_list)
			#return min(curr_taxa_pair_list)
	else:
		return 0
