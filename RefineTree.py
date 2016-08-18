#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

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
					' curr_taxa_idx: ' + str(curr_taxa_idx) + ' curr_taxa: ' + \
						str(curr_taxa) + ' supporting tree set: ' + str(curr_taxa_support_tree)) 
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
	
	if (DEBUG_LEVEL >= 2):
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
def Refine_Supertree_Binary_Form(Input_Treelist, Curr_Suptree, Output_Text_File):

	# contains all the taxa list of all multifurcating nodes
	Global_Clust_Species_List = []
	
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
			"""
			append the clust_species_list in the list "Global_Clust_Species_List"
			"""
			Global_Clust_Species_List.append(clust_species_list)
	
	"""
	now navigate through individual elements of Global_Clust_Species_List
	"""
	for nnode in range(len(Global_Clust_Species_List)):
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n\n *** Examining the index ' + str(nnode) + '  of Global_Clust_Species_List')
			fp.write('\n The set of taxa for analysis: ' + str(Global_Clust_Species_List[nnode]))
			fp.write('\n Degree of multifurcation: ' + str(len(Global_Clust_Species_List[nnode])))
			fp.close()
		
		clust_species_list = Global_Clust_Species_List[nnode]
		taxa_list = []
		for i in range(len(clust_species_list)):
			taxa_list.extend(clust_species_list[i])

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n taxa_list: ' + str(taxa_list))
			fp.close()
		
		LCA_node = Curr_Suptree.mrca(taxon_labels=taxa_list)

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n Label of LCA_node: ' + str(Node_Label(LCA_node)))
			fp.close()
			
		Curr_Suptree = ResolveMultifurcation_Latest(LCA_node, Input_Treelist, Curr_Suptree, \
			Global_Clust_Species_List[nnode], len(Global_Clust_Species_List[nnode]), Output_Text_File, nnode)
		
	"""
	now delete all the files within the current directory which ends with '_nexus.tre' and 'temp.txt'
	"""
	for f in os.listdir("."):
		if f.endswith("_nexus.tre"):
			sys_cmd = 'rm ' + str(f)
			os.system(sys_cmd)

	return Curr_Suptree

