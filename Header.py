#!/usr/bin/env python

import dendropy
from dendropy import Tree, Taxon, TaxonSet, Node
import numpy
import time
import os
from cStringIO import StringIO
from optparse import OptionParser

# we define custom edge types
RELATION_R3 = 0	# equality relationship
RELATION_R1 = 1
RELATION_R2 = 2
RELATION_R4 = 3	# no relationship
UNDEFINED_RELATION = 4

# variables used to denote whether we use traditional NJ method
# or use a variant of it, namely the agglomerative clustering
TRADITIONAL_NJ = 1
AGGLO_CLUST = 2

#---------------------------------------------
""" 
this is a dictionary for storing information about individual taxa clusters
each cluster is basically a collection of taxa related via relation r3
"""
Cluster_Info_Dict = dict()

""" 
the dictionary defines one particular taxa and its associated information
"""
Taxa_Info_Dict = dict()

""" 
this dictionary defines the taxa pair (couplet) relations and associated operations
each entry of this dictionary is indexed by a pair of taxon labels 
"""
TaxaPair_Reln_Dict = dict()

"""
queue storing relations of conflicting couplets (supporting more than one type of relation 
in the input trees)
"""
Cost_List_Taxa_Pair_Multi_Reln = []

"""
queue storing relations of non-conflicting couplets (supporting only one type of relation 
in the input trees)
"""
Cost_List_Taxa_Pair_Single_Reln = [] 

""" 
this list contains the complete set of taxa present in the input trees 
"""
COMPLETE_INPUT_TAXA_LIST = []

""" 
this list contains the current set of active taxa cluster (indices)
"""
CURRENT_CLUST_IDX_LIST = []

"""
this is the debug level
set for printing the necessary information
"""
DEBUG_LEVEL = 2

# add - sourya
MODE_PERCENT = 0.25	#0.4 #0.35	#0.5 
MODE_BIN_COUNT = 40

##-----------------------------------------------------
"""
this class defines a leaf node of input candidate source trees
that is, it corresponds to one taxa 
"""
class Single_Taxa(object):
	def __init__(self):
		""" 
		this variable signifies the cluster index that the taxa belongs 
		if the value is -1 (< 0) then it is still not inserted in a cluster 
		otherwise it is part of a valid cluster 
		"""
		self.clust_idx_part = -1
		"""
		these lists are constructed using the original allowed relations between individual couplets
		when individual gene trees are analyzed
		"""
		self.Allowed_Reln_List = [[] for i in range(4)]

	def _Get_Taxa_Part_Clust_Idx(self):
		return self.clust_idx_part
		
	def _Set_Clust_Idx_taxa_Part(self, inp_clust_idx):
		self.clust_idx_part = inp_clust_idx
	
	"""
	appends one taxa for one target relation (which is allowed)
	"""
	def _AddAllowedReln(self, reln_type, inp_taxa):
		self.Allowed_Reln_List[reln_type].append(inp_taxa)
		
	def _GetAllowedRelnSet(self, reln_type):
		return self.Allowed_Reln_List[reln_type]
	
	# this function is called after formation of consensus tree
	def _PrintFinalTaxaInfo(self, key, Output_Text_File):
		fp = open(Output_Text_File, 'a')    
		fp.write('\n taxa key: ' + str(key))
		fp.write('\n taxa label: ' + str(COMPLETE_INPUT_TAXA_LIST[key]))
		fp.write('\n taxa is part of the cluster ID: ' + str(self.clust_idx_part))
		fp.close()

##-----------------------------------------------------
""" 
this class defines a couplet, according to the information obtained from input trees
key of this class --- taxa1, taxa2  
the consensus relation, frequencies, priority, and support scores 
in the class, the edge type signifies the relation between a pair of taxa
"""
class Reln_TaxaPair(object):
	def __init__(self):    
		""" 
		frequencies of individual relations 
		there are 4 types of edges (relationship) between a pair of taxa 
		"""
		self.freq_count = [0] * 4    
		""" 
		a connection priority value is defined as the 
		no of occurrences of this particular relation between this pair of taxa 
		minus the sum of no of occurrences of other relation types between this couplet
		"""
		self.priority_reln = [0] * 4    
		""" 
		this is the support score for different types of relations between a couplet
		"""
		self.support_score = [0] * 4
		""" 
		For this couplet, it stores the extra gene count with respect to all the gene trees
		"""
		self.XL_sum_gene_trees = []
		""" 
		this variable stores the no of trees supporting the taxa pair 
		"""
		self.supporting_trees = 0
		"""
		this list contains the union of taxa list underlying the LCA of this couplet
		for individual input trees
		"""
		self.LCA_Underlying_Taxa_List = []
		"""
		this is a variable containing the binned average of the XL values
		of very high frequency
		initially the value is set as -1, to signify that the computation is not done
		once the computation (for a couplet) is done, the value is subsequently used and returned
		"""
		self.binned_avg_XL = -1
		"""
		this variable contains the sum of tree internodees between individual couplets
		computed for all the gene trees supported by this couplet
		"""
		self.sum_internode_count = 0
		#""" 
		#frequencies of individual relations 
		#there are 4 types of edges (relationship) between a pair of taxa 
		#"""
		#self.freq_instance = [0] * 4
		#"""
		#this is a collection of 4 different lists of variable size
		#individual lists represent R1 to R4 relation
		#contents of a list are the tree indices supporting the occurrence of one particular relation
		#*** Note: initialization of the list with empty sublists - this is the standard syntax ****
		#"""
		#self.Tree_Reln_List = [[] for i in range(4)]
		"""
		allowed relation list which are included in the support score queue
		"""
		self.Allowed_Reln_List = []
	
	#----------------------------------
	def _AddAllowedReln(self, inp_reln):
		self.Allowed_Reln_List.append(inp_reln)
	
	def _Check_Single_Reln_OnlyAllowed(self, inp_reln):
		if (len(self.Allowed_Reln_List) == 1) and (inp_reln in self.Allowed_Reln_List):
			return True
		return False
	
	def _GetAllowedRelnList(self):
		return self.Allowed_Reln_List
	#----------------------------------
	def _AddLevel(self, val):
		self.sum_internode_count = self.sum_internode_count + val

	def _GetAvgSumLevel(self):
		return (self.sum_internode_count * 1.0) / self.supporting_trees
	#----------------------------------
	"""
	Appends underlying taxon set of the LCA node for this couplet
	Note: we store the index of individual taxon, with respect to the COMPLETE_INPUT_TAXA_LIST
	"""
	def _AppendUnderlyingTaxonList(self, inp_list):
		inp_taxa_list_index = [COMPLETE_INPUT_TAXA_LIST.index(x) for x in inp_list]
		self.LCA_Underlying_Taxa_List = list(set(self.LCA_Underlying_Taxa_List) | set(inp_taxa_list_index))

	"""
	Get the complete set of taxon, underlying the LCA nodes, for all input trees, 
	corresponding to this couplet
	"""
	def _GetUnderlyingTaxonList(self):
		return self.LCA_Underlying_Taxa_List

	#----------------------------------
	"""
	this function adds one supporting tree for this couplet
	"""
	def _AddSupportingTree(self):
		self.supporting_trees = self.supporting_trees + 1

	"""
	this function returns the number of input trees supporting this couplet
	"""
	def _GetNoSupportTrees(self):
		return self.supporting_trees
	
	#----------------------------------
	"""
	this function adds one XL value, computed for a particular input tree
	to the list of XL values for this couplet
	"""
	def _AddXLVal(self, val):
		self.XL_sum_gene_trees.append(val)
	
	def _GetXLList(self):
		return self.XL_sum_gene_trees
	
	def _GetXLSumGeneTrees(self):
		return sum(self.XL_sum_gene_trees)
		
	"""
	this function computes the average of XL measures
	"""
	def _GetAvgXLGeneTrees(self):
		return (sum(self.XL_sum_gene_trees) * 1.0) / self.supporting_trees

	"""
	function to return the average of XL values for this couplet
	depending on the user parameters, average, median, or binned average XL is returned
	"""
	def _GetNormalizedXLSumGeneTrees(self, dist_type):
		if (dist_type == 1) or (dist_type == 3):
			return self._GetAvgXLGeneTrees()
		elif (dist_type == 2):
			# average of mean and mode
			return (self._GetAvgXLGeneTrees() + self._GetMultiModeXLVal()) / 2.0

	#------------------------------------------------
	"""
	this function computes the binned average of XL values associated for this couplet
	"""
	def _GetMultiModeXLVal(self, Output_Text_File=None):
		if (self.binned_avg_XL == -1):
			
			Bin_Width = (1.0 / MODE_BIN_COUNT)
			len_list = [0] * MODE_BIN_COUNT
			
			if Output_Text_File is not None:
				fp = open(Output_Text_File, 'a') 
			
			# sort the XL list
			self.XL_sum_gene_trees.sort()
			
			for j in range(len(self.XL_sum_gene_trees)):
				curr_xl_val = self.XL_sum_gene_trees[j]
				bin_idx = int(curr_xl_val / Bin_Width)
				if (bin_idx == MODE_BIN_COUNT):
					bin_idx = bin_idx - 1
				len_list[bin_idx] = len_list[bin_idx] + 1
			
			if Output_Text_File is not None:
				for i in range(MODE_BIN_COUNT):
					fp.write('\n bin idx: ' + str(i) + ' len:  ' + str(len_list[i]))
			
			# this is the maximum length of a particular bin
			# corresponding to max frequency
			max_freq = max(len_list)
			
			if Output_Text_File is not None:
				fp.write('\n Max freq: ' + str(max_freq))
			
			num = 0
			denom = 0
			for i in range(MODE_BIN_COUNT):
				if (len_list[i] >= (MODE_PERCENT * max_freq)):
					list_start_idx = sum(len_list[:i])
					list_end_idx = list_start_idx + len_list[i] - 1
					value_sum = sum(self.XL_sum_gene_trees[list_start_idx:(list_end_idx+1)])
					num = num + value_sum
					denom = denom + len_list[i]
					if Output_Text_File is not None:
						fp.write('\n Included bin idx: ' + str(i) + ' starting point: ' + str(list_start_idx) \
							+ 'ending point: ' + str(list_end_idx) + ' sum: ' + str(value_sum))
			
			self.binned_avg_XL = (num / denom)
			
			if Output_Text_File is not None:
				fp.write('\n Final binned average XL: ' + str(self.binned_avg_XL))
				fp.close()
			
		return self.binned_avg_XL
	
	#----------------------------------
	"""
	this function returns the frequency of the consensus relation
	"""
	def _GetConsensusFreq(self):
		return max(self.freq_count)

	"""
	this function checks whether the input relation type is a consensus relation
	among this couplet
	"""
	def _CheckTargetRelnConsensus(self, inp_reln_type):
		if (self.freq_count[inp_reln_type] == max(self.freq_count)):
			return True
		return False

	"""
	this function checks whether the input relation type is a consensus relation
	among this couplet
	"""
	def _CheckTargetRelnMajorityConsensus(self, inp_reln_type):
		if (self.freq_count[inp_reln_type] >= (0.5 * sum(self.freq_count))):
			return True
		return False

	"""
	this function returns the frequency of the input relation
	specified by the variable 'reln_type'
	"""
	def _GetEdgeWeight(self, reln_type):
		return self.freq_count[reln_type]      
	
	"""
	this function adds a specified frequency count (default 1)
	corresponding to the relation specified by the variable 'reln_type'
	"""
	def _AddEdgeCount(self, reln_type, tr_idx, val=1):
		self.freq_count[reln_type] = self.freq_count[reln_type] + val
		## add - sourya
		#self.freq_instance[reln_type] = self.freq_instance[reln_type] + 1
		## add - sourya
		#self.Tree_Reln_List[reln_type].append(tr_idx)
	
	#def _NormalizeFreq(self):
		#for reln_type in range(4):
			#if (self.freq_count[reln_type] > 0):
				#self.freq_count[reln_type] = (self.freq_count[reln_type] * 1.0) / self.freq_instance[reln_type]
	
	#----------------------------------
	""" 
	this function computes the support score value associated with individual couplet
	for all different relations
	"""
	def _SetCostMetric(self):
		for reln_type in range(4):
			# assign the score metric for this relation type
			self.support_score[reln_type] = self.freq_count[reln_type] * self.priority_reln[reln_type]

	"""
	this function returns the support score of the input relation
	specified by the variable 'reln_type'
	"""
	def _GetEdgeCost_ConnReln(self, reln_type):
		return self.support_score[reln_type]
	
	"""
	this function updates (increments) the support score of the input relation
	specified by the variable 'reln_type'
	and by the amount 'incr_cost'
	"""
	def _IncrEdgeCost_ConnReln(self, reln_type, incr_cost):
		self.support_score[reln_type] = self.support_score[reln_type] + incr_cost

	#----------------------------------	
	""" 
	this function calculates connection priority value for 
	each of the relation types
	"""
	def _SetConnPrVal(self):
		"""
		this is the sum of frequencies for all the relation types
		"""
		listsum = sum(self.freq_count)
		"""
		now determine the connection priority of a 
		particular relation type with respect to other relations     
		"""
		for reln_type in range(4):
			"""
			here we use the difference of current relation type 
			frequency with the frequencies of all other relations
			"""
			self.priority_reln[reln_type] = 2 * self.freq_count[reln_type] - listsum

	"""
	this function returns the priority value for a given input relation 
	specified by the variable 'reln_type'
	"""
	def _GetConnPrVal(self, reln_type):
		return self.priority_reln[reln_type]
	
	"""
	this function checks whether a couplet is non-conflicting
	that is, only one relation between them exists throughout all the gene trees
	in such a case, a binary variable 1 and the corresponding 
	relation type is returned in the form of a list
	otherwise, a value 0 and a defaukt relation R4 is returned
	"""
	def _CheckNonConflictingCouplet(self):
		# this is the sum of frequencies for all the relation types
		listsum = sum(self.freq_count)
		""" 
		this code section is used when there exists an unique relation (non-conflicting couplet)
		and we try to detect it
		"""
		outlist = [0, RELATION_R4]
		for reln_type in range(4):
			if (self.freq_count[reln_type] == listsum) and (listsum > 0):
				outlist = [1, reln_type]
				break
			elif (self.freq_count[reln_type] > 0) and (self.freq_count[reln_type] < listsum):
				break
		return outlist
	
	#----------------------------------	
	"""
	this function prints all the information associated with a couplet
	"""
	def _PrintRelnInfo(self, key, Output_Text_File):
		fp = open(Output_Text_File, 'a')    
		fp.write('\n\n\n taxa pair key: ' + str(key) + ' couplet:  ' + str(COMPLETE_INPUT_TAXA_LIST[key[0]]) + '  and ' + str(COMPLETE_INPUT_TAXA_LIST[key[1]]))
		fp.write('\n relations [type/count/priority_reln/score]: ')
		for i in range(4):
			fp.write('\n [' + str(i) + '/' + str(self.freq_count[i]) + '/' + str(self.priority_reln[i]) + '/' + str(self.support_score[i]) + ']')
		fp.write('\n AVERAGE Sum of excess gene **** : ' + str(self._GetAvgXLGeneTrees()))
		fp.write('\n No of supporting trees : ' + str(self.supporting_trees))
		fp.write('\n Average sum of internode count : ' + str(self._GetAvgSumLevel()))    
		#fp.write('\n R1 relation occurring trees: ' + str(self.Tree_Reln_List[RELATION_R1]))
		#fp.write('\n R2 relation occurring trees: ' + str(self.Tree_Reln_List[RELATION_R2]))
		#fp.write('\n R3 relation occurring trees: ' + str(self.Tree_Reln_List[RELATION_R3]))
		#fp.write('\n R4 relation occurring trees: ' + str(self.Tree_Reln_List[RELATION_R4]))
		fp.close()
	
#-----------------------------------------------------
""" 
this class is representative of a cluster of taxa that are related via equality relationship 
according to the rule of equivalence partition 
"""
class Cluster_node(object):
	def __init__(self, inp_taxa=None):
		"""
		taxa list of the current cluster
		"""
		self.Species_List = [] 
		"""
		# set to 1 once the cluster is traversed during DFS order of traversing the clusters
		this is required in printing the supertree in newick format 
		"""
		self.explored = 0   
		"""
		stores the indices of clusters cy, such that, depending on the relation type:
		curr_clust->cy /  cy->curr_clust / R3 (cy, curr_clust) / R4 (cy, curr_clust) are present
		"""
		self.Reln_List = [[] for i in range(4)]
		"""
		during initialization, append one tuple to this cluster
		"""
		if inp_taxa is not None:
			self._Append_taxa(inp_taxa)    

	#def _ResetItems(self):
		#self.explored = 0
		#self.out_edge_list[:] = []
		#self.out_edge_list = []
		#self.in_edge_list[:] = []
		#self.in_edge_list = []
		#self.no_edge_list[:] = []
		#self.no_edge_list = []
		#self.eq_edge_list[:] = []
		#self.eq_edge_list = []
	
	#---------------------------------------------
	"""
	these functions keep track whether the cluster node is used 
	during newick string formation for supertree construction
	each of the clusters (containing a set of taxa) should be visited 
	exactly once for supertree generation
	"""
	def _SetExploredStatus(self):
		self.explored = 1

	def _ResetExploredStatus(self):
		self.explored = 0
		
	def _GetExploredStatus(self):
		return self.explored

	#---------------------------------------------
	"""
	following function modifies and returns (if required) 
	the constituent taxa list within this cluster
	"""
	def _GetSpeciesList(self):
		return self.Species_List
	
	def _GetCardinality(self):
		return len(self.Species_List)
				
	# append one species information in this cluster
	def _Append_taxa(self, inp_taxa):
		if inp_taxa not in self.Species_List:
			self.Species_List.append(inp_taxa)

	def _Remove_taxa(self, inp_taxa):
		if inp_taxa in self.Species_List:
			self.Species_List.remove(inp_taxa)    

	#---------------------------------------------
	"""
	returns the number of clusters cy such that cy->curr_clust connection is present
	"""
	def _Get_Indegree(self):
		#return len(self.in_edge_list)
		return len(self.Reln_List[RELATION_R2])

	"""
	returns the number of clusters cy such that curr_clust->cy connection is present
	"""
	def _Get_Outdegree(self):
		#return len(self.out_edge_list)
		return len(self.Reln_List[RELATION_R1])
	
	"""
	returns the list of cluster indices connected with the current cluster
	via the relation "reln_type"
	"""
	def _GetClustRelnList(self, reln_type):
		return self.Reln_List[reln_type]
	
	"""
	appends one cluster index to the list of the specified "reln_type"
	"""
	def _AddRelnInstance(self, reln_type, dest_clust_idx):
		if (dest_clust_idx not in self.Reln_List[reln_type]):
			self.Reln_List[reln_type].append(dest_clust_idx)
	
	"""
	removes one cluster index to the list of the specified "reln_type"
	"""
	def _RemoveRelnInstance(self, reln_type, dest_clust_idx):
		if dest_clust_idx in self.Reln_List[reln_type]:
			self.Reln_List[reln_type].remove(dest_clust_idx)    

	#--------------------------------------------------------
	def _PrintClusterInfo(self, key, Output_Text_File):
		fp = open(Output_Text_File, 'a')    
		fp.write('\n cluster key: ' + str(key))
		fp.write('\n species list: ' + str(self.Species_List))
		fp.write('\n out edge list (R1): ' + str(self.Reln_List[RELATION_R1]))
		fp.write('\n in edge list (R2): ' + str(self.Reln_List[RELATION_R2]))
		fp.write('\n No edge list (R4): ' + str(self.Reln_List[RELATION_R4]))
		fp.write('\n Eq edge list (R3): ' + str(self.Reln_List[RELATION_R3]))
		fp.close()
