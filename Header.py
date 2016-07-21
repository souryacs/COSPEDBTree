#!/usr/bin/env python

import dendropy
from dendropy import TreeList, Tree, Taxon, TaxonSet, Node
import numpy
import time
import os
from cStringIO import StringIO
from optparse import OptionParser
import copy
import networkx as nx
#import matplotlib.pyplot as plt

"""
following variables denote different relations between a couplet
"""
"""
for a couplet (x,y), R3(x,y) means x and y are siblings
"""
RELATION_R3 = 0	
"""
for a couplet (x,y), R1(x,y) means LCA(x,y) = parent node of x
"""
RELATION_R1 = 1
"""
for a couplet (x,y), R2(x,y) means R1(y,x)
"""
RELATION_R2 = 2
"""
for a couplet (x,y), R4(x,y) means x and y have independent evolution
"""
RELATION_R4 = 3
UNDEFINED_RELATION = 4

"""
Two types of NJ methods are explored:
1) traditional NJ method with difference of distance values
2) Normalized / ratio of distance values
"""
TRADITIONAL_NJ = 1
AGGLO_CLUST = 2

#---------------------------------------------
""" 
this is a dictionary for storing information about individual taxa clusters
each cluster is basically a collection of taxa related via relation r3
"""
Cluster_Info_Dict = dict()

"""
this is a dictionary containing information of a pair of taxa clusters
"""
Cluster_Pair_Info_Dict = dict()

""" 
the dictionary defines one particular taxa and its associated information
"""
Taxa_Info_Dict = dict()

""" 
this dictionary defines the taxa pair (couplet) relations and associated operations
each entry of this dictionary is indexed by a pair of taxon labels 
"""
TaxaPair_Reln_Dict = dict()

#--------------------------------------------------
"""
queue storing R3 relations of couplets
provided that the couplet is non conflicting (only R3 relation is present)
"""
Queue_Score_R3_SingleReln = []

"""
queue storing R3 relations of couplets
provided that the couplet is conflicting but the R3 relation is majority consensus
"""
Queue_Score_R3_MajCons = []

"""
queue storing relations of individual cluster pairs
"""
Queue_Score_Cluster_Pair = []

"""
queue storing relations of individual cluster pairs 
only when one type of relation is allowed between them
"""
Queue_Score_Cluster_Pair_NonConflict = []

#--------------------------------------------------
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
DEBUG_LEVEL = 0

"""
this is the percentage of consensus relation frequency
that a relation must possess to be included as a candidate relation between a cluster pair
"""
PERCENT_MAX_FREQ1 = 0.5
PERCENT_MAX_FREQ2 = 0.7

#-----------------------------------------------------
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
		this list stores the indices of input trees having this taxon
		"""
		self.Support_Tree_List = []

	def _AddSupportTree(self, idx):
		self.Support_Tree_List.append(idx)
		
	def _GetSupportTreeList(self):
		return self.Support_Tree_List

	def _Get_Taxa_Part_Clust_Idx(self):
		return self.clust_idx_part
		
	def _Set_Clust_Idx_taxa_Part(self, inp_clust_idx):
		self.clust_idx_part = inp_clust_idx
	
	# this function is called after formation of consensus tree
	def _PrintFinalTaxaInfo(self, key, Output_Text_File):
		fp = open(Output_Text_File, 'a')    
		fp.write('\n taxa key: ' + str(key))
		fp.write('\n taxa label: ' + str(COMPLETE_INPUT_TAXA_LIST[key]))
		fp.write('\n taxa is part of the cluster ID: ' + str(self.clust_idx_part))
		fp.close()

#-----------------------------------------------------
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
		this variable contains the sum of tree internodees between individual couplets
		computed for all the gene trees supported by this couplet
		"""
		self.sum_internode_count = 0
		"""
		allowed relation list which are included in the support score queue
		"""
		self.Allowed_Reln_List = []
	
	#----------------------------------
	def _AddAllowedReln(self, inp_reln):
		if inp_reln not in self.Allowed_Reln_List:
			self.Allowed_Reln_List.append(inp_reln)
	
	def _Check_Single_Reln_OnlyAllowed(self, inp_reln):
		if (len(self.Allowed_Reln_List) == 1) and (inp_reln in self.Allowed_Reln_List):
			return True
		return False
	
	def _GetAllowedRelnList(self):
		return self.Allowed_Reln_List
	
	def _CheckInpRelnAllowed(self, inp_reln):
		if inp_reln in self.Allowed_Reln_List:
			return True
		return False
	
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
		# comment - sourya
		#inp_taxa_list_index = [COMPLETE_INPUT_TAXA_LIST.index(x) for x in inp_list]
		#self.LCA_Underlying_Taxa_List = list(set(self.LCA_Underlying_Taxa_List) | set(inp_taxa_list_index))
		# add - sourya
		self.LCA_Underlying_Taxa_List = list(set(self.LCA_Underlying_Taxa_List) | set(inp_list))

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
	def _AddEdgeCount(self, reln_type, val=1):
		self.freq_count[reln_type] = self.freq_count[reln_type] + val
	
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
	def _PrintRelnInfo(self, key, Output_Text_File, WEIGHT_TAXA_SUBSET):
		fp = open(Output_Text_File, 'a')    
		fp.write('\n\n\n taxa pair key: ' + str(key) + ' couplet:  ' + str(COMPLETE_INPUT_TAXA_LIST[key[0]]) \
			+ ', ' + str(COMPLETE_INPUT_TAXA_LIST[key[1]]))
		fp.write('\n relations [type/count]: ')
		for i in range(4):
			fp.write('\t [' + str(i) + '/' + str(self.freq_count[i]) + ']')
		fp.write('\n AVERAGE Sum of excess gene **** : ' + str(self._GetAvgXLGeneTrees()))
		fp.write('\n No of supporting trees : ' + str(self.supporting_trees))
		fp.write('\n Average sum of internode count : ' + str(self._GetAvgSumLevel()))    
		fp.write('\n Allowed reln list: ' + str(self.Allowed_Reln_List))
		if (WEIGHT_TAXA_SUBSET == True):
			fp.write('\n Underlying taxa indices list: ' + str(self.LCA_Underlying_Taxa_List))
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
		set to 1 once the cluster is traversed during DFS order of traversing the clusters
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
	def _PrintClusterInfo(self, key, Output_Text_File, suppress=False):
		fp = open(Output_Text_File, 'a')    
		fp.write('\n\n cluster key: ' + str(key))
		fp.write('\n species list: ' + str(self.Species_List))
		if (suppress == False):
			fp.write('\n out edge list (R1): ' + str(self.Reln_List[RELATION_R1]))
			fp.write('\n in edge list (R2): ' + str(self.Reln_List[RELATION_R2]))
		fp.close()

#-----------------------------------------------------
""" 
this class is representative of a pair of taxa clusters
"""
class Cluster_Pair(object):
	def __init__(self, R1_freq, R2_freq, R3_freq, R4_freq):
		"""
		frequencies of individual relations
		"""
		self.freq_count = [0] * 4
		"""
		priorities of individual relations
		"""
		self.priority_measure = [0] * 4
		"""
		support scores of individual relations
		"""
		self.support_score_measure = [0] * 4
		"""
		set of relations that can exist between this pair of clusters
		in the created DAG
		"""
		self.possible_reln_list = []
		"""
		initialize the frequencies and the pseudo R4 relation frequencies for this pair of cluster
		"""
		self._InitRelnFreq(R1_freq, R2_freq, R3_freq, R4_freq)
		"""
		initialize the priority and support score measures for individual relations
		"""
		self._ComputePriority()
		self._ComputeSupportScore()
	
	#-----------------------------------------------------------
	def _ComputePriority(self):
		self.priority_measure[RELATION_R1] = self.freq_count[RELATION_R1] + self.freq_count[RELATION_R3] \
			- self.freq_count[RELATION_R2] - self.freq_count[RELATION_R4]
		self.priority_measure[RELATION_R2] = self.freq_count[RELATION_R2] + self.freq_count[RELATION_R3] \
			- self.freq_count[RELATION_R1] - self.freq_count[RELATION_R4]
		self.priority_measure[RELATION_R4] = self.freq_count[RELATION_R4] \
			- self.freq_count[RELATION_R1] - self.freq_count[RELATION_R2]
		self.priority_measure[RELATION_R3] = self.freq_count[RELATION_R3] - self.freq_count[RELATION_R1] \
			- self.freq_count[RELATION_R2] - self.freq_count[RELATION_R4]
	
	def _ComputeSupportScore(self):
		for r in range(4):
			# modify - sourya
			#self.support_score_measure[r] = self.freq_count[r] * self.priority_measure[r]
			self.support_score_measure[r] = self.freq_count[r] + self.priority_measure[r]
	
	def _GetPriority(self, inp_reln):
		return self.priority_measure[inp_reln]
	
	def _GetSupportScore(self, inp_reln):
		return self.support_score_measure[inp_reln]
	
	def _InitRelnFreq(self, R1_freq, R2_freq, R3_freq, R4_freq):
		self.freq_count[RELATION_R1] = R1_freq
		self.freq_count[RELATION_R2] = R2_freq
		self.freq_count[RELATION_R3] = R3_freq
		self.freq_count[RELATION_R4] = R4_freq
	
	def _AddPossibleReln(self, inp_reln):
		if inp_reln not in self.possible_reln_list:
			self.possible_reln_list.append(inp_reln)
	
	def _GetPossibleRelnList(self):
		return self.possible_reln_list
	
	def _RemovePossibleReln(self, inp_reln):
		if inp_reln in self.possible_reln_list:
			self.possible_reln_list.remove(inp_reln)
	
	def _GetFreq(self, inp_reln):
		return self.freq_count[inp_reln]
		
	def _PrintClusterPairInfo(self, key, outfile):
		fp = open(outfile, 'a')    
		fp.write('\n\n *** cluster pair key: ' + str(key[0]) + ',' + str(key[1]))
		fp.write('\n relations [freq]: ')
		fp.write('\n [' + ' R1 ' + '/' + str(self.freq_count[RELATION_R1]) + \
			' R2 ' + '/' + str(self.freq_count[RELATION_R2]) + ' R3 ' + '/' + str(self.freq_count[RELATION_R3]) + \
				' R4 ' + '/' + str(self.freq_count[RELATION_R4]) + ']')
		fp.write('\n support scores: ')
		fp.write('\n [' + ' R1 ' + '/' + str(self.support_score_measure[RELATION_R1]) \
			+ ' R2 ' + '/' + str(self.support_score_measure[RELATION_R2]) \
			+ ' R3 ' + '/' + str(self.support_score_measure[RELATION_R3]) \
				+ ' R4 ' + '/' + str(self.support_score_measure[RELATION_R4]) + ']')
		fp.write('\n Allowed relation list: ' + str(self.possible_reln_list))
		fp.close()
		
	