#!/usr/bin/env python

import Header
from Header import *

##-------------------------------------------
#"""
#this function prints the tree in Newick format 
#latest version - 16th May 2016
#@parameters:
	#root_clust_node_idx: The cluster from which the node is traversed
	#queue_cluster: an array of node indices, which stores the current list of queued (but not explored) nodes
	
#In subsequent comments, we denote by cx, the current cluster (root_clust_node_idx)
#"""
#def PrintNewick(root_clust_node_idx, queue_cluster):
	#if 1:
		#print 'in function printnewick:   root_clust_node_idx: ', root_clust_node_idx
	#if 0:
		#print 'taxa set: ', Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()  
		#print 'out clust list: ', Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1)

	#Tree_Str_List = ''
	#"""
	#process the node provided it has not been explored yet
	#"""
	#if (Cluster_Info_Dict[root_clust_node_idx]._GetExploredStatus() == 0):  
		#"""
		#set the explored status of the current node to true
		#"""
		#Cluster_Info_Dict[root_clust_node_idx]._SetExploredStatus()
		#"""
		#this is the list of taxa of this cluster
		#"""
		#spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
		
		#"""
		#this is the parent list of the root cluster
		#"""
		#root_clust_R2_list = Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R2)
		#"""
		#this is the R1 list of the root cluster
		#"""
		#root_clust_R1_list = Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1)
		
		#"""
		#Here we form three different lists:
		#1) outnodes: contains clusters cy such that Cx->Cy holds
		#2) R4_R1_outnodes: contains clusters cy such that Cx--->Cy holds
		#3) R4_R1R2_outnodes: contains clusters cy such that both Cx--->Cy and Cx<--->Cy holds
		#Note that individual Cy clusters should not be explored
		#"""
		#outnodes = []
		#R4_R1_outnodes = []
		#R4_R1R2_outnodes = []
		
		#"""
		#obtain the list of clusters cy such that Cx->Cy holds
		#condition: Cy should not be already explored
		#"""
		#for l in root_clust_R1_list:
			#if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				#if l not in queue_cluster:
					#outnodes.append(l)
					#queue_cluster.append(l)
		
		#if 1:
			#print 'list of outnodes of the cluster : ',root_clust_node_idx, '  is:  ', str(outnodes)
		
		#"""
		#prepare following two lists:
		#1) R4_R1_outnodes: it contains the list of clusters Cy such that Cx--->Cy holds 
				#conditions: 1) Cy should not be already explored
										#2) Cy should have zero indegree
										#3) The final possible R2 list entry of Cy should contain Cx
		#2) R4_R1R2_outnodes: it contains the list of clusters Cy such that both Cx--->Cy and Cx<---Cy hold
				#conditions: 1) Cy should not be already explored
										#2A) Either Cy and Cx should have the same parent
										#2B) Or, Cy should have zero indegree
		#"""
		#for l in Cluster_Info_Dict[root_clust_node_idx]._GetPossibleR1List():
			#"""
			#check if the cluster l also has root_clust_node_idx as its final possible R2 cluster
			#"""
			#if (root_clust_node_idx in Cluster_Info_Dict[l]._GetFinalPossibleR2List()):
				#if 1:
					#print ' Possible R1 list entry: ', l, '  which contains final possible R2 list entry as well'
				#if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
					#if 1:
						#print ' The entry : ', l, '  is still not explored'
					#"""
					#explore the R2 list of this cluster
					#"""
					#l_R2_list = Cluster_Info_Dict[l]._GetClustRelnList(RELATION_R2)
					#if l in Cluster_Info_Dict[root_clust_node_idx]._GetPossibleR2List():	# Both Cx-->Cy and Cx<---Cy hold
						#if 1:
							#print ' The entry : ', l, '  also belongs to the possible R2 list of the cluster: ', root_clust_node_idx
						#flag = False
						#if (len(l_R2_list) == 0):	#l has zero indegree
							#flag = True
						#elif (len(root_clust_R2_list) > 0):	# the current cluster has > 0 indegree
							#if (root_clust_R2_list[0] == l_R2_list[0]):	# parents match for both clusters
								#flag = True
						#if (flag == True):
							#if l not in queue_cluster:
								#queue_cluster.append(l)
								#R4_R1R2_outnodes.append(l)
					#else:	# Only Cx-->Cy holds
						#if 1:
							#print ' The entry : ', l, '  only belongs to the possible R1 list of the cluster: ', root_clust_node_idx
						#if (len(l_R2_list) == 0):	#l has zero indegree
							#if 1:
								#print ' The entry : ', l, '  has zero indegree'
							#if l not in queue_cluster:
								#queue_cluster.append(l)
								#R4_R1_outnodes.append(l)
				#else:
					#if 1:
						#print ' The entry : ', l, '  is already explored'
					
				
		#if 1:
			#print 'list of R4_R1_outnodes: ', str(R4_R1_outnodes)
			#print 'list of R4_R1R2_outnodes: ', str(R4_R1R2_outnodes)
		
		##-------------------------------------------------------------------------
		#"""
		#if R4_R1R2_outnodes list is non empty, first start with an opening bracket 
		#since here two or more taxa clusters will be embedded
		#"""
		#if (len(R4_R1R2_outnodes) > 0):	# and (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + '('
		
		##-------------------------------------------------------------------------
		#"""
		#if R4_R1_outnodes list is non empty, first start with an opening bracket 
		#since here two or more taxa clusters will be embedded
		#"""
		#if (len(R4_R1_outnodes) > 0):	# and (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + '('
		
		##-------------------------------------------------------------------------
		#"""
		#if the cluster has one or more out edges
		#then enclose the string with opening and closing brackets
		#"""
		#if (len(outnodes) > 0):	# and (len(spec_list) == 1):
			#Tree_Str_List = Tree_Str_List + '('
		
		##-------------------------------------------------------------------------
		#""" 
		#at first, print the contents of this taxa cluster
		#if the cluster has more than one taxon, then use ( and ) to enclose the taxa list
		#"""
		#if (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + '('
		#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
		#if (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + ')'
		
		##-------------------------------------------------------------------------
		#"""
		#here we check if the cluster has one or more out edges
		#then recursively traverse all the out edge clusters
		#"""
		#if (len(outnodes) > 0):
			#"""
			#first add one comma
			#"""
			#Tree_Str_List = Tree_Str_List + ','
			#"""
			#then add one opening bracket, within which, all the out edge cluster contents will reside
			#"""
			#if (len(outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			
			#for i in range(len(outnodes)):
				#if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):  
					#if 1:
						#print 'Printing the cluster within outnodes: ', str(outnodes[i])
					#Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i], queue_cluster)
					#if (i < (len(outnodes) - 1)):
						#"""
						#we check whether any subsequent node belonging to the outnodes list
						#is left for traverse
						#"""
						#j = i + 1
						#while (j < len(outnodes)):
							#if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1
						#"""
						#in this case, we append one comma
						#"""
						#if (j < len(outnodes)):
							#Tree_Str_List = Tree_Str_List + ','      
			
			#"""
			#at last, append one closing bracket, signifying the end of out edge cluster contents
			#"""
			#if (len(outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + ')'

		##-------------------------------------------------------------------------
		#"""
		#if the cluster has one or more out edges
		#then enclose the string with opening and closing brackets
		#"""
		#if (len(outnodes) > 0):
			#Tree_Str_List = Tree_Str_List + ')'
		##-------------------------------------------------------------------------
		#"""
		#if R4_R1_outnodes list is non empty, print the contents of clusters 
		#which are connected by out edges to the current cluster
		#"""
		#if (len(R4_R1_outnodes) > 0):
			
			#Tree_Str_List = Tree_Str_List + ','
			#"""
			#opening bracket before printing the contents of R4 out edge clusters of the current cluster
			#"""
			#if (len(R4_R1_outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			
			#"""
			#navigate through individual out edge clusters connected to this cluster
			#and print their contents via recursive calling of this function
			#"""
			#for i in range(len(R4_R1_outnodes)):
				#if (Cluster_Info_Dict[R4_R1_outnodes[i]]._GetExploredStatus() == 0):  
					#if 1:
						#print 'Printing the cluster within R4_R1_outnodes: ', str(R4_R1_outnodes[i])
					#Tree_Str_List = Tree_Str_List + PrintNewick(R4_R1_outnodes[i], queue_cluster)
					#if (i < (len(R4_R1_outnodes) - 1)):
						#"""
						#we check whether any subsequent node belonging to the R4_outnodes list
						#is left for traverse
						#"""
						#j = i + 1
						#while (j < len(R4_R1_outnodes)):
							#if (Cluster_Info_Dict[R4_R1_outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1
						#"""
						#in this case, we append one comma
						#"""
						#if (j < len(R4_R1_outnodes)):
							#Tree_Str_List = Tree_Str_List + ','
		
			#"""
			#closing bracket after printing the contents of R4 out edge clusters of the current cluster
			#"""
			#if (len(R4_R1_outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + ')'
		
		##-------------------------------------------------------------------------
		#"""
		#if R4_R1_outnodes list is non empty, end with a closing bracket 
		#"""
		#if (len(R4_R1_outnodes) > 0):	# and (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + ')'
		##-------------------------------------------------------------------------
		#"""
		#if R4_R1R2_outnodes list is non empty, print the contents of clusters 
		#which are connected by out edges to the current cluster
		#Note: all these out edge connections were originally initiated by R4 consensus relation
		#"""
		#if (len(R4_R1R2_outnodes) > 0):
			#Tree_Str_List = Tree_Str_List + ','
			#"""
			#opening bracket before printing the contents of R4 out edge clusters of the current cluster
			#"""
			#if (len(R4_R1R2_outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			
			#"""
			#navigate through individual out edge clusters connected to this cluster
			#and print their contents via recursive calling of this function
			#"""
			#for i in range(len(R4_R1R2_outnodes)):
				#if (Cluster_Info_Dict[R4_R1R2_outnodes[i]]._GetExploredStatus() == 0):  
					#if 1:
						#print 'Printing the cluster within R4_R1R2_outnodes: ', str(R4_R1R2_outnodes[i])
					#Tree_Str_List = Tree_Str_List + PrintNewick(R4_R1R2_outnodes[i], queue_cluster)
					#if (i < (len(R4_R1R2_outnodes) - 1)):
						#"""
						#we check whether any subsequent node belonging to the R4_outnodes list
						#is left for traverse
						#"""
						#j = i + 1
						#while (j < len(R4_R1R2_outnodes)):
							#if (Cluster_Info_Dict[R4_R1R2_outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1
						#"""
						#in this case, we append one comma
						#"""
						#if (j < len(R4_R1R2_outnodes)):
							#Tree_Str_List = Tree_Str_List + ','
		
			#"""
			#closing bracket after printing the contents of R4 out edge clusters of the current cluster
			#"""
			#if (len(R4_R1R2_outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + ')'
		
		##-------------------------------------------------------------------------
		#"""
		#if R4_R1R2_outnodes list is non empty, end with a closing bracket 
		#"""
		#if (len(R4_R1R2_outnodes) > 0):	# and (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + ')'
		##--------------------------------------------------------------------------------------------------

	#if 1:
		#print 'Out of the function printnewick'
	
	#return Tree_Str_List    

##-------------------------------------------
# comment - sourya
##-------------------------------------------
#"""
#this function prints the tree in Newick format 
#adapted from COSPEDSPEC
#latest version - 2nd May 2016
#"""
#def PrintNewick(root_clust_node_idx, queue_cluster):
	#if 1:
		#print 'in function printnewick:   root_clust_node_idx: ', root_clust_node_idx
	#if 0:
		#print 'taxa set: ', Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()  
		#print 'out clust list: ', Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1)

	#Tree_Str_List = ''
	#"""
	#process the node provided it has not been explored yet
	#"""
	#if (Cluster_Info_Dict[root_clust_node_idx]._GetExploredStatus() == 0):  
		#"""
		#set the explored status of the current node to true
		#"""
		#Cluster_Info_Dict[root_clust_node_idx]._SetExploredStatus()
		#"""
		#this is the list of taxa of this cluster
		#"""
		#spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
		
		#"""
		#this is the parent list of the root cluster
		#"""
		#root_clust_R2_list = Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R2)
		#"""
		#this is the R1 list of the root cluster
		#"""
		#root_clust_R1_list = Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1)
		
		#"""
		#Here we form three different lists:
		#1) outnodes: contains clusters cy such that Cx->Cy holds
		#2) R4_R1_outnodes: contains clusters cy such that Cx--->Cy holds
		#3) R4_R1R2_outnodes: contains clusters cy such that both Cx--->Cy and Cx<--->Cy holds
		#Note that individual Cy clusters should not be explored
		#"""
		#outnodes = []
		#R4_R1_outnodes = []
		#R4_R1R2_outnodes = []
		
		#"""
		#obtain the list of clusters cy such that Cx->Cy holds
		#condition: Cy should not be already explored
		#"""
		#for l in root_clust_R1_list:
			#if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				#if l not in queue_cluster:
					#outnodes.append(l)
					#queue_cluster.append(l)
		
		#if 1:
			#print 'list of outnodes of the cluster : ',root_clust_node_idx, '  is:  ', str(outnodes)
		
		#"""
		#prepare following two lists:
		#1) R4_R1_outnodes: it contains the list of clusters Cy such that Cx--->Cy holds 
				#conditions: 1) Cy should not be already explored
										#2) Cy should have zero indegree
										#3) The final possible R2 list entry of Cy should contain Cx
		#2) R4_R1R2_outnodes: it contains the list of clusters Cy such that both Cx--->Cy and Cx<---Cy hold
				#conditions: 1) Cy should not be already explored
										#2A) Either Cy and Cx should have the same parent
										#2B) Or, Cy should have zero indegree
		#"""
		#for l in Cluster_Info_Dict[root_clust_node_idx]._GetPossibleR1List():
			#"""
			#check if the cluster l also has root_clust_node_idx as its final possible R2 cluster
			#"""
			#if (root_clust_node_idx in Cluster_Info_Dict[l]._GetFinalPossibleR2List()):
				#if 1:
					#print ' Possible R1 list entry: ', l, '  which contains final possible R2 list entry as well'
				#if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
					#if 1:
						#print ' The entry : ', l, '  is still not explored'
					#"""
					#explore the R2 list of this cluster
					#"""
					#l_R2_list = Cluster_Info_Dict[l]._GetClustRelnList(RELATION_R2)
					#if l in Cluster_Info_Dict[root_clust_node_idx]._GetPossibleR2List():	# Both Cx-->Cy and Cx<---Cy hold
						#if 1:
							#print ' The entry : ', l, '  also belongs to the possible R2 list of the cluster: ', root_clust_node_idx
						#flag = False
						#if (len(l_R2_list) == 0):	#l has zero indegree
							#flag = True
						#elif (len(root_clust_R2_list) > 0):	# the current cluster has > 0 indegree
							#if (root_clust_R2_list[0] == l_R2_list[0]):	# parents match for both clusters
								#flag = True
						#if (flag == True):
							#if l not in queue_cluster:
								#queue_cluster.append(l)
								#R4_R1R2_outnodes.append(l)
					#else:	# Only Cx-->Cy holds
						#if 1:
							#print ' The entry : ', l, '  only belongs to the possible R1 list of the cluster: ', root_clust_node_idx
						#if (len(l_R2_list) == 0):	#l has zero indegree
							#if 1:
								#print ' The entry : ', l, '  has zero indegree'
							#if l not in queue_cluster:
								#queue_cluster.append(l)
								#R4_R1_outnodes.append(l)
				#else:
					#if 1:
						#print ' The entry : ', l, '  is already explored'
					
				
		#if 1:
			#print 'list of R4_R1_outnodes: ', str(R4_R1_outnodes)
			#print 'list of R4_R1R2_outnodes: ', str(R4_R1R2_outnodes)
		
		##"""
		##if R4_R1_outnodes list is non empty, first start with an opening bracket 
		##since here two or more taxa clusters will be embedded
		##"""
		##if (len(R4_R1_outnodes) > 0) and (len(spec_list) > 1):
			##Tree_Str_List = Tree_Str_List + '('
		
		##"""
		##if R4_R1R2_outnodes list is non empty, first start with an opening bracket 
		##since here two or more taxa clusters will be embedded
		##"""
		##if (len(R4_R1R2_outnodes) > 0) and (len(spec_list) > 1):
			##Tree_Str_List = Tree_Str_List + '('

		#R4_outnodes_flag = False
		#if ((len(R4_R1_outnodes) > 0) or (len(R4_R1R2_outnodes) > 0)) and (len(spec_list) > 1):
			#R4_outnodes_flag = True

		#if (R4_outnodes_flag == True):
			#Tree_Str_List = Tree_Str_List + '('

		#"""
		#if the cluster has one or more out edges
		#then enclose the string with opening and closing brackets
		#"""
		#if (len(outnodes) > 0):	# and (R4_outnodes_flag == False):	# and (len(spec_list) == 1):
			#Tree_Str_List = Tree_Str_List + '('
		
		#""" 
		#at first, print the contents of this taxa cluster
		#if the cluster has more than one taxon, then use ( and ) to enclose the taxa list
		#"""
		#if (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + '('
		#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
		#if (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + ')'
		
		#"""
		#here we check if the cluster has one or more out edges
		#then recursively traverse all the out edge clusters
		#"""
		#if (len(outnodes) > 0):
			#"""
			#first add one comma
			#"""
			#Tree_Str_List = Tree_Str_List + ','
			
			#"""
			#then add one opening bracket, within which, all the out edge cluster contents will reside
			#"""
			#if (len(outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			
			#for i in range(len(outnodes)):
				#if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):  
					#if 1:
						#print 'Printing the cluster within outnodes: ', str(outnodes[i])
					#Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i], queue_cluster)
					#if (i < (len(outnodes) - 1)):
						#"""
						#we check whether any subsequent node belonging to the outnodes list
						#is left for traverse
						#"""
						#j = i + 1
						#while (j < len(outnodes)):
							#if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1
						#"""
						#in this case, we append one comma
						#"""
						#if (j < len(outnodes)):
							#Tree_Str_List = Tree_Str_List + ','      
			
			#"""
			#at last, append one closing bracket, signifying the end of out edge cluster contents
			#"""
			#if (len(outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + ')'

		#"""
		#if the cluster has one or more out edges
		#then enclose the string with opening and closing brackets
		#"""
		#if (len(outnodes) > 0):	# and (R4_outnodes_flag == False):	# and (len(spec_list) == 1):
			#Tree_Str_List = Tree_Str_List + ')'
	
		##-------------------------------------------------------------------------
		
		## add - sourya
		#R4_club_flag = False
		
		#"""
		#if R4_R1R2_outnodes list is non empty, print the contents of clusters 
		#which are connected by out edges to the current cluster
		#Note: all these out edge connections were originally initiated by R4 consensus relation
		#"""
		#if (len(R4_R1R2_outnodes) > 0):
			
			#Tree_Str_List = Tree_Str_List + ','
			
			## add - sourya
			#if (R4_club_flag == False) and ((len(R4_R1R2_outnodes) + len(R4_R1_outnodes)) > 1):
				#Tree_Str_List = Tree_Str_List + '('
				#R4_club_flag = True
			
			#"""
			#opening bracket before printing the contents of R4 out edge clusters of the current cluster
			#"""
			#if (len(R4_R1R2_outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			
			#"""
			#navigate through individual out edge clusters connected to this cluster
			#and print their contents via recursive calling of this function
			#"""
			#for i in range(len(R4_R1R2_outnodes)):
				#if (Cluster_Info_Dict[R4_R1R2_outnodes[i]]._GetExploredStatus() == 0):  
					#if 1:
						#print 'Printing the cluster within R4_R1R2_outnodes: ', str(R4_R1R2_outnodes[i])
					#Tree_Str_List = Tree_Str_List + PrintNewick(R4_R1R2_outnodes[i], queue_cluster)
					#if (i < (len(R4_R1R2_outnodes) - 1)):
						#"""
						#we check whether any subsequent node belonging to the R4_outnodes list
						#is left for traverse
						#"""
						#j = i + 1
						#while (j < len(R4_R1R2_outnodes)):
							#if (Cluster_Info_Dict[R4_R1R2_outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1
						#"""
						#in this case, we append one comma
						#"""
						#if (j < len(R4_R1R2_outnodes)):
							#Tree_Str_List = Tree_Str_List + ','
		
			#"""
			#closing bracket after printing the contents of R4 out edge clusters of the current cluster
			#"""
			#if (len(R4_R1R2_outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + ')'
		
		##"""
		##if R4_R1R2_outnodes list is non empty, end with a closing bracket 
		##"""
		##if (len(R4_R1R2_outnodes) > 0) and (len(spec_list) > 1):
			##Tree_Str_List = Tree_Str_List + ')'
		
		##--------------------------------------------------------------------------------------------------
		
		#"""
		#if R4_R1_outnodes list is non empty, print the contents of clusters 
		#which are connected by out edges to the current cluster
		#Note: all these out edge connections were originally initiated by R4 consensus relation
		#"""
		#if (len(R4_R1_outnodes) > 0):
			
			#Tree_Str_List = Tree_Str_List + ','
			
			## add - sourya
			#if (R4_club_flag == False) and ((len(R4_R1R2_outnodes) + len(R4_R1_outnodes)) > 1):
				#Tree_Str_List = Tree_Str_List + '('
				#R4_club_flag = True
			
			#"""
			#opening bracket before printing the contents of R4 out edge clusters of the current cluster
			#"""
			#if (len(R4_R1_outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			
			#"""
			#navigate through individual out edge clusters connected to this cluster
			#and print their contents via recursive calling of this function
			#"""
			#for i in range(len(R4_R1_outnodes)):
				#if (Cluster_Info_Dict[R4_R1_outnodes[i]]._GetExploredStatus() == 0):  
					#if 1:
						#print 'Printing the cluster within R4_R1_outnodes: ', str(R4_R1_outnodes[i])
					#Tree_Str_List = Tree_Str_List + PrintNewick(R4_R1_outnodes[i], queue_cluster)
					#if (i < (len(R4_R1_outnodes) - 1)):
						#"""
						#we check whether any subsequent node belonging to the R4_outnodes list
						#is left for traverse
						#"""
						#j = i + 1
						#while (j < len(R4_R1_outnodes)):
							#if (Cluster_Info_Dict[R4_R1_outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1
						#"""
						#in this case, we append one comma
						#"""
						#if (j < len(R4_R1_outnodes)):
							#Tree_Str_List = Tree_Str_List + ','
		
			#"""
			#closing bracket after printing the contents of R4 out edge clusters of the current cluster
			#"""
			#if (len(R4_R1_outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + ')'
		
		##"""
		##if R4_R1_outnodes list is non empty, end with a closing bracket 
		##"""
		##if (len(R4_R1_outnodes) > 0) and (len(spec_list) > 1):
			##Tree_Str_List = Tree_Str_List + ')'

		## add - sourya
		#if ((len(R4_R1R2_outnodes) + len(R4_R1_outnodes)) > 1):
			#Tree_Str_List = Tree_Str_List + ')'
		##--------------------------------------------------------------------------------------------------

		#if (R4_outnodes_flag == True):
			#Tree_Str_List = Tree_Str_List + ')'

	#if 1:
		#print 'Out of the function printnewick'
	
	#return Tree_Str_List    

##-------------------------------------------------
##************** comment - sourya
##-------------------------------------------
#"""
#this function prints the tree in Newick format - old function
#adapted from COSPEDSPEC
#"""
#def PrintNewick(root_clust_node_idx, queue_cluster):
	#if 1:
		#print 'in function printnewick:   root_clust_node_idx: ', root_clust_node_idx
	#if 0:
		#print 'taxa set: ', Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()  
		#print 'out clust list: ', Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1)

	#Tree_Str_List = ''
	#"""
	#process the node provided it has not been explored yet
	#"""
	#if (Cluster_Info_Dict[root_clust_node_idx]._GetExploredStatus() == 0):  
		#"""
		#set the explored status of the current node to true
		#"""
		#Cluster_Info_Dict[root_clust_node_idx]._SetExploredStatus()
		#"""
		#this is the list of taxa of this cluster
		#"""
		#spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
		#"""
		#Here we form three different lists:
		#1) outnodes: contains clusters cy such that Cx->Cy holds
		#2) R4_R1_outnodes: contains clusters cy such that Cx--->Cy holds
		#3) R4_R1R2_outnodes: contains clusters cy such that Cx<--->Cy holds
		#Note that individual Cy clusters should not be explored
		#"""
		#outnodes = []
		#R4_R1_outnodes = []
		#R4_R1R2_outnodes = []
		##R4_outnodes = []
		
		#"""
		#get the out edge list of the current cluster which are not explored yet 
		#"""
		#for l in Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1):
			#if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				#if l not in queue_cluster:
					#outnodes.append(l)
					#queue_cluster.append(l)
		
		#for l in Cluster_Info_Dict[root_clust_node_idx]._GetPossibleR1List():
			#if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				#if (len(Cluster_Info_Dict[l]._GetClustRelnList(RELATION_R2)) == 0):
					#if l not in queue_cluster:
						#queue_cluster.append(l)
						#if l in Cluster_Info_Dict[root_clust_node_idx]._GetPossibleR2List():
							#R4_R1R2_outnodes.append(l)
						#else:
							#R4_R1_outnodes.append(l)
						##R4_outnodes.append(l)
		
		#if 1:
			#print 'list of outnodes: ', str(outnodes)
			#print 'list of R4_R1_outnodes: ', str(R4_R1_outnodes)
			#print 'list of R4_R1R2_outnodes: ', str(R4_R1R2_outnodes)
			##print 'list of R4_outnodes: ', str(R4_outnodes)
		
		#"""
		#if R4_R1_outnodes list is non empty, first start with an opening bracket 
		#since here two or more taxa clusters will be embedded
		#"""
		#if (len(R4_R1_outnodes) > 0):	# and (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + '('
		
		#"""
		#if R4_R1R2_outnodes list is non empty, first start with an opening bracket 
		#since here two or more taxa clusters will be embedded
		#"""
		#if (len(R4_R1R2_outnodes) > 0):	# and (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + '('

		##"""
		##if R4_outnodes list is non empty, first start with an opening bracket 
		##since here two or more taxa clusters will be embedded
		##"""
		##if (len(R4_outnodes) > 0):	# and (len(spec_list) > 1):
			##Tree_Str_List = Tree_Str_List + '('

		#"""
		#if the cluster has one or more out edges
		#then enclose the string with opening and closing brackets
		#"""
		#if (len(outnodes) > 0):	# and (len(spec_list) == 1):
			#Tree_Str_List = Tree_Str_List + '('
		
		#""" 
		#at first, print the contents of this taxa cluster
		#if the cluster has more than one taxon, then use ( and ) to enclose the taxa list
		#"""
		#if (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + '('
		#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
		#if (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + ')'
		
		#"""
		#here we check if the cluster has one or more out edges
		#then recursively traverse all the out edge clusters
		#"""
		#if (len(outnodes) > 0):
			#"""
			#first add one comma
			#"""
			#Tree_Str_List = Tree_Str_List + ','
			
			#"""
			#then add one opening bracket, within which, all the out edge cluster contents will reside
			#"""
			#if (len(outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			
			#for i in range(len(outnodes)):
				#if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):  
					#if 1:
						#print 'Printing the cluster within outnodes: ', str(outnodes[i])
					#Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i], queue_cluster)
					#if (i < (len(outnodes) - 1)):
						#"""
						#we check whether any subsequent node belonging to the outnodes list
						#is left for traverse
						#"""
						#j = i + 1
						#while (j < len(outnodes)):
							#if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1
						#"""
						#in this case, we append one comma
						#"""
						#if (j < len(outnodes)):
							#Tree_Str_List = Tree_Str_List + ','      
			
			#"""
			#at last, append one closing bracket, signifying the end of out edge cluster contents
			#"""
			#if (len(outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + ')'

		#"""
		#if the cluster has one or more out edges
		#then enclose the string with opening and closing brackets
		#"""
		#if (len(outnodes) > 0):	# and (len(spec_list) == 1):
			#Tree_Str_List = Tree_Str_List + ')'
	
		###-------------------------------------------------------------------------
		##"""
		##if R4_outnodes list is non empty, print the contents of clusters 
		##which are connected by out edges to the current cluster
		##Note: all these out edge connections were originally initiated by R4 consensus relation
		##"""
		##if (len(R4_outnodes) > 0):
			
			##Tree_Str_List = Tree_Str_List + ','
			
			##"""
			##opening bracket before printing the contents of R4 out edge clusters of the current cluster
			##"""
			##if (len(R4_outnodes) > 1):
				##Tree_Str_List = Tree_Str_List + '('
			
			##"""
			##navigate through individual out edge clusters connected to this cluster
			##and print their contents via recursive calling of this function
			##"""
			##for i in range(len(R4_outnodes)):
				##if (Cluster_Info_Dict[R4_outnodes[i]]._GetExploredStatus() == 0):  
					##if 0:
						##print 'Printing the cluster within R4_outnodes: ', str(R4_outnodes[i])
					##Tree_Str_List = Tree_Str_List + PrintNewick(R4_outnodes[i], queue_cluster)
					##if (i < (len(R4_outnodes) - 1)):
						##"""
						##we check whether any subsequent node belonging to the R4_outnodes list
						##is left for traverse
						##"""
						##j = i + 1
						##while (j < len(R4_outnodes)):
							##if (Cluster_Info_Dict[R4_outnodes[j]]._GetExploredStatus() == 0):  
								##break
							##j = j + 1
						##"""
						##in this case, we append one comma
						##"""
						##if (j < len(R4_outnodes)):
							##Tree_Str_List = Tree_Str_List + ','
		
			##"""
			##closing bracket after printing the contents of R4 out edge clusters of the current cluster
			##"""
			##if (len(R4_outnodes) > 1):
				##Tree_Str_List = Tree_Str_List + ')'
		
		##"""
		##if R4_outnodes list is non empty, end with a closing bracket 
		##"""
		##if (len(R4_outnodes) > 0):	# and (len(spec_list) > 1):
			##Tree_Str_List = Tree_Str_List + ')'
	
		##-------------------------------------------------------------------------
		#"""
		#if R4_R1R2_outnodes list is non empty, print the contents of clusters 
		#which are connected by out edges to the current cluster
		#Note: all these out edge connections were originally initiated by R4 consensus relation
		#"""
		#if (len(R4_R1R2_outnodes) > 0):
			
			#Tree_Str_List = Tree_Str_List + ','
			
			#"""
			#opening bracket before printing the contents of R4 out edge clusters of the current cluster
			#"""
			#if (len(R4_R1R2_outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			
			#"""
			#navigate through individual out edge clusters connected to this cluster
			#and print their contents via recursive calling of this function
			#"""
			#for i in range(len(R4_R1R2_outnodes)):
				#if (Cluster_Info_Dict[R4_R1R2_outnodes[i]]._GetExploredStatus() == 0):  
					#if 1:
						#print 'Printing the cluster within R4_R1R2_outnodes: ', str(R4_R1R2_outnodes[i])
					#Tree_Str_List = Tree_Str_List + PrintNewick(R4_R1R2_outnodes[i], queue_cluster)
					#if (i < (len(R4_R1R2_outnodes) - 1)):
						#"""
						#we check whether any subsequent node belonging to the R4_outnodes list
						#is left for traverse
						#"""
						#j = i + 1
						#while (j < len(R4_R1R2_outnodes)):
							#if (Cluster_Info_Dict[R4_R1R2_outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1
						#"""
						#in this case, we append one comma
						#"""
						#if (j < len(R4_R1R2_outnodes)):
							#Tree_Str_List = Tree_Str_List + ','
		
			#"""
			#closing bracket after printing the contents of R4 out edge clusters of the current cluster
			#"""
			#if (len(R4_R1R2_outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + ')'
		
		#"""
		#if R4_R1R2_outnodes list is non empty, end with a closing bracket 
		#"""
		#if (len(R4_R1R2_outnodes) > 0):	# and (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + ')'
		
		##--------------------------------------------------------------------------------------------------
		
		#"""
		#if R4_R1_outnodes list is non empty, print the contents of clusters 
		#which are connected by out edges to the current cluster
		#Note: all these out edge connections were originally initiated by R4 consensus relation
		#"""
		#if (len(R4_R1_outnodes) > 0):
			
			#Tree_Str_List = Tree_Str_List + ','
			
			#"""
			#opening bracket before printing the contents of R4 out edge clusters of the current cluster
			#"""
			#if (len(R4_R1_outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			
			#"""
			#navigate through individual out edge clusters connected to this cluster
			#and print their contents via recursive calling of this function
			#"""
			#for i in range(len(R4_R1_outnodes)):
				#if (Cluster_Info_Dict[R4_R1_outnodes[i]]._GetExploredStatus() == 0):  
					#if 1:
						#print 'Printing the cluster within R4_R1_outnodes: ', str(R4_R1_outnodes[i])
					#Tree_Str_List = Tree_Str_List + PrintNewick(R4_R1_outnodes[i], queue_cluster)
					#if (i < (len(R4_R1_outnodes) - 1)):
						#"""
						#we check whether any subsequent node belonging to the R4_outnodes list
						#is left for traverse
						#"""
						#j = i + 1
						#while (j < len(R4_R1_outnodes)):
							#if (Cluster_Info_Dict[R4_R1_outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1
						#"""
						#in this case, we append one comma
						#"""
						#if (j < len(R4_R1_outnodes)):
							#Tree_Str_List = Tree_Str_List + ','
		
			#"""
			#closing bracket after printing the contents of R4 out edge clusters of the current cluster
			#"""
			#if (len(R4_R1_outnodes) > 1):
				#Tree_Str_List = Tree_Str_List + ')'
		
		#"""
		#if R4_R1_outnodes list is non empty, end with a closing bracket 
		#"""
		#if (len(R4_R1_outnodes) > 0):	# and (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + ')'

		##--------------------------------------------------------------------------------------------------

	#if 1:
		#print 'Out of the function printnewick'
	
	#return Tree_Str_List    

##-------------------------------------------
##************** end comment - sourya
##-----------------------------------------------------
#"""
#this is a modified function - sourya - April 7, 2016
#this function prints the tree in Newick format
#"""
#def PrintNewick(root_clust_node_idx):
	#if 1:
		#print 'in function printnewick:   root_clust_node_idx: ', str(root_clust_node_idx)
	#if 0:
		#print 'taxa set: ', Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()  
		#print 'out clust list: ', Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1)

	#"""
	#initialize the tree representative newick string
	#"""
	#Tree_Str_List = ''
	
	#"""
	#process the node provided it has not been explored yet
	#"""
	#if (Cluster_Info_Dict[root_clust_node_idx]._GetExploredStatus() == 0):  
		#"""
		#set the explored status of the current node to true
		#"""
		#Cluster_Info_Dict[root_clust_node_idx]._SetExploredStatus()
		
		#"""
		#here we explore the out edge list of the current cluster Cx to other clusters Cy
		#however, we distinguish between two types of out edges
		#1) where the edge Cx->Cy was formed due to the consensus R1 relation between Cx and Cy
		#these clusters Cy will go in a list of nodes "outnodes"
		#2) where the edge Cx->Cy was formed due to the consensus R4 relation between Cx and Cy
		#these clusters Cy were originally placed in the candidate out edge list of Cx
		#here these clusters are stored in a list called "R4_outnodes"
		
		#Note that individual Cy clusters should not be explored
		#"""
		#outnodes = []
		#R4_outnodes = []
		
		#for l in Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1):
			#if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				#outnodes.append(l)
				
		#for l in Cluster_Info_Dict[root_clust_node_idx]._GetPossibleR1List():
			#if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				#if (len(Cluster_Info_Dict[l]._GetClustRelnList(RELATION_R2)) == 0):
					#R4_outnodes.append(l)
		
		#print 'list of outnodes: ', str(outnodes)
		#print 'list of R4_outnodes: ', str(R4_outnodes)
		
		#"""
		#taxa list of the current cluster
		#"""
		#spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
		
		#"""
		#if R4_outnodes list is non empty, first start with an opening bracket 
		#since here two or more taxa clusters will be embedded
		#"""
		#if (len(R4_outnodes) > 0):
			#Tree_Str_List = Tree_Str_List + '('
		
		#"""
		#if the current cluster has no direct descendant, we just print its constituent species list
		#otherwise, we enclose the newick string with additional set of opening and closing bracket
		#"""
		#if (len(outnodes) == 0):
			#if (len(spec_list) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
			#if (len(spec_list) > 1):
				#Tree_Str_List = Tree_Str_List + ')'
		#else:
			#"""
			#this is the opening bracket of additional set 
			#since the cluster has one or more out clusters
			#"""
			#Tree_Str_List = Tree_Str_List + '('
			
			#"""
			#print the constituent species list
			#as there is one or more out edge clusters, we do not need to check about the cardinality of the species set 
			#or insert additional set of brackets
			#"""
			#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
			#Tree_Str_List = Tree_Str_List + ','  
			
			#"""
			#opening bracket before printing the contents of out edge clusters of the current cluster
			#"""
			#Tree_Str_List = Tree_Str_List + '('
			
			#"""
			#navigate through individual out edge clusters connected to this cluster
			#and print their contents via recursive calling of this function
			#"""
			#for i in range(len(outnodes)):
				#if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):
					#print 'Printing the cluster within outnodes: ', str(outnodes[i])
					#Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i])
					#if (i < (len(outnodes) - 1)):
						#"""
						#we check whether any subsequent node belonging to the outnodes list
						#is left for traverse
						#"""
						#j = i + 1
						#while (j < len(outnodes)):
							#if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1	      
						#"""
						#in this case, we append one comma
						#"""
						#if (j < len(outnodes)):
							#Tree_Str_List = Tree_Str_List + ','
			
			#"""
			#closing bracket after printing the contents of out edge clusters of the current cluster
			#"""
			#Tree_Str_List = Tree_Str_List + ')'
			
			#"""
			#this is the closing bracket of additional set 
			#since the cluster has one or more out clusters
			#"""
			#Tree_Str_List = Tree_Str_List + ')'
		
		#"""
		#if R4_outnodes list is non empty, print the contents of clusters 
		#which are connected by out edges to the current cluster
		#Note: all these out edge connections were originally initiated by R4 consensus relation
		#"""
		#if (len(R4_outnodes) > 0):
			
			#Tree_Str_List = Tree_Str_List + ','
			
			#"""
			#opening bracket before printing the contents of R4 out edge clusters of the current cluster
			#"""
			#Tree_Str_List = Tree_Str_List + '('
			
			#"""
			#navigate through individual out edge clusters connected to this cluster
			#and print their contents via recursive calling of this function
			#"""
			#for i in range(len(R4_outnodes)):
				#if (Cluster_Info_Dict[R4_outnodes[i]]._GetExploredStatus() == 0): 
					#print 'Printing the cluster within R4 outnodes: ', str(R4_outnodes[i])
					#Tree_Str_List = Tree_Str_List + PrintNewick(R4_outnodes[i])
					#if (i < (len(R4_outnodes) - 1)):
						#"""
						#we check whether any subsequent node belonging to the R4_outnodes list
						#is left for traverse
						#"""
						#j = i + 1
						#while (j < len(R4_outnodes)):
							#if (Cluster_Info_Dict[R4_outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1
						#"""
						#in this case, we append one comma
						#"""
						#if (j < len(R4_outnodes)):
							#Tree_Str_List = Tree_Str_List + ','
		
			#"""
			#closing bracket after printing the contents of R4 out edge clusters of the current cluster
			#"""
			#Tree_Str_List = Tree_Str_List + ')'
		
		#"""
		#if R4_outnodes list is non empty, end with a closing bracket 
		#"""
		#if (len(R4_outnodes) > 0):
			#Tree_Str_List = Tree_Str_List + ')'
	
	#if 1:
		#print 'Out of the function printnewick'
	
	#return Tree_Str_List    

#--------------------------------------------------------
#"""
#this function prints the tree in Newick format - old function
#this was used in latest code - 13.06.2016
#"""
#def PrintNewick(root_clust_node_idx):
	#if 0:
		#print 'in function printnewick:   root_clust_node_idx: ', root_clust_node_idx
		#print 'taxa set: ', Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()  
		#print 'out clust list: ', Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1)

	#Tree_Str_List = ''
	#"""
	#process the node provided it has not been explored yet
	#"""
	#if (Cluster_Info_Dict[root_clust_node_idx]._GetExploredStatus() == 0):  
		#"""
		#set the explored status of the current node to true
		#"""
		#Cluster_Info_Dict[root_clust_node_idx]._SetExploredStatus()
		#"""
		#get the out edge list of the current node which are not explored yet 
		#"""
		#outnodes = []
		#for l in Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1):
			#if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				#outnodes.append(l)
		## comment - sourya
		#if (len(outnodes) == 0):
		## add - sourya
		##if (len(outnodes) <= 1):
			#spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
			#if (len(spec_list) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
			#if (len(spec_list) > 1):
				#Tree_Str_List = Tree_Str_List + ')'
		#else:
			#Tree_Str_List = Tree_Str_List + '('
			#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList())
			#Tree_Str_List = Tree_Str_List + ','    
			#Tree_Str_List = Tree_Str_List + '('
			#for i in range(len(outnodes)):
				#if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):  
					#Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i])
					#if (i < (len(outnodes) - 1)):
						#"""
						#we check whether any subsequent node belonging to the outnodes list
						#is left for traverse
						#"""
						#j = i + 1
						#while (j < len(outnodes)):
							#if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1
						#"""
						#in this case, we append one comma
						#"""
						#if (j < len(outnodes)):
							#Tree_Str_List = Tree_Str_List + ','
			
			#Tree_Str_List = Tree_Str_List + ')'
			#Tree_Str_List = Tree_Str_List + ')'
		
	#return Tree_Str_List    

#-------------------------------------------
"""
this function prints the tree in Newick format - 
this was used in the code of COSPEDSPEC
"""
def PrintNewick(root_clust_node_idx):
	if 0:
		print 'in function printnewick:   root_clust_node_idx: ', root_clust_node_idx
		print 'taxa set: ', Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()  
		print 'out clust list: ', Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1)

	Tree_Str_List = ''
	"""
	process the node provided it has not been explored yet
	"""
	if (Cluster_Info_Dict[root_clust_node_idx]._GetExploredStatus() == 0):  
		"""
		set the explored status of the current node to true
		"""
		Cluster_Info_Dict[root_clust_node_idx]._SetExploredStatus()
		"""
		this is the list of taxa of this cluster
		"""
		spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
		"""
		get the out edge list of the current cluster which are not explored yet 
		"""
		outnodes = []
		for l in Cluster_Info_Dict[root_clust_node_idx]._GetClustRelnList(RELATION_R1):
			if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				outnodes.append(l)
		
		""" 
		at first, print the contents of this taxa cluster
		if the cluster has more than one taxon, then use ( and ) to enclose the taxa list
		"""
		if (len(outnodes) > 0):	# and (len(spec_list) == 1):
			Tree_Str_List = Tree_Str_List + '('
			
		if (len(spec_list) > 1):
			Tree_Str_List = Tree_Str_List + '('
		Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
		if (len(spec_list) > 1):
			Tree_Str_List = Tree_Str_List + ')'
		"""
		here we check if the cluster has one or more out edges
		then recursively traverse all the out edge clusters
		"""
		if (len(outnodes) > 0):
			"""
			first add one comma
			"""
			Tree_Str_List = Tree_Str_List + ','
			
			# then add one opening bracket, within which, all the out edge cluster contents will reside
			Tree_Str_List = Tree_Str_List + '('
			
			for i in range(len(outnodes)):
				if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):  
					Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i])
					if (i < (len(outnodes) - 1)):
						"""
						we check whether any subsequent node belonging to the outnodes list
						is left for traverse
						"""
						j = i + 1
						while (j < len(outnodes)):
							if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
								break
							j = j + 1
						"""
						in this case, we append one comma
						"""
						if (j < len(outnodes)):
							Tree_Str_List = Tree_Str_List + ','      
			
			"""
			at last, append one closing bracket, signifying the end of out edge cluster contents
			"""
			Tree_Str_List = Tree_Str_List + ')'

		if (len(outnodes) > 0):	# and (len(spec_list) == 1):
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
"""
def DefineLeafPairReln(xl_val, lca_level, lca_rank, node1, node2, reln_type, curr_tree_taxa, wt_taxa_subset):

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
		correct_order = True
	else:
		target_key = (node2_idx, node1_idx)
		target_reln_type = Complementary_Reln(reln_type)
		correct_order = False
		
	"""
	check if the key exists - else first create the key
	"""
	if target_key not in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict.setdefault(target_key, Reln_TaxaPair())

	"""
	THE VARIABLE "intersect_ratio" is used for the weighted frequency value
	sourya - change the variable value for explore
	"""
	if (wt_taxa_subset == True):
		und_tax_list = TaxaPair_Reln_Dict[target_key]._GetUnderlyingTaxonList()
		intersect_taxa_set = set(und_tax_list) & set(curr_tree_taxa)
		intersect_ratio = (len(intersect_taxa_set) * 1.0) / len(und_tax_list)
	else:
		intersect_ratio = 1	#(len(curr_tree_taxa) * 1.0) / len(COMPLETE_INPUT_TAXA_LIST)	#1	#lca_rank

	TaxaPair_Reln_Dict[target_key]._AddSupportingTree()
	TaxaPair_Reln_Dict[target_key]._AddXLVal(xl_val / intersect_ratio)
	TaxaPair_Reln_Dict[target_key]._AddEdgeCount(target_reln_type, intersect_ratio)
	TaxaPair_Reln_Dict[target_key]._AddLevel(sum_of_branch_count)
	TaxaPair_Reln_Dict[target_key]._AddLCARank(lca_rank)
	
	#-----------------------
	if (target_reln_type == RELATION_R4):
		if ((node1_level - lca_level) == 2) and ((node2_level - lca_level) > 2):
			if (correct_order == False):
				TaxaPair_Reln_Dict[target_key]._AddFreqPseudoR1(1, intersect_ratio)	#1)
			else:
				TaxaPair_Reln_Dict[target_key]._AddFreqPseudoR1(0, intersect_ratio)	#1)
		elif ((node1_level - lca_level) > 2) and ((node2_level - lca_level) == 2):
			if (correct_order == False):
				TaxaPair_Reln_Dict[target_key]._AddFreqPseudoR1(0, intersect_ratio)	#1)
			else:
				TaxaPair_Reln_Dict[target_key]._AddFreqPseudoR1(1, intersect_ratio)	#1)
		elif ((node1_level - lca_level) == 2) and ((node2_level - lca_level) == 2):
			TaxaPair_Reln_Dict[target_key]._AddFreqPseudoR1(2, intersect_ratio)	#1)
	#-----------------------
	
	return

#--------------------------------------------------------
"""
this function derives couplet relations belonging to one tree
that is provided as an input argument to this function
@parameters:  
	WEIGHT_TAXA_SUBSET: If True, the relation takes care of the 
											set of taxa underlying the LCA node for this couplet
"""
def DeriveCoupletRelations(Curr_tree, WEIGHT_TAXA_SUBSET):
  
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
		
		# modified - sourya
		curr_node_rank = no_of_taxa - curr_node_level
		
		# we have used normalized LCA rank
		#curr_node_rank = ((no_of_taxa - curr_node_level) * 1.0) / no_of_taxa
		
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
						RELATION_R3, curr_tree_taxa, WEIGHT_TAXA_SUBSET)
		
		"""
		one leaf node (direct descendant) and another leaf node (under one internal node)
		will be related by ancestor / descendant relations
		"""
		if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
			for p in curr_node_child_leaf_nodes:
				for q in curr_node_child_internal_nodes:
					for r in q.leaf_nodes():
						DefineLeafPairReln(xl_val, curr_node_level, curr_node_rank, p, r, RELATION_R1, \
							curr_tree_taxa, WEIGHT_TAXA_SUBSET)
		
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
								curr_tree_taxa, WEIGHT_TAXA_SUBSET)

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
		modified - sourya
		we list the indices of the taxa, not their complete names
		to save the memory
		"""
		#taxa_under_curr_node = GetTaxaUnderInternalNode(curr_node)
		taxa_under_curr_node = GetTaxa_IDX_UnderInternalNode(curr_node)

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

##--------------------------------------------------------
#""" 
#this function prints the elements of the queue (which stores the couplet scores 
#for individual relations 
#"""
#def PrintQueueInfo(inp_queue, Output_Text_File):
	#fp = open(Output_Text_File, 'a')
	#for elem in inp_queue:
		#fp.write('\n' + str(elem))
	#fp.close()

#-----------------------------------------------------
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

#-----------------------------------------------------
"""
this function returns the list of taxa (in terms of their indices in the COMPLETE_INPUT_TAXA_LIST) 
underlying the given internal node
@param: curr_node: Input node under which the taxa set will be explored
"""
def GetTaxa_IDX_UnderInternalNode(curr_node):
	taxa_list_from_curr_internal_node = []
	for n in curr_node.leaf_nodes():
		label = n.taxon.label
		idx = COMPLETE_INPUT_TAXA_LIST.index(label)
		taxa_list_from_curr_internal_node.append(idx)
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
this function computes average distance measure between a pair of taxa clusters
used for binary refinement of the supertree
@param: 
1) taxa_clust1: first taxa list
2) taxa_clust2: second taxa list
3) DIST_MAT_TYPE: Type of distance employed
		1) mean of XL measure
		2) mean (Average, Mode based filtered average) of XL measure
		3) mean of internode count measure
		4) mean of coalescence rank measure
4) single_elem: can contain one of possible three values
		0: only one element of taxa_clust1 and one element of taxa_clust2 will be compared
		1: cluster containing taxa_clust1[0] and cluster containing taxa_clust2[0] will be compared
		2: All pairs of elements of taxa_clust1 and taxa_clust2 will be compared
		
@returns:
1) Average / max distance measure
2) 1 / 0 depending on whether the cluster pair is supported by at least one tree

"""
def FindAvgDistanceMeasure(taxa_clust1, taxa_clust2, DIST_MAT_TYPE, single_elem=2, type_of_output=0):
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
		#clust1 = Taxa_Info_Dict[taxa1]._Get_Taxa_Part_Clust_Idx()
		#clust2 = Taxa_Info_Dict[taxa2]._Get_Taxa_Part_Clust_Idx()
		taxa1_idx = COMPLETE_INPUT_TAXA_LIST.index(taxa1)
		taxa2_idx = COMPLETE_INPUT_TAXA_LIST.index(taxa2)
		clust1 = Taxa_Info_Dict[taxa1_idx]._Get_Taxa_Part_Clust_Idx()
		clust2 = Taxa_Info_Dict[taxa2_idx]._Get_Taxa_Part_Clust_Idx()
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
		x1_idx = COMPLETE_INPUT_TAXA_LIST.index(x1)
		for x2 in taxa_list2:  
			x2_idx = COMPLETE_INPUT_TAXA_LIST.index(x2)
			if (x1_idx < x2_idx):
				target_key = (x1_idx, x2_idx)
			elif (x1_idx > x2_idx):
				target_key = (x2_idx, x1_idx)
			else:
				continue
			if target_key in TaxaPair_Reln_Dict:
				if (DIST_MAT_TYPE == 1) or (DIST_MAT_TYPE == 2):
					val = TaxaPair_Reln_Dict[target_key]._GetNormalizedXLSumGeneTrees(DIST_MAT_TYPE)
				elif (DIST_MAT_TYPE == 3):
					val = TaxaPair_Reln_Dict[target_key]._GetAvgSumLevel()
				elif (DIST_MAT_TYPE == 4):
					#val = TaxaPair_Reln_Dict[target_key]._GetAvgLCARank()
					val = TaxaPair_Reln_Dict[target_key]._GetMaxLCARank()
				else:
					val = TaxaPair_Reln_Dict[target_key]._GetPseudoR4RelnRatio()
				"""
				append the value in the list
				"""
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

##----------------------------------------------------------------
#"""
#this function computes average XL information between a pair of taxa clusters
#@param: taxa_clust1: first taxa list
				#taxa_clust2: second taxa list
#"""
#def AvgR1RelnPriority(taxa_clust1, taxa_clust2):
	#score_val = 0
	#couplet_count = 0
	#for t1 in taxa_clust1:
		#t1_idx = COMPLETE_INPUT_TAXA_LIST.index(t1)
		#for t2 in taxa_clust2:
			#t2_idx = COMPLETE_INPUT_TAXA_LIST.index(t2)
			#"""
			#formation of the couplet key
			#"""
			#if (t1_idx < t2_idx):
				#target_key = (t1_idx, t2_idx)
				#target_reln = RELATION_R1
			#else:
				#target_key = (t2_idx, t1_idx)
				#target_reln = RELATION_R2
			
			#"""
			#accumulate the score
			#"""
			#if target_key in TaxaPair_Reln_Dict:
				#couplet_count = couplet_count + 1
				#pr_val = TaxaPair_Reln_Dict[target_key]._GetConnPrVal(target_reln)
				##ntree = TaxaPair_Reln_Dict[target_key]._GetNoSupportTrees()
				##freq_val = TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(target_reln)
				#score_val = score_val + (pr_val * 1.0) #/ ntree
	
	#if (couplet_count > 0):
		#return couplet_count, ((score_val * 1.0) / couplet_count)
		##return score_val
	
	#return 0, 0

#---------------------------------------------------------------------
"""
this function takes input of two sets of taxa (from two different taxa clusters)
and computes the relative frequency and support score measures
"""
def GetFreqScore_ClusterPair(taxa_list1, taxa_list2):
	"""
	for individual pair of clusters, first initialize the variables containing different statistics
	"""
	R1_freq = R2_freq = R3_freq = R4_freq = 0
	#R1_score = R2_score = R4_score = 0
	#R3_score = 0
	couplet_count = 0
	pseudo_r4_r1_count = pseudo_r4_r2_count = pseudo_r4_r3_count = 0
	
	#"""
	#these flags help to identify whether R1 / R2 relations are allowed between these pair of clusters
	#"""
	#r1_reln_allowed = True
	#r2_reln_allowed = True
	
	"""
	check individual couplets and accumulate the statistics
	"""
	for x1 in taxa_list1:
		x1_idx = COMPLETE_INPUT_TAXA_LIST.index(x1)
		for x2 in taxa_list2:  
			x2_idx = COMPLETE_INPUT_TAXA_LIST.index(x2)
			if (x1_idx < x2_idx):
				target_key = (x1_idx, x2_idx)
				compl_reln = False
			else:
				target_key = (x2_idx, x1_idx)
				compl_reln = True
			
			if target_key in TaxaPair_Reln_Dict:
				"""
				update the couplet count
				"""
				couplet_count = couplet_count + 1
				"""
				obtain the allowed relation list for this couplet
				"""
				allowed_reln_list = TaxaPair_Reln_Dict[target_key]._GetAllowedRelnList()
				
				"""
				accumulate the frequency and support score corresponding to the relation R3
				"""
				R3_freq = R3_freq + TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R3)
				#R3_score = R3_score + TaxaPair_Reln_Dict[target_key]._GetEdgeCost_ConnReln(RELATION_R3)

				"""
				accumulate the frequency and support score corresponding to the relation R4
				"""
				R4_freq = R4_freq + TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R4)
				#R4_score = R4_score + TaxaPair_Reln_Dict[target_key]._GetEdgeCost_ConnReln(RELATION_R4)

				"""
				two cases corresponding to the couplet indexed in the dictionary
				"""
				if (compl_reln == False):
					"""
					accumulate the frequency and support score corresponding to the relation R1
					"""
					R1_freq = R1_freq + TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R1)
					#R1_score = R1_score + TaxaPair_Reln_Dict[target_key]._GetEdgeCost_ConnReln(RELATION_R1)
					
					#if (RELATION_R1 not in allowed_reln_list) and (RELATION_R4 in allowed_reln_list):	# or RELATION_R2 in allowed_reln_list):
						#r1_reln_allowed = False

					"""
					accumulate the frequency and support score corresponding to the relation R2
					"""
					R2_freq = R2_freq + TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R2)
					#R2_score = R2_score + TaxaPair_Reln_Dict[target_key]._GetEdgeCost_ConnReln(RELATION_R2)

					#if (RELATION_R2 not in allowed_reln_list) and (RELATION_R4 in allowed_reln_list):	# or RELATION_R1 in allowed_reln_list):
						#r2_reln_allowed = False
					
					"""
					R4 relation and within it, pseudo R1, R2, and R3 relations
					"""
					pseudo_r4_r1_count = pseudo_r4_r1_count + TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(0)
					pseudo_r4_r2_count = pseudo_r4_r2_count + TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(1)
					pseudo_r4_r3_count = pseudo_r4_r3_count + TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(2)
					
				else:
					"""
					accumulate the frequency and support score corresponding to the relation R1
					here we check the score of the relation R2 
					since the couplet order is reversed
					"""
					R1_freq = R1_freq + TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R2)
					#R1_score = R1_score + TaxaPair_Reln_Dict[target_key]._GetEdgeCost_ConnReln(RELATION_R2)
					
					#if (RELATION_R2 not in allowed_reln_list) and (RELATION_R4 in allowed_reln_list):	# or RELATION_R1 in allowed_reln_list):
						#r1_reln_allowed = False
						
					"""
					accumulate the frequency and support score corresponding to the relation R2
					here we check the score of the relation R1
					since the couplet order is reversed
					"""
					R2_freq = R2_freq + TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(RELATION_R1)
					#R2_score = R2_score + TaxaPair_Reln_Dict[target_key]._GetEdgeCost_ConnReln(RELATION_R1)
						
					#if (RELATION_R1 not in allowed_reln_list) and (RELATION_R4 in allowed_reln_list):	# or RELATION_R2 in allowed_reln_list):
						#r2_reln_allowed = False

					"""
					R4 relation and within it, pseudo R1, R2, and R3 relations
					"""
					pseudo_r4_r1_count = pseudo_r4_r1_count + TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(1)
					pseudo_r4_r2_count = pseudo_r4_r2_count + TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(0)
					pseudo_r4_r3_count = pseudo_r4_r3_count + TaxaPair_Reln_Dict[target_key]._GetFreqPseudoR1(2)

	return R1_freq, R2_freq, R3_freq, R4_freq, couplet_count, pseudo_r4_r1_count, pseudo_r4_r2_count, pseudo_r4_r3_count

##----------------------------------------------------
#"""
#this function returns the frequency of 'inp_reln' from idx1 to idx2
#"""
#def GetRelnFreq(idx1, idx2, inp_reln):
	#if (idx1 < idx2):
		#target_key = (idx1, idx2)
		#if target_key in TaxaPair_Reln_Dict:
			#return TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(inp_reln)
	#else:
		#target_key = (idx2, idx1)
		#if target_key in TaxaPair_Reln_Dict:
			#TaxaPair_Reln_Dict[target_key]._GetEdgeWeight(Complementary_Reln(inp_reln))
	#"""
	#for unsupported couplets, we return frequency of -1
	#"""
	#return -1
	
#----------------------------------------------------
"""
this function checks whether the 'inp_reln' from idx1 to idx2 is a consensus relation or not
"""
def CheckConsensusReln(idx1, idx2, inp_reln):
	if (idx1 < idx2):
		target_key = (idx1, idx2)
		if target_key in TaxaPair_Reln_Dict:
			return TaxaPair_Reln_Dict[target_key]._CheckTargetRelnConsensus(inp_reln)
	else:
		target_key = (idx2, idx1)
		if target_key in TaxaPair_Reln_Dict:
			return TaxaPair_Reln_Dict[target_key]._CheckTargetRelnConsensus(Complementary_Reln(inp_reln))

	"""
	for unsupported couplets, we return False
	"""
	return False
	
##----------------------------------------------------
#"""
#this function checks whether the 'inp_reln' from idx1 to idx2 is an allowed relation or not
#"""
#def CheckAllowedReln(idx1, idx2, inp_reln):
	#if (idx1 < idx2):
		#target_key = (idx1, idx2)
		#if target_key in TaxaPair_Reln_Dict:
			#return TaxaPair_Reln_Dict[target_key]._CheckInpRelnAllowed(inp_reln)
	#else:
		#target_key = (idx2, idx1)
		#if target_key in TaxaPair_Reln_Dict:
			#return TaxaPair_Reln_Dict[target_key]._CheckInpRelnAllowed(Complementary_Reln(inp_reln))
	
	#return False
	
##----------------------------------------------------
#"""
#this function checks whether the 'inp_reln' from idx1 to idx2 is an allowed relation and also the 
#only allowed relation or not
#"""
#def CheckSingleAllowedReln(idx1, idx2, inp_reln):
	#if (idx1 < idx2):
		#target_key = (idx1, idx2)
		#if target_key in TaxaPair_Reln_Dict:
			#return TaxaPair_Reln_Dict[target_key]._Check_Single_Reln_OnlyAllowed(inp_reln)
	#else:
		#target_key = (idx2, idx1)
		#if target_key in TaxaPair_Reln_Dict:
			#return TaxaPair_Reln_Dict[target_key]._Check_Single_Reln_OnlyAllowed(Complementary_Reln(inp_reln))
	
	#return False

#----------------------------------------------------------
"""
this function checks whether the 'inp_reln' from clust1 to clust2 is an allowed relation and also the 
only allowed relation or not
"""
def CheckAllowedRelnClusterPair(clust1, clust2, inp_reln):
	if (clust1 < clust2):
		target_key = (clust1, clust2)
		if target_key in Cluster_Pair_Info_Dict:
			return Cluster_Pair_Info_Dict[target_key]._CheckInpRelnAllowed(inp_reln)
	else:
		target_key = (clust2, clust1)
		if target_key in Cluster_Pair_Info_Dict:
			return Cluster_Pair_Info_Dict[target_key]._CheckInpRelnAllowed(Complementary_Reln(inp_reln))

#----------------------------------------------------------
"""
this function returns the allowed relation list between this cluster pair
"""
def GetAllowedRelnClusterPair(clust1, clust2):
	if (clust1 < clust2):
		target_key = (clust1, clust2)
		if target_key in Cluster_Pair_Info_Dict:
			return Cluster_Pair_Info_Dict[target_key]._GetPossibleRelnList()
	else:
		target_key = (clust2, clust1)
		if target_key in Cluster_Pair_Info_Dict:
			return Cluster_Pair_Info_Dict[target_key]._GetPossibleRelnList()
	
	return []
	
#----------------------------------------------------
"""
this function returns the frequency of 'inp_reln' from clust1 to clust2
"""
def GetRelnFreq_ClusterPair(clust1, clust2, inp_reln):
	if (clust1 < clust2):
		target_key = (clust1, clust2)
		if target_key in Cluster_Pair_Info_Dict:
			return Cluster_Pair_Info_Dict[target_key]._GetFreq(inp_reln)
	else:
		target_key = (clust2, clust1)
		if target_key in Cluster_Pair_Info_Dict:
			return Cluster_Pair_Info_Dict[target_key]._GetFreq(Complementary_Reln(inp_reln))
	
	"""
	for unsupported cluster pairs, return 0
	"""
	return 0

##----------------------------------------------------
#"""
#this function returns the frequency of 'inp_reln' from clust1 to clust2
#"""
#def GetRelnTransitiveFreq_ClusterPair(clust1, clust2, inp_reln):
	#if (clust1 < clust2):
		#target_key = (clust1, clust2)
		#if target_key in Cluster_Pair_Info_Dict:
			#return Cluster_Pair_Info_Dict[target_key]._GetTransitiveFreq(inp_reln)
	#else:
		#target_key = (clust2, clust1)
		#if target_key in Cluster_Pair_Info_Dict:
			#return Cluster_Pair_Info_Dict[target_key]._GetTransitiveFreq(Complementary_Reln(inp_reln))
	
	#"""
	#for unsupported cluster pairs, return 0
	#"""
	#return 0
	
#----------------------------------------------------
"""
this function returns the frequency of 'inp_reln' from clust1 to clust2
"""
def GetRelnScore_ClusterPair(clust1, clust2, inp_reln):
	if (clust1 < clust2):
		target_key = (clust1, clust2)
		if target_key in Cluster_Pair_Info_Dict:
			return Cluster_Pair_Info_Dict[target_key]._GetSupportScore(inp_reln)
	else:
		target_key = (clust2, clust1)
		if target_key in Cluster_Pair_Info_Dict:
			return Cluster_Pair_Info_Dict[target_key]._GetSupportScore(Complementary_Reln(inp_reln))
	
	"""
	for unsupported cluster pairs, return 0
	"""
	return 0

##-----------------------------------------------------------------
#"""
#this function checks whether there is an existing relation from "clust1" to "clust2"
#"""
#def Check_Existing_Relation(ReachMat, clust1, clust2, inp_reln):
	#clust1_reachmat_idx = CURRENT_CLUST_IDX_LIST.index(clust1)
	#clust2_reachmat_idx = CURRENT_CLUST_IDX_LIST.index(clust2)
	#if (inp_reln == RELATION_R1):
		#if (ReachMat[clust1_reachmat_idx][clust2_reachmat_idx] == 1):
			#return True
	#if (inp_reln == RELATION_R2):
		#if (ReachMat[clust2_reachmat_idx][clust1_reachmat_idx] == 1):
			#return True
	#if (inp_reln == RELATION_R4):
		#if (ReachMat[clust2_reachmat_idx][clust1_reachmat_idx] == 2) and (ReachMat[clust1_reachmat_idx][clust2_reachmat_idx] == 2):
			#return True
		
	#return False

#-------------------------------------------------------------
"""
this function checks whether there exists a pair of cluster (at least one relation between them 
with respect to the input trees)
"""
def ClusterPairExists(clust1, clust2):
	if (clust1 < clust2):
		target_key = (clust1, clust2)
	else:
		target_key = (clust2, clust1)
	
	if target_key in Cluster_Pair_Info_Dict:
		return True
	
	return False
		
#-------------------------------------------------------------
"""
this function checks whether the R1 relation clust list and R2 relation clust list of "clust2" 
is subsets of corresponding lists in "clust1"
"""
def CheckSubsetClust(clust1, clust2):
	clust1_R1 = Cluster_Info_Dict[clust1]._GetClustRelnList(RELATION_R1)
	clust1_R2 = Cluster_Info_Dict[clust1]._GetClustRelnList(RELATION_R2)
	clust2_R1 = Cluster_Info_Dict[clust2]._GetClustRelnList(RELATION_R1)
	clust2_R2 = Cluster_Info_Dict[clust2]._GetClustRelnList(RELATION_R2)
	if (set(clust2_R1).issubset(set(clust1_R1)) == True):
		if (set(clust2_R2).issubset(set(clust1_R2)) == True):	#(len(clust2_R2) == 0):
			return True
	return False
	
#---------------------------------------------
"""
the function is used for checking equality of two floating point numbers
"""
def FlEq(a, b, eps=0.000001):
	#return (abs(math.log( a ) - math.log(b)) <= eps)
	return (abs(a - b) <= eps)

#-----------------------------------------------------
# this function finds the MRCA of this two input taxa labels
# this is a custom function
# without using standard dendropy routine
def Find_MRCA(Inp_Tree, spec_list):
	node1 = Inp_Tree.find_node_with_taxon_label(spec_list[0])
	pn = node1.parent_node
	while (pn is not None):
		leaf_labels = []
		for n in pn.leaf_nodes():
			leaf_labels.append(n.taxon.label)
		if set(spec_list).issubset(set(leaf_labels)):
			return pn
		pn = pn.parent_node
			
	return None

	
	
	
