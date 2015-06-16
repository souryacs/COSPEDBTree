#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#------------------------------------------------------
""" this function sorts the input list containing the costs of different relations between individual taxa pairs
we don't use python standard sort routine
rather we use sorting giving different weightage to different edge types, in case of identical cost values """
def Sort_Cost_List_Initial(Cost_List_Node_Pair):
  Sort_Priority_Queue(Cost_List_Node_Pair)
  
#------------------------------------------------------
""" this code section implements the max priority queue """
#------------------------------------------------------
# parent node of current node
def Parent(idx):
  return int((idx-1) / 2)

# left child node of current node  
def Left(idx):
  return (2*idx+1)

# right child node of current node
def Right(idx):
  return (2*idx+2)
  
#------------------------------------------------------------------------
# version 4
#-------------------------
"""
def Lower_Score_Value(inp_queue, i, j):
  #---------------------
  # two cases 
  #---------------------
  # case 1 - when both of the scores are non-positive
  if (inp_queue[i][3] <= 0) and (inp_queue[j][3] <= 0):
    if (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
      return 1      
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])):
      return 0      
    elif (inp_queue[i][3] < inp_queue[j][3]):
      return 1
    elif (inp_queue[i][3] > inp_queue[j][3]):
      return 0
    elif ((inp_queue[i][2] == NO_EDGE) and (inp_queue[j][2] != NO_EDGE)):
      # no edge has high priority
      return 0
    elif ((inp_queue[j][2] == NO_EDGE) and (inp_queue[i][2] != NO_EDGE)):
      return 1
    elif ((inp_queue[i][2] == BI_DIRECTED_EDGE) and ((inp_queue[j][2] == DIRECTED_OUT_EDGE) \
      or (inp_queue[j][2] == DIRECTED_IN_EDGE))):
      return 1      
    elif ((inp_queue[j][2] == BI_DIRECTED_EDGE) and ((inp_queue[i][2] == DIRECTED_OUT_EDGE) \
      or (inp_queue[i][2] == DIRECTED_IN_EDGE))):
      return 0
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
      return 1
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])):
      return 0    
  else:
    # case 2 - otherwise we follow previous algorithm 
    if (inp_queue[i][3] < inp_queue[j][3]):
      return 1
    elif (inp_queue[i][3] == inp_queue[j][3]):
      # for tie case of cost
      if ((inp_queue[i][2] == NO_EDGE) and (inp_queue[j][2] != NO_EDGE)):
	# change - sourya
	# no edge has high priority
	return 0
      elif ((inp_queue[j][2] == NO_EDGE) and (inp_queue[i][2] != NO_EDGE)):
	return 1
      # end change - sourya
      elif ((inp_queue[i][2] == BI_DIRECTED_EDGE) and ((inp_queue[j][2] == DIRECTED_OUT_EDGE) or (inp_queue[j][2] == DIRECTED_IN_EDGE))):
	# bi directed edge is set to low priority, compared to other directed edges - add - sourya
	return 1      
      elif ((inp_queue[j][2] == BI_DIRECTED_EDGE) and ((inp_queue[i][2] == DIRECTED_OUT_EDGE) or (inp_queue[i][2] == DIRECTED_IN_EDGE))):
	return 0
      elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	    < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
	# higher edge priority have greater priority
	return 1
      elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])\
	    < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])):
	return 0    
      elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	    < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
	# higher edge frequency have greater priority - add - sourya
	return 1      
      elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])\
	    < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])):
	return 0      
      
  return 0
"""
#------------------------------------------------------------------------
# version 1 
#-------------------------

# this function is used for cost based sorting of the inp_queue in the unweighted supertree
# this function compares two elements of the heap
# returns 1 if element in the i'th index is lower than element in the j'th index

def Lower_Score_Value(inp_queue, i, j):
  if (inp_queue[i][3] < inp_queue[j][3]):
    return 1
  elif (inp_queue[i][3] > inp_queue[j][3]):
    return 0
  else:	#if (inp_queue[i][3] == inp_queue[j][3]):
    # for tie case of cost
    if ((inp_queue[i][2] == NO_EDGE) and (inp_queue[j][2] != NO_EDGE)):
      # no edge has low priority
      return 1
    elif ((inp_queue[j][2] == NO_EDGE) and (inp_queue[i][2] != NO_EDGE)):	#added
      return 0
    elif ((inp_queue[i][2] == BI_DIRECTED_EDGE) and ((inp_queue[j][2] == DIRECTED_OUT_EDGE) or (inp_queue[j][2] == DIRECTED_IN_EDGE))):
      # bi directed edge is set to low priority, compared to other directed edges - add - sourya
      return 1      
    elif ((inp_queue[j][2] == BI_DIRECTED_EDGE) and ((inp_queue[i][2] == DIRECTED_OUT_EDGE) or (inp_queue[i][2] == DIRECTED_IN_EDGE))):	#added
      return 0
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
      # higher edge priority have greater priority
      return 1
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])):	#added
      return 0    
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
      # higher edge frequency have greater priority - add - sourya
      return 1      
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])):	#added
      return 0      
    
  return 0

#------------------------------------------------------------------------
# old version
"""
# this function compares two elements of the heap
# returns 1 if element in the i'th index is lower than element in the j'th index
def Lower_Score_Value(inp_queue, i, j):
  if (inp_queue[i][3] < inp_queue[j][3]):
    return 1
  elif (inp_queue[i][3] == inp_queue[j][3]):
    # for tie case of cost
    if ((inp_queue[i][2] == NO_EDGE) and (inp_queue[j][2] != NO_EDGE)):
      # no edge has low priority
      return 1
    elif ((inp_queue[i][2] == BI_DIRECTED_EDGE) and ((inp_queue[j][2] == DIRECTED_OUT_EDGE) or (inp_queue[j][2] == DIRECTED_IN_EDGE))):	# add - sourya
      # bi directed edge is set to low priority, compared to other directed edges - add - sourya
      return 1      
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
      # higher edge priority have greater priority
      return 1
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
      # higher edge frequency have greater priority - add - sourya
      return 1      
  return 0
"""
#------------------------------------------------------------------------
  
# this function exchanges two elements in the heap
def Exchange_Elem(inp_queue, i, j):
  temp_key = inp_queue[i][0]
  inp_queue[i][0] = inp_queue[j][0]
  inp_queue[j][0] = temp_key
  temp_key = inp_queue[i][1]
  inp_queue[i][1] = inp_queue[j][1]
  inp_queue[j][1] = temp_key
  temp_edge = inp_queue[i][2]
  inp_queue[i][2] = inp_queue[j][2]
  inp_queue[j][2] = temp_edge
  temp_val = inp_queue[i][3]
  inp_queue[i][3] = inp_queue[j][3]
  inp_queue[j][3] = temp_val

# maintain max heap property
# note: heap_size may not be the actual length of the queue
# but the working length (on which the remaining sorting operation will take place)
def Max_Heapify(inp_queue, idx, heap_size):
  l = Left(idx)
  r = Right(idx)
  if (l < heap_size) and Lower_Score_Value(inp_queue, idx, l):
    largest_idx = l
  else:
    largest_idx = idx
  if (r < heap_size) and Lower_Score_Value(inp_queue, largest_idx, r):
    largest_idx = r
  if (largest_idx != idx):
    Exchange_Elem(inp_queue, idx, largest_idx)
    Max_Heapify(inp_queue, largest_idx, heap_size)

# extract the current maximum and also pop the element from the heap
def Heap_Extract_Max(inp_queue):
  if (len(inp_queue) < 1):
    print 'underflow of max priority queue'
  # 1st element is the maximum
  # comment - sourya
  # this only assigns the pointer which will lead incorrect result
  # when the heap is modified
  # if the 
  #max_elem = inp_queue[0]
  # add - sourya
  max_elem = list(inp_queue[0])
  # end add - sourya
  # replace the first element with the last element of the queue
  Exchange_Elem(inp_queue, 0, len(inp_queue) - 1)
  # delete the last element of the queue
  del inp_queue[len(inp_queue) - 1]
  heap_size = len(inp_queue)
  # call the max_heapify function to maintain the heap property
  # 0 is the starting index of the list storing the heap structure
  Max_Heapify(inp_queue, 0, heap_size)
  return max_elem
  
# this function builds the priority queue (max heap property)
def Build_Max_Heap(inp_queue):
  heap_size = len(inp_queue)
  for idx in range(int(len(inp_queue) / 2), -1, -1):
    Max_Heapify(inp_queue, idx, heap_size)
    
# this is the heap sort algorithm
def Sort_Priority_Queue(inp_queue):
  Build_Max_Heap(inp_queue)
  heap_size = len(inp_queue)
  ## should be commented - sourya
  #for idx in range((len(inp_queue) - 1), 0, -1):
    #Exchange_Elem(inp_queue, 0, idx)
    #heap_size = heap_size - 1
    #Max_Heapify(inp_queue, idx, heap_size)
  
def Queue_Increase_Score(inp_queue, idx, target_val):
  inp_queue[idx][3] = target_val
  while (idx > 0) and Lower_Score_Value(inp_queue, Parent(idx), idx):
    Exchange_Elem(inp_queue, idx, Parent(idx))
    idx = Parent(idx)
    
def Queue_Decrease_Score(inp_queue, idx, target_val):
  inp_queue[idx][3] = target_val
  Max_Heapify(inp_queue, idx, len(inp_queue))
      