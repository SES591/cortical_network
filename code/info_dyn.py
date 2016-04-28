#!/usr/bin/python
#bionetworks.py
#last update : 14 Aug 2014

__author__ = '''Hyunju Kim'''


import networkx as nx
import os
import sys
import random as ran
from math import log
from optparse import OptionParser, OptionGroup
from scipy import *
import numpy as np
from collections import defaultdict

from time_evol import n_states

#########################################################################
# def build_historyList_new(node_timeseries, historylength):
#     n_slices = len(node_timeseries) - historylength + 1
#     historylist = []
#     for i in range(n_slices):
#         list_slice = node_timeseries[i:(historylength+i)]
#         joined_slice = ''.join(map(str, list_slice))
#         dec_slice = int(joined_slice, 2)
#         historylist.append(dec_slice)
#     return historylist


#########################################################################
def build_historyList(node_timeseries, historylength):
    '''

    Description:
        -- To generate a list of decimal states for k-previous states of a node with a given list of dynamical states of the nodes
     
    Arguments:
        -- 1. node_timeseries (used to be aList)= a given list of dynamical states of a node (this is the list of the timeseries of a single node from a single initital state)
        -- 2. historylength = the number of previous states converted to a decimal state
        -- 3. n_states = the number of possible states for each node (2 by default)
          
    Return:
        -- a list of decimal states for k-previous states of a node
    '''
    # print "node_timeseries", node_timeseries

    historylist = []
    # print historylist
    historyunit = node_timeseries[:historylength] 
    # print "hist-unit",historyunit
    node_timeseries[:historylength] = [] # remove length of history from front of node_timeseries
    # print node_timeseries
    #print "node_timeseries - Unit", node_timeseries

    # exit()
    
    historystate = 0
    for s in range(historylength):
        historystate += historyunit[s] * (n_states**s)
        # print "hist-state", historystate


    #print "unit", historyunit[s]
    #print historystate

    historylist.append(historystate)

    for x in node_timeseries:
        historystate = historystate / n_states + x * np.power(n_states, historylength - 1)
        historylist.append(historystate)
        #print x, historystate


        
    return historylist
#########################################################################
def compute_AI(timeSeriesNode, historyLength, Nbr_Initial_States):
    '''
    Description:
    -- compute AI for every node using distribution from all possible initial conditions or an arbitrary set of initial conditions
    
    Arguments:
        timeSeriesNode -- the time series for a single node (so this is a dict)
    Note:
    -- for an arbitrary set of initial conditions, one need to specify the number of initial conditions, Nbr_Initial_States, and name every initial condition from 0 to Nbr_Initial_States - 1
    '''
    count_currState_hiState = defaultdict(int)
    count_hiState = defaultdict(int)
    count_currState = defaultdict(int)


    # for si in range(Nbr_Initial_States):
    for si in timeSeriesNode:
        aList = list(timeSeriesNode[si])
        historyList = build_historyList(aList, historyLength) # aList becomes aList[historyLength:] after historyList function

        

        #print "test", len(aList)
        #*** To obtain the distribution for each pattern ***#

        for s in range(len(aList)):
            count_currState_hiState[(aList[s], historyList[s])] += 1
            count_hiState[historyList[s]] += 1
            count_currState[aList[s]] += 1

    #        print count_currState_hiState
    #        print count_currState
    #        print count_hiState
    AI = 0
    for si in timeSeriesNode:
        aList = list(timeSeriesNode[si])
        historyList = build_historyList(aList, historyLength) # aList becomes aList[historyLength:] after historyList function


        #print "after counting aList", aList
        #print "after counting history", historyList

        
        sampleLength = len(aList) * Nbr_Initial_States
        for s in range(len(aList)):
            prob_currState_hiState = float(count_currState_hiState[(aList[s], historyList[s])]) / float(sampleLength)
            #print "joint", prob_currState_hiState
            prob_hiState = float(count_hiState[historyList[s]]) / float(sampleLength)
            #print "his", prob_hiState, historyList[s]
            prob_currState = float(count_currState[aList[s]]) / float(sampleLength)
            #print "curr", prob_currState
            AI = AI + log( prob_currState_hiState / ( prob_currState * prob_hiState)) / log(2.0) # since the summation is over not all possible pattern of currState_hiState
    AI = AI / float(sampleLength)
    return AI

#########################################################################
def compute_TE(timeSeriesNodeA, timeSeriesNodeB, historyLength, Nbr_Initial_States):

    '''
    Description:
    -- compute TE for every pair of nodes using distribution from all possible initial conditions or an arbitrary set of initial conditions
    
    Note:
    -- for an arbitrary set of initial conditions, one need to specify the number of initial conditions, Nbr_Initial_States, and name every initial condition from 0 to Nbr_Initial_States - 1
    '''

    

    #*** declare dic for distribution to compute Transfer Entropy ***#
    count_tarCurrState_tarHiState_sourPrevState = defaultdict(int)
    count_tarCurrState_tarHiState = defaultdict(int)
    count_tarHiState_sourPrevState = defaultdict(int)
    count_tarHiState = defaultdict(int)
    
    for si in timeSeriesNodeA: # doesn't matter if i'm looping through keys of timeSeriesNodeA or B because they are the same
        sourList = list(timeSeriesNodeA[si])
        tarList = list(timeSeriesNodeB[si])
        historyList = build_historyList(tarList, historyLength) # tarList becomes tarList[historyLength:] after historyList function
        sourList[:historyLength - 1] = [] # sourList becomes sourList[historyLength-1:]
        #*** To obtain the distribution for each pattern to compute Transfer Entropy ***#
        for s in range(len(tarList)):
            count_tarCurrState_tarHiState_sourPrevState[(tarList[s], historyList[s], sourList[s])] += 1
            count_tarCurrState_tarHiState[(tarList[s], historyList[s])] += 1
            count_tarHiState_sourPrevState[(historyList[s], sourList[s])] += 1
            count_tarHiState[historyList[s]] += 1

    #*** obtain the distribution for each pattern to compute Active Information ***#
    TE = 0
    for si in timeSeriesNodeA:
        sourList = list(timeSeriesNodeA[si])
        tarList = list(timeSeriesNodeB[si])
        historyList = build_historyList(tarList, historyLength) # tarList becomes tarList[historyLength:] after historyList function
        sourList[:historyLength - 1] = [] # sourList becomes sourList[historyLength-1:]
        sampleLength = len(tarList) * Nbr_Initial_States
        for s in range(len(tarList)):
            prob_tarCurrState_tarHiState_sourPrevState = float(count_tarCurrState_tarHiState_sourPrevState[(tarList[s], historyList[s], sourList[s])]) / float(sampleLength)
            prob_tarCurrState_tarHiState = float(count_tarCurrState_tarHiState[(tarList[s], historyList[s])]) / float(sampleLength)
            prob_tarHiState_sourPrevState = float(count_tarHiState_sourPrevState[(historyList[s], sourList[s])]) / float(sampleLength)
            prob_tarHiState = float(count_tarHiState[historyList[s]]) / float(sampleLength)
            TE = TE + log( (prob_tarCurrState_tarHiState_sourPrevState * prob_tarHiState) / ( prob_tarHiState_sourPrevState * prob_tarCurrState_tarHiState)) / log(2.0) # since the summation is over not all possible pattern of tarCurrState_tarHiState_sourPrevState but tarList, there is no prob_tarCurrState_tarHiState_sourPrevState multiplied by the log term.
    TE = TE / float(sampleLength)
    return TE

#########################################################################
def main():

    # timeSeriesNodeA = {0: [0,0,1, 0, 0, 1, 1, 1]}#, 1:[0,0,1, 0, 0, 0, 1, 0] }
    # timeSeriesNodeB = {0: [0,0,1, 0, 1, 0, 0, 1]}#, 1:[1,1,1, 0, 0, 0, 1, 1] }

    timeSeriesNodeA = {0: [0,0,1, 0, 0, 1, 1, 1], 1:[0,0,1, 0, 0, 0, 1, 0] }
    historylength = 5

    n_slices = len(timeSeriesNodeA[0]) - historylength + 1
    for i in range(n_slices):
        list_slice = timeSeriesNodeA[0][i:(historylength+i)]
        joined_slice = ''.join(map(str, list_slice))
        print joined_slice
        dec_slice = int(joined_slice, 2)
        print dec_slice
        # print int('00001',2)

    print build_historyList(timeSeriesNodeA[0], historylength)

if __name__=='__main__':
    main()

#print compute_local_TE(timeSeriesNodeA[0], timeSeriesNodeB[1], 2, 1,2)

#########################################################################
#timeSeriesNode = {0: [0,0,1, 0, 0, 1, 1, 1], 1:[0,0,1, 0, 0, 0, 1, 0] }


#hList = build_historyList(aList, 2)
#print "history list in decimal states", hList

#print compute_AI(timeSeriesNode, 2, 2)

#tarList = [1,1,1,0,1,1,0,1,0,0,1,1,0,0,1,1,1,0,1,1,0,1]
#sourList = [1,0,1,0,1,0,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,1]
#print comTE(sourList, tarList, 1)

