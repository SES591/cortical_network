#!/usr/bin/python
#bioinfo.py

## this is what we'll run the most

__author__ = '''Hyunju Kim/Harrison Smith'''

import os
import sys
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

import input_net as inet
import updating_rule as ur

## Global Params ##
n_states = 2 # number of states for a node (same as Nbr_States)

###################


#########################################
def decimal_to_binary(nodes_list, network_ID): # more left in the nodes list means higher order of 2 in binary
    """:
        Converts network_ID (eg 1024) 
        to 
        network_bID (eg. {'pS': 1, 'pP': 1, 'gS': 1, 'gF': 1, 'gP': 1, 'pC': 1, 'gE': 1, 'pF': 1, 'pE': 1, 'gC': 1})

        Arguments:
            nodes_list [list of strs]
                list of nodes in the network

            network_ID [int]
                network_ID specifiying unique initial network configuration

        Returns:
            network_bID [dict of booleans]
                Dictionary of key=node, value= 0 or 1 specifying node state
    """
    network_bID = {}
    x = len(nodes_list) -1 # x==len(n_nodes)-1? why is this the number of nodes - 1?
    for node in nodes_list:
        network_bID[node] = network_ID / (n_states**x)
        network_ID = network_ID % (n_states**x)
        x = x - 1

    return network_bID

#########################################
def binary_to_decimal(nodes_list, network_bID):  # more left in the nodes list means higher order of 2 in binary
    """:
        Converts network_bID (eg. {'pS': 1, 'pP': 1, 'gS': 1, 'gF': 1, 'gP': 1, 'pC': 1, 'gE': 1, 'pF': 1, 'pE': 1, 'gC': 1})
        to 
        network_ID (eg 1024) 

        Arguments:
            nodes_list [list of strs]
                list of nodes in the network

            network_ID [int]
                network_ID specifiying unique initial network configuration

        Returns:
            network_bID [dict of booleans]
                Dictionary of key=node, value= 0 or 1 specifying node state
    """
    network_ID = 0
    x = len(nodes_list) -1
    for node in nodes_list:
        network_ID = network_ID + network_bID[node] * (n_states**x)
        x = x - 1
    return network_ID

#########################################
def ensemble_time_series(G, nodes_list, timesteps=20):
    ''':
        Applies boolean updating rule to the network nodes for the specified number of timesteps.
        AND
        Applies boolean updating rule to the network nodes accross ALL possible initial network states.

        Arguments:
            G [networkx Graph object]

            nodes_list [list of strs]
                list of nodes in the network

            timesteps [int]
                number of timesteps for network to evolve
        
        Return:
            timeseriesdata [dict of (dict of list of ints)]
                ex. {'gC': {0: [0, 1, 0, 0, 1, 0, 0, 1, 0, 0], 1: [1, 0, 0, 1, 0, 0, 1, 0, 0, 1],...}, 'gF': {...}} 
                innermost list (ex. [0, 1, 0, 0, 1, 0, 0, 1, 0, 0]) - timeseries of states of the node
                innermost key (ex. 0)- network_ID specificing the initial state of all nodes in the network
                outermost key (ex. 'gC')- node_ID specifing the node whose states are being listed
    '''
    
    n_nodes = len(G.nodes()) # number of nodes in network
    n_init_networks = n_states**n_nodes # number of possible initial network states
    
    timeseriesdata = {}
    ## Create empty time series list for every combination 
    ## of node and initial network state
    for node in G.nodes():
        timeseriesdata[node] = {}
        for network_ID in range(0, n_init_networks):
            timeseriesdata[node][network_ID] = [] 
    
    ## Fill time series list with individual node states for 
    ## every combination of node and initial network state
    for network_ID in range(0, n_init_networks):
        network_bID = decimal_to_binary(nodes_list, network_ID)
        for step in range(0, timesteps):
            prev_network_bID = network_bID.copy()
            for node in nodes_list:
                timeseriesdata[node][network_ID].append(prev_network_bID[node])
            network_bID = ur.boolean_updating(G, prev_network_bID)

    return timeseriesdata
#########################################
def net_state_transition(G, nodes_list):
    ''':
        Arguments:
            G [networkx Graph object]
            
            nodes_list [list of strs]
                list of nodes in the network
        Return:
            G_transition_graph [networkx Graph object]
                networkx graph directed graph showing how network configurations map to one another
    '''
    
    n_nodes = len(G.nodes())
    n_init_networks = n_states**n_nodes
    G_transition_graph = nx.DiGraph()

    for prev_network_ID in range(0, n_init_networks):

        prev_network_bID = decimal_to_binary(nodes_list, prev_network_ID)
        network_bID = ur.boolean_updating(G, prev_network_bID)
        network_ID = binary_to_decimal(nodes_list, network_bID)
        G_transition_graph.add_edge(prev_network_ID, network_ID) #this is a directed graph showing you attractor landscape for how each of the states trasitions to each other
    
    return G_transition_graph
#########################################
def find_attractor(G_transition_graph):
    
    ''':
        Arguments:
            G_transition_graph [networkx Graph object]
                networkx graph directed graph showing how network configurations map to one another
        Return:
            attractors [dict of list of lists of ints]
                ['fixed'] = [[532][948]]
                ['cycle'] = []
    '''
    attractor_list = nx.simple_cycles(G_transition_graph) #in case of deterministic system, any cycle without considering edge direction will be directed cycle.
    attractors = {}
    attractors['fixed'] = []
    attractors['cycle'] = []

    for network_ID in attractor_list:
        if len(network_ID) == 1:
            attractors['fixed'].append(network_ID)
        else:
            attractors['cycle'].append(network_ID)

    return attractors #this outputs decID of attractor states (fixed and cyclic)
#########################################
def main():
    print "time_evol module is the main code."
    edge_file = '../data/inputs/edges-init.dat' 
    node_file = '../data/inputs/ant-nodes-init.dat'
    
    G = inet.read_network_from_file(edge_file, node_file)
    nodes_list = G.nodes() # print graph nodes without states. Could also do nx.nodes(G)

    print G.nodes()

    #####

    timeseriesdata = ensemble_time_series(G, nodes_list, 10)#, timesteps=20)

    network_ID = 1 #this is the ID of the init state
    biStates = decimal_to_binary(nodes_list, network_ID)
    print 'initial state', biStates
    for node in G.nodes():
        print node, timeseriesdata[node][1]

    G_transition_graph = net_state_transition(G, nodes_list)
    # nx.draw(G_transition_graph)
    # plt.show()

    attractors = find_attractor(G_transition_graph)
    print attractors





if __name__=='__main__':
    main()
