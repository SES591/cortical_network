#!/usr/bin/python
#bioinfo.py

## this is what we'll run the most

__author__ = '''Hyunju Kim/Harrison Smith'''

import os
import sys
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import OrderedDict

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

            network_bID [dict of booleans]
                Dictionary of key=node, value= 0 or 1 specifying node state

        Returns:
            network_ID [int]
                network_ID specifiying unique initial network configuration
    """
    network_ID = 0
    x = len(nodes_list) -1
    for node in nodes_list:
        network_ID = network_ID + network_bID[node] * (n_states**x)
        x = x - 1
    return network_ID

#########################################
def time_series_all(G, nodes_list, timesteps=20):
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
                networkx graph directed graph showing how network configurations map to one another in a deterministic system
                XXXTHIS TRANSITION GRAPH IS DEPENDENT ON THE STARTING STATE OF G BECAUSE IT ONLY MAPS THE STATE TRANSTITIONS
                XXXFROM THAT STARTING STATE. IT SHOULD BE SHOWING THE ATTRACTOR LANDSCAPE REGARDLESS OF YOUR STARTING POSITION.
                XXXI THINK IT IS MESSING UP BECAUSE I AM PASSING G THROUGH THE BOOLEAN UPDATING FUNCTION.
                4/13/16 4:44pm- FIXED. now this should show the complete transition graph. Still need to check to make sure
                that the attractors being output are accessible from the specified initial state though.

    '''
    
    n_nodes = len(G.nodes())
    n_init_networks = n_states**n_nodes
    G_transition_graph = nx.DiGraph()

    for prev_network_ID in range(0, n_init_networks):

        prev_network_bID = decimal_to_binary(nodes_list, prev_network_ID)
        network_bID = ur.boolean_updating(G, prev_state=prev_network_bID)
        network_ID = binary_to_decimal(nodes_list, network_bID)
        G_transition_graph.add_edge(prev_network_ID, network_ID) #this is a directed graph showing you attractor landscape for how each of the states trasitions to each other
    
    return G_transition_graph
#########################################
def control_state_1_transition(G, nodes_list, control_node):
    ''':
        Arguments:
            G [networkx Graph object]
            
            nodes_list [list of strs]
                list of nodes in the network
        Return:
            G_transition_graph [networkx Graph object]
                networkx graph directed graph showing how network configurations map to one another in a deterministic system
                XXXTHIS TRANSITION GRAPH IS DEPENDENT ON THE STARTING STATE OF G BECAUSE IT ONLY MAPS THE STATE TRANSTITIONS
                XXXFROM THAT STARTING STATE. IT SHOULD BE SHOWING THE ATTRACTOR LANDSCAPE REGARDLESS OF YOUR STARTING POSITION.
                XXXI THINK IT IS MESSING UP BECAUSE I AM PASSING G THROUGH THE BOOLEAN UPDATING FUNCTION.
                4/13/16 4:44pm- FIXED. now this should show the complete transition graph. Still need to check to make sure
                that the attractors being output are accessible from the specified initial state though.

    '''
    
    n_nodes = len(G.nodes())
    n_init_networks = n_states**n_nodes
    G_transition_graph = nx.DiGraph()

    for prev_network_ID in range(0, n_init_networks):

        prev_network_bID = decimal_to_binary(nodes_list, prev_network_ID)
        network_bID = ur.control_1_updating(G, prev_state=prev_network_bID, control_node=control_node)
        network_ID = binary_to_decimal(nodes_list, network_bID)
        G_transition_graph.add_edge(prev_network_ID, network_ID) #this is a directed graph showing you attractor landscape for how each of the states trasitions to each other
    
    return G_transition_graph
#########################################
def control_state_0_transition(G, nodes_list, control_node):
    ''':
        Arguments:
            G [networkx Graph object]
            
            nodes_list [list of strs]
                list of nodes in the network
        Return:
            G_transition_graph [networkx Graph object]
                networkx graph directed graph showing how network configurations map to one another in a deterministic system
                XXXTHIS TRANSITION GRAPH IS DEPENDENT ON THE STARTING STATE OF G BECAUSE IT ONLY MAPS THE STATE TRANSTITIONS
                XXXFROM THAT STARTING STATE. IT SHOULD BE SHOWING THE ATTRACTOR LANDSCAPE REGARDLESS OF YOUR STARTING POSITION.
                XXXI THINK IT IS MESSING UP BECAUSE I AM PASSING G THROUGH THE BOOLEAN UPDATING FUNCTION.
                4/13/16 4:44pm- FIXED. now this should show the complete transition graph. Still need to check to make sure
                that the attractors being output are accessible from the specified initial state though.

    '''
    
    n_nodes = len(G.nodes())
    n_init_networks = n_states**n_nodes
    G_transition_graph = nx.DiGraph()

    for prev_network_ID in range(0, n_init_networks):

        prev_network_bID = decimal_to_binary(nodes_list, prev_network_ID)
        network_bID = ur.control_0_updating(G, prev_state=prev_network_bID, control_node=control_node)
        network_ID = binary_to_decimal(nodes_list, network_bID)
        G_transition_graph.add_edge(prev_network_ID, network_ID) #this is a directed graph showing you attractor landscape for how each of the states trasitions to each other
    
    return G_transition_graph

#########################################
def find_attractor_old(G_transition_graph):
    
    ''':
        Arguments:
            G_transition_graph [networkx Graph object]
                networkx graph directed graph showing how network configurations map to one another
        Return:
            attractors [dict of list of lists of ints]
                ['fixed'] = [[532][948]]
                ['cycle'] = []
        --> The output of this tells me all the cycles in the network, but it does not tell me whether they
            are accessible from the initial state that i'm interested in. need to figure this out.
    '''
    attractor_list = nx.simple_cycles(G_transition_graph) #in case of deterministic system, any cycle without considering edge direction will be directed cycle.
    attractors = {}
    attractors['fixed'] = []
    attractors['cycle'] = []

    for network_ID in attractor_list:
        # print network_ID
        if len(network_ID) == 1:
            attractors['fixed'].append(network_ID)
        else:
            attractors['cycle'].append(network_ID)

    return attractors #this outputs decID of attractor states (fixed and cyclic)

#########################################
def find_attractor(G_transition_graph):
    
    ''':
        Arguments:
            G_transition_graph [networkx Graph object]
                networkx graph directed graph showing how network configurations map to one another
        Return:
            attractors ordered by basin size (largest is first)
    '''
    attractor_list = nx.simple_cycles(G_transition_graph) #in case of deterministic system, any cycle without considering edge direction will be directed cycle.
    
    # attractlist = [i for i in copy.copy(attractor_list)]
    # print attractlist

    attractors = {}
    
    undirectedmap = nx.DiGraph.to_undirected(G_transition_graph)
    

    for u in attractor_list:
        
        if len(u) == 1:
            attractors[u[0]] = {}
            attractors[u[0]]['type'] = 'fixed'
        else:
            ## added the for statement so that the cyclic attractor  would print all attractors
            for each in range(len(u)):
                attractors[u[each]] = {}
                attractors[u[each]]['type'] = 'cycle'

    # print [(k,v) for k,v in attractors.items()]

    for v in attractors.iterkeys():
        basin = nx.node_connected_component(undirectedmap, v)
        attractors[v]['basin'] = basin
        attractors[v]['basin-size'] = len(basin)
    
    sorted_attractors = OrderedDict(sorted(attractors.items(), key=lambda kv: kv[1]['basin-size'], reverse=True))
    return sorted_attractors

#########################################
def time_series_pa(G, nodes_list, Initial_States_List, timesteps=20):
    
    ''':
        Description:
        -- compute timeseries for only a subset of nodes (Initial_States_List)
        
        Return:
        -- 1. timeSeriesData (only for primary attractor, if PA is given as the Initial_States_List)
    '''

    timeseriesdata = time_series_all(G, nodes_list, timesteps=20)
    #------------------------------------------------------------------------------
    reduced_timeseriesdata = dict()

    for node in G.nodes():
        reduced_timeseriesdata[node] = {}
        for network_ID in Initial_States_List:
            reduced_timeseriesdata[node][network_ID] = timeseriesdata[node][network_ID]


    return reduced_timeseriesdata

#########################################
def time_series_one(G, nodes_list, Initial_State, timesteps=20):
    
    '''
        Description:
        -- compute timeseries from one initial state
        
        Return:
        -- 1. timeSeriesData (only for given initial_state)
    '''

    timeseriesdata = time_series_all(G, nodes_list, timesteps=20)
    #------------------------------------------------------------------------------
    reduced_timeseriesdata = dict()

    for node in G.nodes():
        reduced_timeseriesdata[node] = {}
        reduced_timeseriesdata[node][Initial_State] = timeseriesdata[node][Initial_State]

    return reduced_timeseriesdata

#########################################
def main():
    print "time_evol module is the main code."
    edge_file = '../data/inputs/edges-init.dat' 
    # edge_file = '../data/inputs/edges_1-init.dat' 
    node_file = '../data/inputs/ant-nodes-init.dat'
    
    G = inet.read_network_from_file(edge_file, node_file)
    nodes_list = G.nodes() # print graph nodes without states. Could also do nx.nodes(G)

    # print G.nodes()
    initial_state = nx.get_node_attributes(G,'state')
    # print initial_state #OrderedDict(sorted(initial_state.items(), key=lambda t: t[0]))

    #####

    # timeseriesdata = time_series_all(G, nodes_list, 10)#, timesteps=20)

    # network_ID = 1 #this is the ID of the init state
    # biStates = decimal_to_binary(nodes_list, network_ID)
    # print 'initial state', biStates
    # for node in G.nodes():
    #     print node, timeseriesdata[node][1]

    G_transition_graph = net_state_transition(G, nodes_list) 
    
    ####################################################
    #### calculate the max number of steps to get to any given attractor in my system
    # target_from_source_dict = nx.single_source_shortest_path(G_transition_graph.to_undirected(), 43)
    # print max([len(target_from_source_dict[i]) for i in target_from_source_dict])
    # exit()
    ####################################################
    # nx.draw(G_transition_graph)
    # plt.show()

    attractors = find_attractor(G_transition_graph) # why does this change if I update G?

    print attractors.keys()
    # print nx.number_weakly_connected_components(G_transition_graph)
    # print [i for i in nx.weakly_connected_components(G_transition_graph)]
    print [len(i) for i in nx.weakly_connected_components(G_transition_graph)]

    # attractor_bIDs = [decimal_to_binary(nodes_list,k) for k in attractors.keys()]

    # for k in attractors.keys():
    #     print k, decimal_to_binary(nodes_list,k)

    ###########################
    ## this is the graph with state transitions when you add a control rule of state=1
    # control_node='pS'

    # G_control_graph = control_state_1_transition(G, nodes_list, control_node)

    # attractors = find_attractor(G_control_graph) # why does this change if I update G?

    # print attractors

    # # nx.write_graphml(G_control_graph,'G_1_control_graph_%s.xml'%control_node)

    # print attractors.keys() 
    # print nx.number_weakly_connected_components(G_control_graph)
    # print [len(i) for i in nx.weakly_connected_components(G_control_graph)]

    ###########################
    ## this is the graph with state transitions when you add a control rule of state=0
    # control_node='gC'

    # G_control_graph = control_state_0_transition(G, nodes_list, control_node)

    # attractors = find_attractor(G_control_graph) # why does this change if I update G?

    # nx.write_graphml(G_control_graph,'G_0_control_graph_%s.xml'%control_node)

    # print attractors.keys() 
    # print nx.number_weakly_connected_components(G_control_graph)
    # print [len(i) for i in nx.weakly_connected_components(G_control_graph)]
    ###########################

    

    # print decimal_to_binary(nodes_list,attractors.keys()[0])

    # print attractor_bIDs

    # print attractors.keys()[0] # this prints the decID of the primary attractor
    # print attractors.keys()
    # for i in attractors.keys()[0]:
    #     print i[0]

    # exit()

    # print [i for i in nx.attracting_components(G_transition_graph)]
    # print [i.nodes() for i in nx.attracting_component_subgraphs(G_transition_graph)]

    # print nx.number_weakly_connected_components(G_transition_graph)
    # print [i for i in nx.weakly_connected_components(G_transition_graph)]

    # print nx.number_weakly_connected_components(G_control_graph)
    # print [i for i in nx.weakly_connected_components(G_control_graph)]
    # print [len(i) for i in nx.weakly_connected_components(G_control_graph)]

    exit()

    # print nx.node_connected_component(G, 1)
    ## now it's telling me that my graph is basically disconnected?
    # print [n for n in nx.strongly_connected_components(G_transition_graph)]
    # print nx.number_strongly_connected_components(G_transition_graph)
    # print nx.is_strongly_connected(G_transition_graph)






if __name__=='__main__':
    main()
