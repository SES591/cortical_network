#!/usr/bin/python
#updating_rule.py
#last update: 17 DEC 2015

## this is what we'll have to edit the most

__author__ = '''Hyunju Kim/Harrison Smith'''

import os
import sys
import copy
import operator
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import OrderedDict
from itertools import combinations

import input_net as inet
import time_evol

    
#############################################################
def boolean_updating(G, prev_state=None):
    """
    Update according to AND and NOT boolean rules

    Arguments:
        G [networkx Graph object]

    Outputs:
        curr_state [dict]
            Dictionary of keys=nodes, values=states of graph after updating

    NOTES:
        -iterate over nodes u to connected to node v
        -this should activate node only if ALL conditions are met

        Remember: 1=True, 0=False
        
        for edges: False=inhibited, True=activated. 
        for nodes: False=not expressed, True=expressed
        
        Possibilities (node refers to u--the comparison node):
              edge=False, node=False. result == True
              edge=False, node=True. results == False
              edge=True, node=False. results == False
              edge=True, node=True. results == True
        
        Therefore: if edge==node THEN results == True. THUS v node is activated.

    """

    ## Calc prev_state from default network if not provided
    if prev_state == None:

        G_prev = G.copy() # previous state of graph
        prev_state = nx.get_node_attributes(G_prev,'state') # dict of {node:state} in previous state of graph

    ## Update each (v) nodes state, provided v isn't-inhibited/is-induced by each node u
    for v in G.nodes():

        for u in G.predecessors_iter(v):

            # print v, u

            if G[u][v]['weight'] == prev_state[u]:#nx.get_node_attributes(G_prev,'state')[u]: # changing from G to G_prev 4/13/16

                nx.set_node_attributes(G, 'state', {v:1}) #change v node to expressed
            
            else: ## All conditions must be met for node to be expressed
                nx.set_node_attributes(G, 'state', {v:0}) #else, if conditions are ever not met, set node==False and break out of comparing to u nodes
                break

    curr_state = nx.get_node_attributes(G,'state') # dict of key=node, value=state of updated network G

    return curr_state # I'm updating the graph object G based on the previous state, and also outputting G's current state.

#############################################################
def control_1_updating(G, prev_state=None, control_node='gF'):
    """
    Update according to AND and NOT boolean rules

    Arguments:
        G [networkx Graph object]

    Outputs:
        curr_state [dict]
            Dictionary of keys=nodes, values=states of graph after updating

    NOTES:
        -iterate over nodes u to connected to node v
        -this should activate node only if ALL conditions are met

        Remember: 1=True, 0=False
        
        for edges: False=inhibited, True=activated. 
        for nodes: False=not expressed, True=expressed
        
        Possibilities (node refers to u--the comparison node):
              edge=False, node=False. result == True
              edge=False, node=True. results == False
              edge=True, node=False. results == False
              edge=True, node=True. results == True
        
        Therefore: if edge==node THEN results == True. THUS v node is activated.

    """

    ## Calc prev_state from default network if not provided
    if prev_state == None:

        G_prev = G.copy() # previous state of graph
        prev_state = nx.get_node_attributes(G_prev,'state') # dict of {node:state} in previous state of graph

    ## Update each (v) nodes state, provided v isn't-inhibited/is-induced by each node u
    for v in G.nodes():

        for u in G.predecessors_iter(v):

            # print v, u

            if G[u][v]['weight'] == prev_state[u]:#nx.get_node_attributes(G_prev,'state')[u]: # changing from G to G_prev 4/13/16

                nx.set_node_attributes(G, 'state', {v:1}) #change v node to expressed
            
            else: ## All conditions must be met for node to be expressed
                if v != control_node:
                    nx.set_node_attributes(G, 'state', {v:0}) #else, if conditions are ever not met, set node==False and break out of comparing to u nodes
                else:
                    nx.set_node_attributes(G, 'state', {v:1}) ## this is the condition that the control kernal paper said would make everything go to the primary attractor. Let's see if true.
                break

    curr_state = nx.get_node_attributes(G,'state') # dict of key=node, value=state of updated network G

    return curr_state # I'm updating the graph object G based on the previous state, and also outputting G's current state.

#############################################################
def control_0_updating(G, prev_state=None, control_node='pC'):
    """
    Update according to AND and NOT boolean rules

    Arguments:
        G [networkx Graph object]

    Outputs:
        curr_state [dict]
            Dictionary of keys=nodes, values=states of graph after updating

    NOTES:
        -iterate over nodes u to connected to node v
        -this should activate node only if ALL conditions are met

        Remember: 1=True, 0=False
        
        for edges: False=inhibited, True=activated. 
        for nodes: False=not expressed, True=expressed
        
        Possibilities (node refers to u--the comparison node):
              edge=False, node=False. result == True
              edge=False, node=True. results == False
              edge=True, node=False. results == False
              edge=True, node=True. results == True
        
        Therefore: if edge==node THEN results == True. THUS v node is activated.

    """

    ## Calc prev_state from default network if not provided
    if prev_state == None:

        G_prev = G.copy() # previous state of graph
        prev_state = nx.get_node_attributes(G_prev,'state') # dict of {node:state} in previous state of graph

    ## Update each (v) nodes state, provided v isn't-inhibited/is-induced by each node u
    for v in G.nodes():

        for u in G.predecessors_iter(v):

            # print v, u

            if G[u][v]['weight'] == prev_state[u]:#nx.get_node_attributes(G_prev,'state')[u]: # changing from G to G_prev 4/13/16

                if v != control_node:
                    nx.set_node_attributes(G, 'state', {v:1}) #change v node to expressed

                else:
                    nx.set_node_attributes(G, 'state', {v:0}) #change v node to expressed
            
            else: ## All conditions must be met for node to be expressed

                nx.set_node_attributes(G, 'state', {v:0}) #else, if conditions are ever not met, set node==False and break out of comparing to u nodes
                break

    curr_state = nx.get_node_attributes(G,'state') # dict of key=node, value=state of updated network G

    return curr_state # I'm updating the graph object G based on the previous state, and also outputting G's current state.


###############################################################
### function to find differences between dictionaries

def create_state_transitions_dict(G):
    """
    output a dictionary of key=network_ID, value=list of possible next network_IDs
    each of the next possible network_IDs is exactly 1 value different than the previous network ID
    this is working correctly- trust past harrison 4/13/16
    """
    state_trans_dict = {} # key=network_ID, value=list of next network IDs

    n_nodes = len(G.nodes())
    n_init_networks = time_evol.n_states**n_nodes

    for prev_network_ID in range(0, n_init_networks):

        nodes_list = G.nodes()

        prev_network_bID = time_evol.decimal_to_binary(nodes_list, prev_network_ID) #same as 'prev_state' within boolean_updating

        curr_network_bID = boolean_updating(G, prev_state=prev_network_bID)
        
        # get list of nodes that change when applying update rule
        altered_nodes = []
        for node in curr_network_bID:

            if prev_network_bID[node] != curr_network_bID[node]:
                altered_nodes.append(node)

        # print prev_network_ID, altered_nodes

        # look through altered nodes, and make list of all networks with nodes being altered 1 at a time
        state_trans_dict[prev_network_ID] = []
        for node in altered_nodes:
            
            possible_bID = copy.deepcopy(prev_network_bID)
            # print possible_bID #why is this not changing???

            not_flipped = True
            
            # print possible_bID[node]
            # if node was updated, flip it
            if prev_network_bID[node] == 0:
                possible_bID[node] = 1
                not_flipped = False
                # print "0, flipping"
            
            # print possible_bID[node]
             # if node was updated, flip it
            elif (prev_network_bID[node] == 1) and (not_flipped == True):
                possible_bID[node] = 0
                # print "1, flipping"

            # print possible_bID[node]
            else: 
                raise ValueError("Node %s of prev_state must be 0 or 1. Something bad happened."%node)

            # print possible_bID
            # nodes_list = [k for k in possible_bID]

            network_ID = time_evol.binary_to_decimal(nodes_list,possible_bID) #get network_ID of this possible next state

            state_trans_dict[prev_network_ID].append(network_ID)

        # if no nodes are altered, you're in a fixed stable state-- make state_trans_dict[stable_ID]=stable_ID
        if len(altered_nodes) == 0:
            state_trans_dict[prev_network_ID].append(prev_network_ID)
            # print 'prev network ID:',prev_network_ID

    return state_trans_dict
###############################################################
def create_graph_from_transitions(state_trans_dict):
    """ create graph of all possible state transitions for a stochastic system"""

    all_node_paths = []
    for k in state_trans_dict:
        node_path = [k] + state_trans_dict[k]
        all_node_paths.append(node_path)

    G = nx.DiGraph()

    for path in all_node_paths:
        G.add_path(path)

    return G
###############################################################
def create_transition_matrix(state_trans_dict):
    """ calculates transition matrix when given the state_transition_dict"""
    transition_matrix = np.zeros([len(state_trans_dict),len(state_trans_dict)])
    for nID in state_trans_dict:
        for nextID in state_trans_dict[nID]:
            # transition_matrix[nID,nextID] = 1.0/len(state_trans_dict[nID]) # this is equally distributed through rows--> want through columns
            transition_matrix[nextID,nID] = 1.0/len(state_trans_dict[nID])

    for col in range(1024):
        assert sum(transition_matrix[:,col]) > .99 # all cols should sum to 1 (this accts for machine error)

    return transition_matrix

###############################################################
def make_attractor_list(attractors):
    """make flat list of attractors (includes fixed and cyclic) from the time_evol.find_attractor_old 
    function """
    attractor_list = []
    for attracttype in attractors:
        for attractor in attractors[attracttype]:
            for node in attractor:
                attractor_list.append(node)

    return attractor_list

###############################################################
def find_steady_state(T,attractor_list,init_dID,desired_dID):
    """find steady state of transition matrix in stochastic system using init_network ID and known steady states from deterministic state transition dict"""
    from time_evol import n_states
    # T is the transition matrix
    n_nodes = 10
    s = np.zeros([n_states**n_nodes,1]) # vector of probability of states

    s[init_dID] = 1.0
    cutoff = 0.9999
    total_attractor_probabilities = sum([s[a-1] for a in attractor_list])

    n_steps = 0
    while  total_attractor_probabilities < cutoff:
        s = np.dot(T,s)
        assert sum(s) > 0.99 # every s should sum to 1, this accts for machine error
        total_attractor_probabilities = sum([s[a] for a in attractor_list])
        n_steps +=1
        print n_steps, total_attractor_probabilities

    print n_steps

    attractor_dict = {}
    for a in attractor_list:
        attractor_dict[a] = s[a]
    
    return attractor_dict, s[desired_dID]

###############################################################
def full_steady_state_calc(G,nodes_list,init_dID,desired_dID):
    """ run through full steady state calculation from stochastic transition matrix
    to find the probability of being in each attractor during steady state """

    state_trans_dict = create_state_transitions_dict(G)
    deterministic_transition_graph = time_evol.net_state_transition(G, nodes_list) # deterministic transition graph
    attractors = time_evol.find_attractor_old(deterministic_transition_graph) # this must only print out cycles for a subset of the graph...
    attractor_list = make_attractor_list(attractors)
    T = create_transition_matrix(state_trans_dict)
    attractor_dict, prob_desired_dID = find_steady_state(T,attractor_list,init_dID,desired_dID)
    return attractor_dict, prob_desired_dID

###############################################################
def find_len_1_control_kernals(deterministic_transition_graph, nodes_list, attractor_ID):
    """uses the deterministic_transition_graph and specified attractor to find 
    all single control nodes if they exist"""

    subgraphs = [g for g in nx.weakly_connected_components(deterministic_transition_graph)]
    # index_of_largest_subgraph = max(enumerate(all_subgraph_sets), key = lambda tup: len(tup[1]))[0]

    all_subgraph_sets = []
    all_but_attractor_subgraph = []
    for sg in subgraphs:
        if attractor_ID not in sg:
            all_state_sets = []
            for nID in sg:
                all_state_sets.append(set(time_evol.decimal_to_binary(nodes_list,nID).items()))
            all_subgraph_sets.append(all_state_sets)
            all_but_attractor_subgraph.append(all_state_sets)
        else:
            all_state_sets = []
            for nID in sg:
                all_state_sets.append(set(time_evol.decimal_to_binary(nodes_list,nID).items()))
            all_subgraph_sets.append(all_state_sets)

    state_0_list = [0 for i in range(len(nodes_list))]
    state_1_list = [1 for i in range(len(nodes_list))]

    possible_states = zip(nodes_list,state_0_list) + zip(nodes_list,state_1_list) # list of (node,state)

    possible_states_pared = copy.deepcopy(possible_states)

    for sg in all_but_attractor_subgraph:

        for state_set in sg:

            for state in possible_states:

                if (state in state_set) and (state in possible_states_pared):

                    possible_states_pared.remove(state)

    return possible_states_pared

###############################################################
def find_len_2_control_kernals(deterministic_transition_graph, nodes_list, attractor_ID):
    """uses the deterministic_transition_graph and specified attractor to find 
    all pairs of control nodes if they exist

    note: these aren't "strict" control kernels because they specify the states needed to be in the
    main attractor. controlling them doesn't necessarily change what attractor you'll be in.
    """


    subgraphs = [g for g in nx.weakly_connected_components(deterministic_transition_graph)]
    # index_of_largest_subgraph = max(enumerate(all_subgraph_sets), key = lambda tup: len(tup[1]))[0]

    all_subgraph_sets = []
    all_but_attractor_subgraph = []
    for sg in subgraphs:
        if attractor_ID not in sg:
            all_state_sets = []
            for nID in sg:
                all_state_sets.append(set(time_evol.decimal_to_binary(nodes_list,nID).items()))
            all_subgraph_sets.append(all_state_sets)
            all_but_attractor_subgraph.append(all_state_sets)
        else:
            all_state_sets = []
            for nID in sg:
                all_state_sets.append(set(time_evol.decimal_to_binary(nodes_list,nID).items()))
            all_subgraph_sets.append(all_state_sets)

    state_0_list = [0 for i in range(len(nodes_list))]
    state_1_list = [1 for i in range(len(nodes_list))]

    possible_states = zip(nodes_list,state_0_list) + zip(nodes_list,state_1_list) # list of (node,state)

    possible_pairs = [combo for combo in combinations(possible_states, 2)]

    possible_pairs_pared = copy.deepcopy(possible_pairs)

    # remove pairs where both keys are the same (eg. gF and gF)
    for pair in possible_pairs:

        if pair[0][0] == pair[1][0] and (pair in possible_pairs_pared):

            possible_pairs_pared.remove(pair)

    # remove pairs when any of the networks in the non-attractor subgraph contain that pair
    # (ie. that pair can not possibly be a control kernel because it is present in the wrong attractor)
    for sg in all_but_attractor_subgraph:

        for state_set in sg:

            for pair in possible_pairs:

                if (pair[0] in state_set) and (pair[1] in state_set) and (pair in possible_pairs_pared):

                    possible_pairs_pared.remove(pair)

    return possible_pairs_pared # return a list ((node,state),(node,state)) pairs that are control kernels

###############################################################
def find_len_2_control_kernals_fo_real(deterministic_transition_graph, nodes_list, attractor_ID):
    """uses the deterministic_transition_graph and specified attractor to find 
    all pairs of control nodes if they exist

    note: THESE ARE STRICT CONTROL KERNALS--but they don't appear to actually exist
    """


    subgraphs = [g for g in nx.weakly_connected_components(deterministic_transition_graph)]
    # index_of_largest_subgraph = max(enumerate(all_subgraph_sets), key = lambda tup: len(tup[1]))[0]

    # combos_of_all_subgraph_sets = []
    all_subgraph_sets = []
    all_but_attractor_subgraph = []
    for sg in subgraphs:
        if attractor_ID not in sg:
            all_state_sets = []
            for nID in sg:
                all_state_sets.append(set(time_evol.decimal_to_binary(nodes_list,nID).items()))
            all_subgraph_sets.append(all_state_sets)
            all_but_attractor_subgraph.append(all_state_sets)
        else:
            all_state_sets = []
            for nID in sg:
                all_state_sets.append(set(time_evol.decimal_to_binary(nodes_list,nID).items()))
            all_subgraph_sets.append(all_state_sets)

    # print [len(i) for i in all_subgraph_sets]
    # print len(all_subgraph_sets[0])

    state_0_list = [0 for i in range(len(nodes_list))]
    state_1_list = [1 for i in range(len(nodes_list))]

    possible_states = zip(nodes_list,state_0_list) + zip(nodes_list,state_1_list) # list of (node,state)

    possible_pairs = [combo for combo in combinations(possible_states, 2)]

    possible_pairs_pared = copy.deepcopy(possible_pairs)

    # remove pairs where both keys are the same (eg. gF and gF)
    for pair in possible_pairs:

        if pair[0][0] == pair[1][0] and (pair in possible_pairs_pared):

            possible_pairs_pared.remove(pair)


    list_of_all_but_1s = [i for i in combinations(all_subgraph_sets, 2)]

    list_of_control_pairs = []

    # remove pairs when any of the networks in the non-attractor subgraph contain that pair
    # (ie. that pair can not possibly be a control kernel because it is present in the wrong attractor)
    for all_but_1_graph in list_of_all_but_1s:

        temp_possible_pairs_pared = copy.deepcopy(possible_pairs_pared)

        for sg in all_but_1_graph:

            for state_set in sg:

                for pair in possible_pairs:

                    if (pair[0] in state_set) and (pair[1] in state_set) and (pair in temp_possible_pairs_pared):

                        temp_possible_pairs_pared.remove(pair)

        list_of_control_pairs.append(temp_possible_pairs_pared)

    return list_of_control_pairs # return a list ((node,state),(node,state)) pairs that are control kernels



###############################################################
def test_create_state_transition_dict(state_trans_dict):
    """ a function to test if the create_state_trans_dict function is working as intended"""

    # print state_trans_dict
    for nID in state_trans_dict:
        # print nID
        bID = time_evol.decimal_to_binary(nodes_list,nID)
        # print "."*40
        for possible_nID in state_trans_dict[nID]:
            bID_next = time_evol.decimal_to_binary(nodes_list,possible_nID)

            if abs(sum(bID.values())-sum(bID_next.values())) != 1:
                print " more than one thing changed. wtf."

###############################################################
def main():
    print "updating_rule module is the main code."

    edge_file = '../data/inputs/edges-init.dat'
    
    #############################
    ## Read in anterior nodes
    #############################
    node_file = '../data/inputs/ant-nodes-init.dat'
    desired_node_file = '../data/inputs/ant-nodes-final.dat'
    
    G = inet.read_network_from_file(edge_file, node_file)

    nodes_list = G.nodes()
    node_dict = nx.get_node_attributes(G,'state')
    init_dID = time_evol.binary_to_decimal(nodes_list, node_dict)

    #############################
    ## Read in posterior nodes
    #############################
    post_node_file = '../data/inputs/post-nodes-init.dat'
    post_desired_node_file = '../data/inputs/post-nodes-final.dat'

    G_post = inet.read_network_from_file(edge_file, post_node_file)

    post_nodes_list = G_post.nodes()
    post_node_dict = nx.get_node_attributes(G_post,'state')
    post_init_dID = time_evol.binary_to_decimal(post_nodes_list, post_node_dict)
    
    #############################
    # print G.nodes()
    # print G.edges()

    # prev_state = nx.get_node_attributes(G,'state') # dict of key=node, value=state of network G

    # print "network state @ prev step", OrderedDict(sorted(prev_state.items(), key=lambda t: t[0]))
    
    # curr_state = boolean_updating(G)

    # print "network state @ curr step", OrderedDict(sorted(curr_state.items(), key=lambda t: t[0]))

    # print prev_state
    # print curr_state

    ###### calc steady state probabilities #####################################################

    # init_dID,desired_dID = 0, 43
    init_dID,desired_dID = 20, 980
    attractor_dict, prob_desired_dID = full_steady_state_calc(G,nodes_list,init_dID,desired_dID)
    print attractor_dict
    print prob_desired_dID

    exit()

    #############################################################################################


    ## Create deterministic transition graph
    state_trans_dict = create_state_transitions_dict(G)
    G_state_trans = create_graph_from_transitions(state_trans_dict) # stochastic transition graph

    ## Calculate in and out degrees of transition graph
    out_deg_dict = G_state_trans.out_degree(G_state_trans.nodes())
    sorted_out_deg_dict = sorted(out_deg_dict.items(), key=operator.itemgetter(1))

    in_deg_dict = G_state_trans.in_degree(G_state_trans.nodes())
    sorted_in_deg_dict = sorted(in_deg_dict.items(), key=operator.itemgetter(1))

    largest = max(nx.strongly_connected_components(G_state_trans), key=len)
    # print len(largest) # this shows that graph of stochastic state transitions is strongly connected
    # attractors = time_evol.find_attractor_old(G_state_trans)

    create_transition_matrix(state_trans_dict)

    deterministic_transition_graph = time_evol.net_state_transition(G, nodes_list) # deterministic transition graph
    attractors = time_evol.find_attractor_old(deterministic_transition_graph) # this must only print out cycles for a subset of the graph...

    #############################
    ## Find control kernals
    #############################

    attractor_ID = 980

    # print find_len_1_control_kernals(deterministic_transition_graph, nodes_list, attractor_ID)

    # print find_len_2_control_kernals(deterministic_transition_graph, nodes_list, attractor_ID)

    print find_len_2_control_kernals_fo_real(deterministic_transition_graph, nodes_list, attractor_ID)

    exit()

    #############################
    ## Desired node states
    ##      Anterior
    G_desired = inet.read_network_from_file(edge_file, desired_node_file)
    desired_nodes_list = G_desired.nodes()
    desired_node_dict = nx.get_node_attributes(G_desired,'state')
    desired_dID = time_evol.binary_to_decimal(desired_nodes_list,desired_node_dict)
    print "desired_node_dict:", desired_node_dict
    print "desired dID:", desired_dID

    ##      Posterior
    G_desired_post = inet.read_network_from_file(edge_file, post_desired_node_file)
    post_desired_nodes_list = G_desired_post.nodes()
    post_desired_node_dict = nx.get_node_attributes(G_desired_post,'state')
    post_desired_dID = time_evol.binary_to_decimal(post_desired_nodes_list,post_desired_node_dict)
    print "post_desired_node_dict:", post_desired_node_dict
    print "desired dID:", post_desired_dID

    ## Check if path between desired init state and final state
    print nx.has_path(deterministic_transition_graph,init_dID,desired_dID)

    ## Print binary version of attractors
    print "attractors:"
    print attractors

    ## Network 980 is the steady state network, let's see if this is correct based on rules
    # steady_state_network = 
    print "end_ant:", desired_node_dict
    print "beg_ant:", node_dict

    print "end_ant_dID:", desired_dID
    print "beg_ant_dID:", init_dID

    print "end_post:", post_desired_node_dict
    print "beg_post:", post_node_dict

    print "end_post_dID:", post_desired_dID
    print "beg_post_dID:", post_init_dID    




    # print attractors

    # print G_state_trans.edges()
    # for line in nx.generate_edgelist(G_state_trans, data=False):
    #     print line
    # nx.draw(G_state_trans)
    # nx.write_graphml(G_state_trans,'g_state_trans.xml')
    # plt.show()

    # for graph in list(nx.weakly_connected_component_subgraphs(deterministic_transition_graph)):
    #     if 43 in graph.nodes():
    #         print [time_evol.decimal_to_binary(nodes_list, dID) for dID in graph.nodes()]

    # print "  0:", time_evol.decimal_to_binary(nodes_list, 0)
    # print "683:", time_evol.decimal_to_binary(nodes_list, 683)


    # nx.draw(deterministic_transition_graph)
    # nx.write_graphml(deterministic_transition_graph,'deterministic_transition_graph.xml')
    # plt.show()





    exit()

if __name__=='__main__':
    main()
