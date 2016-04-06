#!/usr/bin/python
#updating_rule.py
#last update: 17 DEC 2015

## this is what we'll have to edit the most

__author__ = '''Hyunju Kim/Harrison Smith'''

import os
import sys
import numpy as np
import networkx as nx
from collections import OrderedDict

import input_net as inet
    
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

            if G[u][v]['weight'] == nx.get_node_attributes(G,'state')[u]:

                nx.set_node_attributes(G, 'state', {v:1}) #change v node to expressed
            
            else: ## All conditions must be met for node to be expressed
                nx.set_node_attributes(G, 'state', {v:0}) #else, if conditions are ever not met, set node==False and break out of comparing to u nodes
                break

    curr_state = nx.get_node_attributes(G,'state') # dict of key=node, value=state of updated network G

    return curr_state

###############################################################
def main():
    print "updating_rule module is the main code."
    edge_file = '../data/inputs/edges-init.dat'
    node_file = '../data/inputs/ant-nodes-init.dat'
    
    G = inet.read_network_from_file(edge_file, node_file)

    print G.nodes()
    print G.edges()

    prev_state = nx.get_node_attributes(G,'state') # dict of key=node, value=state of network G

    print "network state @ prev step", OrderedDict(sorted(prev_state.items(), key=lambda t: t[0]))
    
    curr_state = boolean_updating(G)

    print "network state @ curr step", OrderedDict(sorted(curr_state.items(), key=lambda t: t[0]))

    print prev_state
    print curr_state

    exit()

if __name__=='__main__':
    main()
