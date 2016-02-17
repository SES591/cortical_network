#!/usr/bin/python
#input_net.py
#last update: 17 DEC 2015

__author__ = '''Hyunju Kim/Harrison Smith'''


import os
import sys
import numpy as np
import networkx as nx


#################################################################
def read_network_from_file(edge_file, node_file):
    """"""

    '''
    Get edge list, thresholds and nodes from file

    Arguments:
        edge_file [csv delimited .dat file]
            format:
            source_node [str] , target_node [str] , weight [bool]

        node_file [csv delimited .dat file]
            format:
            node [str] , state [bool]

    Outputs:
        G [networkx Graph object]


    NOTES:
        weight = 1 if edge is inductive
        weight = 0 if edge is inhibiting

        state = 1 if expressed
        state = 0 if not expressed 

        weight and state are boolean, but MUST be 0 or 1.
        DO NOT USE 'True' and 'False'. 
        They will not be evaluated correctly.
    '''

    ## read edge list with its weight from edge_file
    G = nx.read_edgelist(edge_file, create_using=nx.DiGraph(), delimiter=',', nodetype=str, data=(('weight', int),))
    ## Raise error if weight is not a 0 or 1
    for edge in nx.get_edge_attributes(G,'weight'):
        weight = nx.get_edge_attributes(G,'weight')[edge]
        if weight != 0 and weight != 1:
            raise ValueError("All weights must be boolean, represented as either 0 or 1")

    ## read node list with its state from node_file
    for line in open(node_file, 'r').readlines():
        items = [x.strip() for x in line.rstrip().split(',')]
        # print items
        if line[0] == '#' or line=='':
            continue
        G.add_node(items[0], state=int(items[1]))
    ## Raise error if state is not a 0 or 1
    for node in nx.get_node_attributes(G,'state'):
        state = nx.get_node_attributes(G,'state')[node]
        if state != 0 and state != 1:
            raise ValueError("All states must be boolean, represented as either 0 or 1")

    return G

#################################################################

def main():
    print "input_net module is the main code."
    edge_file = '../data/inputs/edges-init.dat' #'../data/example/example-net-edges.dat'
    node_file = '../data/inputs/ant-nodes-init.dat' #'../data/example/example-net-nodes.dat'

    G = read_network_from_file(edge_file, node_file)

    print "Edges:\n",G.edges(data=True) # print graph edges with weights
    print "Nodes:\n",nx.get_node_attributes(G,'state') # print nodes with states
    # nx.nodes(G) # print graph nodes without states


if __name__=='__main__':
    main()
