import networkx as nx
import random
import numpy as np

def gnmw(G):
    return(len(G.nodes()),len(G.edges()),np.sum([w['weight'] for (u,v,w) in G.edges(data=True)]))

def gnmRandomWeightedGraph(nodes,edges,weight):
    G = nx.gnm_random_graph(nodes,edges)
    totalWeight = 0
    for (u,v,w) in G.edges(data=True):
        w['weight'] = random.uniform(0,10)
        totalWeight += w['weight']
    
    normFactor = weight/totalWeight
    
    for (u,v,w) in G.edges(data=True):
        w['weight'] = w['weight']*normFactor
        
    return G

def smallWorldNess(G):
    gnmwG = gnmw(G)
    RG = gnmRandomWeightedGraph(gnmwG[0],gnmwG[1],gnmwG[2])
    pLengthG = nx.average_shortest_path_length(G,weight='weight')
    pLengthRG = nx.average_shortest_path_length(RG,weight='weight')
    clusteringG = nx.average_clustering(G,weight='weight')
    clusteringRG = nx.average_clustering(RG,weight='weight')
    SW = (clusteringG/clusteringRG)/(pLengthG/pLengthRG)
    return SW
