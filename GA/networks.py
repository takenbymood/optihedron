import networkx as nx
import pydotplus
from networkx.drawing.nx_pydot import graphviz_layout
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import random

def createMegaStar(l,r,k):
	G = nx.DiGraph()
	G.add_node(1)
	cn = 2
	npl = r+k
	n = l*npl
	for i in range(l):
		for j in range(0,r):
			G.add_node(cn)
			cn += 1
		G.add_node(cn)
		an = cn
		for j in range(0,r):
			G.add_edge(an-j-1,an)
		cn+=1
		kn = cn
		for j in range(0,k):
			G.add_node(cn)
			G.add_edge(an,cn)
			cn += 1
		for j in range(0,k):
			G.add_edge(kn+j,1)
			for m in range(k):
				if m != j:
					G.add_edge(kn+j,kn+m)

	return G

def createBiStar(n):
	G = nx.DiGraph()
	G.add_node(0)
	G.add_edge(0,0)
	for i in range(1,n+1):
		G.add_node(i)
		G.add_edge(i,0)
		G.add_edge(0,i)
	return G

def createStar(n):
	G = nx.DiGraph()
	G.add_node(0)
	G.add_edge(0,0)
	for i in range(1,n+1):
		G.add_node(i)
		G.add_edge(i,0)
	return G

def createIslands(n):
	G = nx.DiGraph()
	for i in range(n):
		G.add_node(i)
	for i in range(n):
		for j in range(n):
			if i!=j:
				G.add_edge(i,j)
	return G

def createSinglets(n):
	G = nx.DiGraph()
	for i in range(n):
		G.add_node(i)
	return G

def drawGraph(G):
	pos=nx.spring_layout(G)
	nx.draw(G)
	plt.savefig("node_colormap.png") # save as png
	plt.show() # display