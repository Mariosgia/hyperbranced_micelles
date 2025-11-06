


#import hoomd
import gsd.hoomd
import math
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import random
import pydot
import os
import re 
import sys
import glob
import heapq
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec
from matplotlib.markers import MarkerStyle
from matplotlib import ticker as mticker
from matplotlib.ticker import LogLocator, MultipleLocator

def distance(x0, x1, dimensions):
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))


#########################
#### START OF MAIN  #####
#########################


#FIGURE PARAMETERS
width = 1.75
height = width/1.4
fig = plt.figure(constrained_layout=False, dpi=300)
fig.set_size_inches(width, height)


plt.rc('font', family='sans-serif', size=8)
plt.rc('text', usetex=True)
plt.rc('mathtext',fontset='dejavusans')
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=7)
plt.rc('legend', fontsize=5,fancybox=False, framealpha=0.0)
plt.rc('axes', linewidth=0.5)
plt.rc("savefig", dpi=300)
plt.rc("lines", linewidth=1., markersize=3, markeredgewidth=2.5)
plt.rcParams['axes.titlesize'] = 8
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{amssymb}\usepackage{stmaryrd}\usepackage{sfmath}'



gen=0

# Initialize the graph
G = nx.Graph()

aa=1
aa1=np.cos(np.pi/4)*aa
aa2=np.sin(np.pi/4)*aa

bb1=np.cos(np.pi/6)*aa
bb2=np.sin(np.pi/6)*aa


num_nodes = 12

pos = {}

for i in range(num_nodes):
    G.add_node(f"R{i}", color="red")
for i in range(num_nodes):
    G.add_node(f"B{i}", color="blue")


if gen==0:
    for i in range(num_nodes - 1):
        G.add_edge(f"R{i}", f"R{i+1}")
    for i in range(num_nodes - 1):
        G.add_edge(f"B{i}", f"B{i+1}")
    G.add_edge(f"R{num_nodes-1}", f"B{0}")
    pos = {}
    for i in range(num_nodes):
        pos[f"B{i}"] = (num_nodes*aa+i*aa, 0)  


if gen==1:

    for i in range(num_nodes-1):
        G.add_edge(f"R{i}", f"R{i+1}")
    for i in range(int(num_nodes/2)-1):
        G.add_edge(f"B{i}", f"B{i+1}")
    for i in range(int(num_nodes/2),num_nodes-1):
        G.add_edge(f"B{i}", f"B{i+1}")

    for i in range(int(num_nodes/2)):
        j=i+int(num_nodes/2)
        pos[f"B{i}"] = (num_nodes*aa-aa+aa1+i*aa1, aa2+i*aa2)
        pos[f"B{j}"] = (num_nodes*aa-aa+aa1+i*aa1, -aa2-i*aa2)  

    G.add_edge(f"R{num_nodes-1}", f"B{0}")
    G.add_edge(f"R{num_nodes-1}", f"B{int(num_nodes/2)}")


if gen==2:
   
    for i in range(num_nodes-1):
        G.add_edge(f"R{i}", f"R{i+1}")

    bl_len=int(num_nodes/6)
    half_len=int(num_nodes/2)

    for i in range(0,bl_len-1):
        G.add_edge(f"B{i}", f"B{i+1}")
    G.add_edge(f"B{bl_len-1}", f"B{bl_len}")
    G.add_edge(f"B{bl_len-1}", f"B{bl_len*2}")
    for i in range(bl_len,bl_len+bl_len-1):
        G.add_edge(f"B{i}", f"B{i+1}")
    for i in range(bl_len*2,bl_len*2+bl_len-1):
        G.add_edge(f"B{i}", f"B{i+1}")



    for i in range(half_len,half_len+bl_len-1):
        G.add_edge(f"B{i}", f"B{i+1}")
    G.add_edge(f"B{half_len+bl_len-1}", f"B{half_len+bl_len}")
    G.add_edge(f"B{half_len+bl_len-1}", f"B{half_len+bl_len*2}")
    for i in range(half_len+bl_len,half_len+bl_len+bl_len-1):
        G.add_edge(f"B{i}", f"B{i+1}")
    for i in range(half_len+bl_len*2,half_len+bl_len*2+bl_len-1):
        G.add_edge(f"B{i}", f"B{i+1}")


    for i in range(0,bl_len):
        pos[f"B{i}"] = (num_nodes*aa-aa+aa1+i*aa1, aa2+i*aa2)

    nod="B"+str(bl_len-1)
    for i in range(0,bl_len):
        j=bl_len+i
        pos[f"B{j}"] = (pos[nod][0]+bb1+i*bb1, pos[nod][1]+bb2+i*bb2)

    nod="B"+str(bl_len-1)
    for i in range(0,bl_len):
        j=2*bl_len+i
        pos[f"B{j}"] = (pos[nod][0]+bb1+i*bb1, pos[nod][1]-bb2-i*bb2)

    for i in range(0,half_len):
        nod="B"+str(i)
        j=half_len+i
        pos[f"B{j}"] = (pos[nod][0], -pos[nod][1])

    G.add_edge(f"R{num_nodes-1}", f"B{0}")
    G.add_edge(f"R{num_nodes-1}", f"B{int(num_nodes/2)}")

for i in range(num_nodes):
    pos[f"R{i}"] = (i*aa, 0) 


  
# Define node colors based on their color attribute
node_colors = [G.nodes[node]["color"] for node in G.nodes]

# Draw the graph


nx.draw(G, pos, node_color=node_colors, with_labels=False, node_size=3, font_weight='bold')
plt.savefig("gen_"+str(gen)+"_graph.pdf",dpi=300)


#intg=individual_graphs[polymer]

#node_colors=[]
#for node in list(intg.nodes()):
#    if np.isin(node, monomer_indices):
#        node_colors.append("blue")
#    elif np.isin(node, ind_F):
#        node_colors.append("red")
#    elif np.isin(node, ind_C):
#        node_colors.append("blue")
## Plot the graph with nodes color-coded by their connected components
#pos = nx.nx_pydot.graphviz_layout(intg)


#nx.draw(intg, pos,node_color=node_colors, with_labels=False, font_weight='bold', node_size=1)
##            plt.show()
#plt.savefig("unzoomed_mol_network.pdf",dpi=300)

#nx.draw(intg, pos,node_color=node_colors, with_labels=False, font_weight='bold', node_size=3)
##            plt.show()

## Get the current axis limits (before zooming)
#xlim = ax.get_xlim()
#ylim = ax.get_ylim()

#xlim1=xlim[0]+3500
#xlim2=xlim[1]  

#ylim1=ylim[0]+100
#ylim2=ylim[1]-1800

#plt.xlim(xlim1, xlim2)  # Adjust x-axis limits to zoom in
#plt.ylim(ylim1, ylim2)  # Adjust y-axis limits to zoom in
#plt.savefig("zoomed_mol_network.pdf",dpi=300)

#ONLY FOR VISUALISATION OF INDIVIDUAL POLYMERS 





