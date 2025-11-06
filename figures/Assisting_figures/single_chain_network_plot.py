########################################################
########################################################
###     SCRIPT USED FOR PLOTTING THE NETWORK     #######
###     OF A SINGLE POLYMER MOLECULE             #######
###        CREATED BY MARIOS GIANNAKOU           #######
########################################################
########################################################



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
width = 3.5
#height = width/1.4
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



gs1=GridSpec(1,1,  left=0.15, bottom=0.14, right=.95, top=.96 )

ax=plt.subplot(gs1[0,0])






pattern = '../../results/micelles/LHBC/dw_1/MD_runs/batches/batch_1/final*gsd'
nam_groups = np.array(glob.glob(pattern))
batch_num=np.array([float(re.findall(r"[-+]?(?:\d*\.*\d+)",i)[0]) for i in nam_groups])


sort_ind=np.argsort(batch_num)
nam_groups=nam_groups[sort_ind]












##Parameters for MD output
cut_distance=0.5
beads_per_mol=8
number_of_molecules=1 #Polymer molecules expected in each batch



#fol1='batch_'+str(batch_start)+'_'+str(batch_end)




mol_counter=0 #Used to count polymer molecules
size_polymer=[]
size_abmonomer=[]
num_mol=len(nam_groups) #Number of molecules in topology file




for batch in nam_groups:

   
    config_file=batch
    print("Batch number "+ str(config_file))
  


    ########################################################
    ## FINDING BONDS MADE AND EXPORTING THEM INTO A GRAPH ##
    ########################################################


    #Importing structure 
    #['F','C','A','B','D','E']
    traj = gsd.hoomd.open(name=config_file, mode='r')

    end_frame = traj.__getitem__(len(traj)-1)
    positions=end_frame.particles.position

    total_N=end_frame.particles.N #Number of beads
    num_F_C=np.sum([end_frame.particles.typeid==0]*1)+np.sum([end_frame.particles.typeid==1]*1) #Number of linear and core beads

    num_F=int(np.sum([end_frame.particles.typeid==0]*1)/number_of_molecules*1.0) #Will be used to set an artificial node later
 
    #We give each monomer an index number that will be a node in the graph and we start enumarating from num_F_C+1 and on
    indices=np.arange(0,total_N)
    mol_N=int(np.sum([end_frame.particles.typeid==2]*1)) #Number of monomers
    monomer_indices=indices[end_frame.particles.typeid==2]
    monomer_indices=np.arange(0,len(monomer_indices))+num_F_C

    
    pos_C=positions[end_frame.particles.typeid==1]
    ind_C=indices[end_frame.particles.typeid==1]
    ind_F=indices[end_frame.particles.typeid==0]
    mol_C=ind_C
    pos_A=positions[end_frame.particles.typeid==2]
    ind_A=monomer_indices
    mol_A=monomer_indices
    pos_B=positions[end_frame.particles.typeid==3]
    ind_B=indices[end_frame.particles.typeid==3]
    mol_B=np.repeat(ind_A, 2)

    size_abmonomer.append(len(pos_A))

    ind_F=indices[end_frame.particles.typeid==0]
    box_dimensions=end_frame.configuration.box[0:3]



    #Molecular indices (1 to number of beads in system) where [[indexA,indexB,distance],...]
    bonded_AB=[] 

    for i in range(0,len(pos_A)):
        dis=distance(pos_A[i],pos_B,box_dimensions)
        dis_argmin=np.argmin(dis)
        dis_min=dis[dis_argmin]

        if dis_min < cut_distance :
            bonded_AB.append([i,dis_argmin,dis_min])
        else:
            bonded_AB.append([i,-1,-1])

    bonded_AC=[]
    for i in range(0,len(pos_A)):
        dis=distance(pos_A[i],pos_C,box_dimensions)
        dis_argmin=np.argmin(dis)
        dis_min=dis[dis_argmin]

        if dis_min < cut_distance :
            bonded_AC.append([i,dis_argmin,dis_min])
        else:
            bonded_AC.append([i,-1,-1])

    bonded_FF=end_frame.bonds.group[end_frame.bonds.typeid==0]
    bonded_FC=end_frame.bonds.group[end_frame.bonds.typeid==1]
    bonded_CC=end_frame.bonds.group[end_frame.bonds.typeid==2]


    bonds=[]
#    colors=[]
    for i in range(0,len(bonded_AB)):
        if bonded_AB[i][1]!=-1: 
            bonds.append([mol_A[bonded_AB[i][0]], mol_B[bonded_AB[i][1]], "A-B" ])
#            colors.append('blue')
    for i in range(0,len(bonded_AC)):
        if bonded_AC[i][1]!=-1: 
            bonds.append([mol_A[bonded_AC[i][0]], mol_C[bonded_AC[i][1]], "A-C" ])
#            colors.append('blue')
    for i in range(0,len(bonded_FF)):
        bonds.append([bonded_FF[i][0],bonded_FF[i][1], "F-F" ])
#        colors.append('red')
    for i in range(0,len(bonded_FC)):
        bonds.append([bonded_FC[i][0],bonded_FC[i][1], "F-C" ])
#        colors.append('red')
    for i in range(0,len(bonded_CC)):
        bonds.append([bonded_CC[i][0],bonded_CC[i][1], "C-C" ])
#        colors.append('yellow')





    ##########################################
    ## PROCESSING BONDS FOR OUTPUT TO SCFT  ##
    ##########################################



###### Will need to determine the end of the block using num_F

    G = nx.Graph()
    for row in bonds:
        G.add_edge(row[0], row[1])

    components = list(nx.connected_components(G))


    gen0_node_start=[min(item) for item in components] #Since we begin with the linear parts in our simulations the minimum would always be the first bead in the network


    individual_graphs = [G.subgraph(component).copy() for component in components] #Seperates each graph


    #Now we seperate for each component
    if len(components) !=  number_of_molecules:
        print("Polymers found differ from number of polymers expected")
    else:

        for polymer in range(0,len(individual_graphs)):

            intg=individual_graphs[polymer]
            

            node_start=gen0_node_start[polymer]
            total_monomers=len(individual_graphs[polymer])
            size_polymer.append(total_monomers)
            
            gen0_node_end=node_start+num_F-1 #This is also determined by prior knowledge of our system i.e. monodispersity of linear block
            #print(gen0_node_end)
        #ONLY FOR VISUALISATION OF INDIVIDUAL POLYMERS

            node_colors=[]
            for node in list(intg.nodes()):
                if np.isin(node, monomer_indices):
                    node_colors.append("blue")
                elif np.isin(node, ind_F):
                    node_colors.append("red")
                elif np.isin(node, ind_C):
                    node_colors.append("blue")
            # Plot the graph with nodes color-coded by their connected components
            pos = nx.nx_pydot.graphviz_layout(intg)


            nx.draw(intg, pos,node_color=node_colors, with_labels=False, font_weight='bold', node_size=1)
#            plt.show()
            plt.savefig("unzoomed_mol_network.pdf",dpi=300)

            nx.draw(intg, pos,node_color=node_colors, with_labels=False, font_weight='bold', node_size=3)
#            plt.show()

            # Get the current axis limits (before zooming)
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            
            xlim1=xlim[0]+3500
            xlim2=xlim[1]  

            ylim1=ylim[0]+100
            ylim2=ylim[1]-1800
  
            plt.xlim(xlim1, xlim2)  # Adjust x-axis limits to zoom in
            plt.ylim(ylim1, ylim2)  # Adjust y-axis limits to zoom in
            plt.savefig("zoomed_mol_network.pdf",dpi=300)

        #ONLY FOR VISUALISATION OF INDIVIDUAL POLYMERS 





