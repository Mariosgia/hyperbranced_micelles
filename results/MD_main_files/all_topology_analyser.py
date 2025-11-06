########################################################
########################################################
###     SCRIPT FOR HYPERBRANCHED POLYMERS        #######
###     ANALYSIS AND PORTING TO SCFT CODE        #######
###        CREATED BY MARIOS GIANNAKOU           #######
########################################################
########################################################

########################################################
### WE BREAK THE LHBC INTO A LINEAR BLOCK (GENERATION=-1)
### AND THE CORE AND AB_2 MONOMERS ARE TREATED THE SAME
### SO IN THE END WE ONLY HAVE TWO TYPES OF MONOMERS
########################################################

########################################################
### INDIVIDUAL BATCH FILES "batch_*/final_configuration.gsd" 
### WHERE * IS AN INTEGER ARE READ AND THE TOPOLOGY IS WRITTEN
### AND EXPORTED INTO "topology.dat" TO BE READ BY THE SCFT CODE
### NOTE THAT ONE CAN HAVE MULTIPLE POLYMERS IN THE GSD FILE 
### BUT THEY SHOULD HAVE THE SAME NUMBER OF MONOMERS IN THE LINEAR BLOCK
########################################################

########################################################
###  THE CODE WORKS BY FINDING FIRST ALL THE BONDS MADE 
###  AND FINIDING OUT THE INDIVIDUAL NETWORKS FORMED
###  THEN IT GOES OVER EACH NETWORK, FINDS THE BRANCHING
###  AND TERMINAL POINTS AND MAKES SIMPLIFIED SUBGRAPHS 
###  OUT OF IT. THEN GENERATIONS AND OTHER RELEVANT
###  INFORMATION ARE EXPORTED FOR USAGE IN SCFT
########################################################

########################################################
### IN THE TOPOLOGY FILE THE MEANING OF EACH COLUMN IS:
### 0. ID FOR EACH SEGMENT (COULD BE ANYTHING EXCEPT 0) 1. LENGTH OF SEGMENT 
### 2. TYPE OF THE SEGMENT -1=STEM 0=TERMINAL 1=INTERNAL
### 4. SOURCE BRANCH NODE 5. TARGET BRACH NODE 6.GENERATION OF SEGMENT
### THE LAST TWO COLUMNS ARE FOR ACCELERATION OF SCFT
### 7. SIMILAR TO ID 0= SAFE CHOICE. IF SET TO SOME ID THEN 
### IT WILL USE IT'S PROPAGATORS TO SAVE UP TIME. 
### 8. 2= USES ONLY QQ FROM SIMILAR ID 3= USES ONLY QQDAG FROM SIMILAR ID
###    1= USES BOTH QQ AND QQDAG FROM SIMILAR ID 0=USES NONE
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

def distance(x0, x1, dimensions):
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))


#########################
#### START OF MAIN  #####
#########################






#Splitting batches into 1 group
pattern = 'MD_runs/batches/batch*/final*gsd'
nam_groups = np.array(glob.glob(pattern))
batch_num=np.array([float(re.findall(r"[-+]?(?:\d*\.*\d+)",i)[0]) for i in nam_groups])


sort_ind=np.argsort(batch_num)
nam_groups=nam_groups[sort_ind]












##Parameters for MD output
cut_distance=0.5
beads_per_mol=8
number_of_molecules=1 #Polymer molecules expected in each batch



#fol1='batch_'+str(batch_start)+'_'+str(batch_end)




fol1='all_batch'
fol2=fol1+'/index_gens_info'

os.makedirs(fol1,exist_ok=True)
os.makedirs(fol1+'/hyperbranched_data',exist_ok=True)
os.makedirs(fol2,exist_ok=True)

mol_counter=0 #Used to count polymer molecules
size_polymer=[]
size_abmonomer=[]
num_mol=len(nam_groups) #Number of molecules in topology file


#Check whether data.dat exists and how many molecules exist already
file_name = fol1+'/hyperbranched_data/topology_data.dat'

if os.path.exists(file_name):
    print("File already exists")
else:
    with open(file_name, 'w') as file:
        file.write(f'# 0.Id 1.Length 2.Type 3.Source node 4.Target node 5.Generation 6. Monomer_type 7.Similar to id 8.Accelerate\n') 
        file.write(f'#Mol-* 1.No. segments \n')  
        file.write(f'Number_Molecules, '+ str(num_mol)+'\n')    

########################################################
## SETTING UP CASES THAT SHOULD BE STUDIED ##
########################################################



for batch in nam_groups:

   
    config_file=batch
    print("Batch number "+ str(config_file))
   # with open('test.dat', 'a') as filetest:
   #     filetest.write("Batch number : " +str(config_file)+"\n")
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
   

    #print("Number of F beads is : " +str(num_F))
    #print(num_F_C)
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

#            node_colors=[]
#            for node in list(intg.nodes()):
#                if np.isin(node, monomer_indices):
#                    node_colors.append("blue")
#                elif np.isin(node, ind_F):
#                    node_colors.append("red")
#                elif np.isin(node, ind_C):
#                    node_colors.append("yellow")
#            # Plot the graph with nodes color-coded by their connected components
#            pos = nx.nx_pydot.graphviz_layout(intg)


#            nx.draw(intg, pos,node_color=node_colors, with_labels=True, font_weight='bold', node_size=100)
#            plt.show()

        #ONLY FOR VISUALISATION OF INDIVIDUAL POLYMERS 



        



            # Calculate the degrees of all nodes
            degrees = intg.degree()
            dic_degrees=dict(degrees)

            

            #### Find nodes with a degree of 2 or more
            nodes_with_1 = [x[0] for x in degrees if x[1]==1]
            nodes_with_2 = [x[0] for x in degrees if x[1]==2]
            nodes_with_3 = [x[0] for x in degrees if x[1]==3]

#            print(nodes_with_1)
#            print(nodes_with_2)
#            print(nodes_with_3)

  

            #Here we artifficially add the end node of the linear block. Comment out if you want the real branch points [gen0_node_end] 
            ############            
            nodes_with_3+=[gen0_node_end]
            dic_degrees[gen0_node_end] = 3
            #############
            nodes_with_both=nodes_with_1+nodes_with_3
            n1=len(nodes_with_1)
            n3=len(nodes_with_3)
            ntot=n1+n3

   
            ##Process for constracting sub graph with degree 1 or 3 nodes only
            path=nx.shortest_path(G, source=nodes_with_3[0], target=nodes_with_1[-1])
            segments={}
            toskip=[]

            for i in range(0,n3):

                src=nodes_with_3[i]

                for j in range(0,ntot):
                    
                    tar=nodes_with_both[j]

                    if src==tar: #Skip the case when source and target is the same node
                        continue

                    if [src,tar] in toskip:
                        continue

                    if (tar in nodes_with_1) : #Find type of segment, internal or terminal
                        typ=0 #Terminal
                    else:
                        typ=1 #Internal


                    path=nx.shortest_path(G, source=src, target=tar)
                    inbetween_path=path[1:-1]
                    degrees_inbetween=[ dic_degrees[i]  for i in inbetween_path ]
                    logical=(1 and 3) in degrees_inbetween
                    if logical==False: #If there are no degree 1 or 3 nodes inbetween then save the connection between src and tar is direct

                        len_path=len(path)-1 #Just by our definition
                        segments[tuple([src,tar])]=[len_path,typ]
                        if j>n1: #i.e. it a degree 3 node
                            toskip.append([tar,src])
                    
            n_segm=len(segments)

#            print(segments)




            #Choosing the linear segment of the polymer as the generation 0 segment
            terminal_segments= [key for key, value in segments.items() if value[1] == 0]
            g0=[tup for tup in terminal_segments if node_start in tup][0] 
            segments[g0][0]+=1 #By our definition 
            segments[g0][1]=-1 #Marking this typ so this is the stem

            # Initialize a dictionary to store the generation of nodes and edges
            generations = {}
            #Now finding generation segments based on generation_0_edge
            generation_0_segment = g0
            inv_g0=tuple(reversed(g0))  

            generations[generation_0_segment] = 0
            subedges=[]
            for i in segments:
                subedges.append(i)

            # Create an empty undirected graph
            subG = nx.Graph()


            # Add edges to the graph
            subG.add_edges_from(subedges)
#                pos = nx.nx_pydot.graphviz_layout(subG)
#                nx.draw(subG,pos, with_labels=True, font_weight='bold', node_size=100)
#                plt.show()

            queue = [generation_0_segment]
            visited_edges = set([generation_0_segment])  # To keep track of visited edges

            while queue:
                current_edge = queue[0]
                queue.pop(0)
                neighbour_edges=[]
                
                for neighbor in subG.neighbors(current_edge[1]):
                    new_edge = (current_edge[1], neighbor)

                    if new_edge not in visited_edges:
                        
                        seg_exist1=segments.get(new_edge)
                        seg_exist2=segments.get(tuple(reversed(new_edge)))
                        
                        #Recordin direction with ascending generation order
                        if seg_exist1!= None:
                            segments[new_edge][:]+=list(new_edge) 
                                
                        if seg_exist2!= None:
                            segments[tuple(reversed(new_edge))][:]+=list(new_edge) 



                        ad=1
                        if inv_g0 in [new_edge]:
                            ad=0

                        generations[new_edge] = generations[current_edge] + ad


                        queue.append(new_edge)
                        visited_edges.add(new_edge)
                        visited_edges.add((new_edge[1], new_edge[0]))  # Mark the reverse edge as visited

            del generations[inv_g0] #Overcounted edge

            # Print the generation of each edge
            for edge, generation in generations.items():

                seg_exist1=segments.get(edge)
                seg_exist2=segments.get(tuple(reversed(edge)))  
              
                if seg_exist1!= None:
                    segments[edge].append(generation)         
                if seg_exist2!= None:
                    segments[tuple(reversed(edge))].append(generation)  

        #Assigning ids to each segment
        i=1
        for key, value in segments.items():
            value.insert(0,i)
            i+=1

        #Assigning monomer type to each segment
        for key, value in segments.items():
            if value[5] == 0:
                value.append(1)
            else:
                value.append(2)

        #Assigning similar to ids for acceleration of SCFT code. Here we simply want to calculate qqdag only once for the longest terminal group
        #First finding the most suitable segment
        len0=0
        for key, value in segments.items():
            if value[2] == 0:
                len1=value[1]
            else:
                len1=0
            if len1>len0:
                len0=len1
                chosen_key=key
        #Assigning the right similar id
        chosen_id=segments[chosen_key][0]
        for key, value in segments.items():
            if value[2] == 0 and key != chosen_key:
                value.append(chosen_id)
            else :
                value.append(0) #Will not consider any other segment and proceed normally


        #Lastly we assign wether to accelerate either qq (2) and qqdag (3) or both (1) or none (0) based on similar id
        for key, value in segments.items():
            if value[2] == 0 and key != chosen_key:
                value.append(3)
            else :
                value.append(0) #Will not consider any other segment and proceed normally


        mol_counter+=1
        with open(file_name, 'a') as file:

            file.write(f'Mol-'+str(mol_counter) + ','+str(n_segm)+'\n')
            for key, data in segments.items():        
                file.write(f'{" , ".join(map(str, data))}\n')
                    




        #Part where we record the generation number of each monomer
        indx=[]
        gens=[]
        for key, data in segments.items():
            
            typ=data[2]
            src=data[3]
            tar=data[4]
            gen=data[5]
            temp_path=nx.shortest_path(intg, source=src, target=tar, weight=None, method='dijkstra') #From detailed path
            inbetween_monomers=temp_path[1:-1]
            
            if typ==-1: #stem
                indx.append(src)
                indx=indx+inbetween_monomers
                indx.append(tar)
                gens=gens+(2+len(inbetween_monomers))*[gen]

            else: #terminal and internal
                
                if len(inbetween_monomers)!=0:
                    indx=indx+inbetween_monomers

                indx.append(tar)

                gens=gens+(1+len(inbetween_monomers))*[gen]
#            elif typ==1: #internal
#                
#                if len(inbetween_monomers)!=0:
#                    indx=indx+inbetween_monomers
#                indx.append(tar)
            

        gens=np.array(gens)
        indx=np.array(indx)
        sort_indices=np.argsort(indx)
        
        indx=indx[sort_indices]
        gens=gens[sort_indices]

        data_indx_gens = np.column_stack((indx, gens))
        fil_nam=fol2+'/mol_id_'+str(mol_counter)+'.dat'
        np.savetxt(fil_nam, data_indx_gens, delimiter=' ', header='# 1. Index 2. Generation', fmt='%i')
        #End of part where we record the generation number of each monomer             










    ########################################
    #### FOR VISUALISATION PURPOSES ONLY ###
    ########################################


    #Create a dictionary to map each number to its corresponding component color
#    component_colors = {}
#    for i, component in enumerate(components):
#        color = plt.cm.tab10(i % 10)  # Choose a color from a colormap
#        for number in component:
#            component_colors[number] = color



#    node_colors = ['blue' for node in G.nodes()]
#    i=0
#    node_colors=[]
#    for node in G.nodes():
#        if np.isin(node, monomer_indices):
#            node_colors.append("blue")
#        elif np.isin(node, ind_F):
#            node_colors.append("red")
#        elif np.isin(node, ind_C):
#            node_colors.append("yellow")
#    # Plot the graph with nodes color-coded by their connected components
#    pos = nx.nx_pydot.graphviz_layout(G)


#    #for i, component in enumerate(components):
#    #    shift = np.array([i*1.5,i*1.5])
#    #    for number in component:
#    #        pos[number] += shift


#    nx.draw(G, pos, node_color=node_colors, with_labels=False, font_weight='bold', node_size=100)
#    plt.savefig('plot.pdf')
#    plt.show()


##########################
##### SIZE ANALYSIS  #####
##########################

size_polymer=np.array(size_polymer)
average_size=np.average(size_polymer)
weight_average=np.average(size_polymer*size_polymer)/average_size

pdi=weight_average/average_size


weights=np.ones(len(size_polymer))/len(size_polymer)

with open(fol1+"/hyperbranched_data/ensemble_details.dat", "w") as file:
    
    file.write("For whole molecule\n")
    file.write("Average size " +str(average_size)+"\n")
    file.write("Weight Average " +str(weight_average)+"\n")
    file.write("Pdi " +str(pdi)+"\n")

    file.write("For ab monomers added\n")
    size=np.array(size_abmonomer)
    average_size=np.average(size)
    weight_average=np.average(size*size)/average_size
    pdi=weight_average/average_size
    file.write("Average size " +str(average_size)+"\n")
    file.write("Weight Average " +str(weight_average)+"\n")
    file.write("Pdi " +str(pdi)+"\n")




header="Weights n_i/n (ratio of number of chains) of each polymer found in the topology file in the same order"
np.savetxt(fol1+'/hyperbranched_data/weights.dat', weights,header=header,fmt='%1.15f')   






