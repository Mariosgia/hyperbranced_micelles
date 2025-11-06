########################################################
########################################################
###  SCRIPT FOR GENERATING HYPERBRANCHED POLYMERS ######
###       INITIAL STATES BASED ON THE PAPER       #######
###  "Linear-hyperbranched block copolymers       #######
###   consisting of polystyrene and dendritic     #######
###         poly(carbosilane) block"              #######
########################################################
########################################################

# F=Linear block C=Core grafted directly to F A= A bead in AB_2 B= B bead in AB_2 
# D= Center Helper bead for AB_2 E= Helper bead for AB_2 

#import hoomd
import gsd.hoomd
import math
import numpy as np
import matplotlib.pyplot as plt
import sys

def attractive_potential(r,d,c):
    
    res=-d*np.cos(r*math.pi/(2*c))*np.heaviside(c-r,0)
    return res

def attractive_force(r,d,c):
    
    res=-((d*math.pi) / (2.0*c) ) * np.sin(r*math.pi/(2*c))*np.heaviside(c-r,0)
    return res

def distance(x0, x1, dimensions):
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))

def helix(t,rr):
    
    x=rr*np.cos(t)
    y=rr*np.sin(t)
   
    return np.column_stack((x,y,t))


arguments = sys.argv
if len(arguments) != 3:
    print("Need two arguments nnF and nnC")


##########################
#### Polymer size parameters
########################

nnF=int(arguments[1]) #Size of linear block
nnC=int(arguments[2]) #Size of core 

num_bonds=nnF+nnC-1
###################
## Helix parameter 
###################
rr=5.0
rad=0.22 #Places a bead every 0.32 rads i.e. 1 distance unit away for rr=3.0
insertion_cutoff=2

#########
number_molecules=1
##########

#Parameters for topology and system

sizex=25 
sizey=sizex
sizez=sizex

ktt=1.0

#Topology
particles_id=np.ones(nnF,dtype='int8').tolist()+(np.ones(nnC,dtype='int8')*2).tolist()

bonds_id=np.ones(nnF-1,dtype='int8').tolist()+[2]+(np.ones(nnC-1,dtype='int8')*3).tolist()

bond_groups=np.column_stack((np.arange(nnF-1), np.arange(1, nnF))).tolist()+[[nnF-1,nnF]]+np.column_stack((np.arange(nnF,nnF+nnC-1), np.arange(nnF+1,nnF+nnC))).tolist()
bond_groups=np.array(bond_groups)



dimensions=np.array([sizex,sizey,sizez])

############################################################
######### CONSTRUCTING ATOMS,BONDS,ANGLES ETC. #############
############################################################

frame = gsd.hoomd.Frame()

# Place a polymer in the box.
frame.particles.N = nnF+nnC

t=np.arange(0, (nnF+nnC) * rad, rad)

particles_position=helix(t,rr)
particles_position=particles_position-np.array([0,0,9])
#print(particles_position)

frame.particles.position=particles_position.tolist()
frame.particles.types = ['F','C','A','B','D','E']
typeid_particles=[0]*nnF+[1]*nnC
frame.particles.typeid = typeid_particles
frame.configuration.box = [sizex, sizey, sizez, 0, 0, 0]



### Connect particles with bonds.
frame.bonds.N = num_bonds
frame.bonds.types = ['F-F','F-C','C-C','A-D','D-E','B-D']
typeid_bonds=[0]*(nnF-1)+[1]+[2]*(nnC-1)
frame.bonds.typeid = typeid_bonds
frame.bonds.group = bond_groups.tolist()



###Create  angle between A-B-C
#frame.angles.N = angles_per_mol*mol_N
frame.angles.types = ['A-B-D','D-B-D']
#frame.angles.typeid = [0,0,1] * mol_N


#frame.angles.group = (angles_group).tolist()



###Create dihedral 
frame.dihedrals.types = ['E-E-E-E']

#Adding the rest of the molecules

for i in range(1,number_molecules):

    frame.particles.N += nnF+nnC

    new_pos_error=1
    
    #Finds a random position to insert the new molecule
    while new_pos_error ==1 :
        random_pos=(np.random.rand(3)-0.5)*sizex/2
        new_pos=np.remainder(particles_position+random_pos+sizex/2.0,sizex)-sizex/2.0
    
        for j in range(0,nnF+nnC-1):
            dis_min=np.amin(distance(new_pos[j],frame.particles.position,dimensions))
            
            if dis_min<=insertion_cutoff:
#                    print("Distance between new molecule and existing molecules is smaller than cutoff")
#                    print("Minimum distance found is : "+str(dis_min))
                new_pos_error=1
                break
            else:
                new_pos_error=0

    frame.particles.typeid=frame.particles.typeid+typeid_particles
    frame.particles.position=frame.particles.position+(new_pos).tolist()

    frame.bonds.N += num_bonds

    frame.bonds.typeid=frame.bonds.typeid+typeid_bonds
    frame.bonds.group = frame.bonds.group+(bond_groups+(num_bonds+1)*i).tolist()                


with gsd.hoomd.open(name='molecular_init.gsd', mode='x') as f:
    f.append(frame)

