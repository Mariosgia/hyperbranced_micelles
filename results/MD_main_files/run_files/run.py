########################################################
########################################################
### SIMULATION SCRIPT FOR HYPERBRANCHED POLYMERS #######
###       FORMATION USING PAIR POTENTIALS        #######
###         CREATED BY MARIOS GIANNAKOU          #######
########################################################
########################################################

# F=Linear block C=Core grafted directly to F A= A bead in AB_2 B= B bead in AB_2 
# D= Center Helper bead for AB_2 E= Helper bead for AB_2 

import hoomd
import gsd.hoomd
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma, factorial
from scipy.interpolate import interp1d
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

def inverse_transform_sampling(x_cumulative, cumulative_probs, size):
    # Create an interpolation function for the inverse CDF
    inverse_cdf = interp1d(cumulative_probs, x_cumulative, kind='linear')

    # Generate uniform random numbers
    u = np.random.uniform(0, 1, size)

    # Apply the inverse transform
    return inverse_cdf(u)

def schulz_zimm(x,k,mu):
    return np.power(k/mu,k)*np.power(x,k-1)/gamma(k)*np.exp(-k*x/mu)



arguments = sys.argv
if len(arguments) != 3:
    print("Need two arguments mu and dw")


#Distribution details
mu=float(arguments[1])
dw=float(arguments[2])
cutoff=0.00001


seedd=np.random.randint(1,20529)
    

if ( abs(dw-1) <0.0001): #Monodisperse
    final_N=int(mu)
else:
    k=1.0/(dw-1)
    xx = np.linspace(0,12*mu,10000)  # Example x values
    probs=schulz_zimm(xx,k,mu)

    xx_temp=xx
    probs_temp=probs

    probs=probs_temp[ ~((xx_temp>mu) & (probs_temp<cutoff)) ]
    xx=xx_temp[ ~((xx_temp>mu) & (probs_temp<cutoff)) ]

    print( " Max Length of distribution : " +str(np.max(xx)))

    cumulative_probs = np.cumsum(probs) # Corresponding cumulative probabilities
    cumulative_probs/=cumulative_probs[-1]
    size = 1  
    custom_samples = inverse_transform_sampling(xx, cumulative_probs, size)
    final_N=int(custom_samples[0])




print("Monomers to be added:" +str(final_N)) #Add this many AB_2 monomers by the end of the simulation


dev=hoomd.device.CPU()
every_timestep=1e4 #Run simulations of so many timesteps until monomer has bonded


mol_N=0

#Parameters for attractive potential
depth=100.0
cut_distance=0.5
ktt=1.0


#Parameters for the topology of AB_2 monomers 
insertion_cutoff=1.5 #Cutoff when adding monomers
beads_per_mol=8
bonds_per_mol=7
angles_per_mol=3
dihedrals_per_mol=1

particles_id=[2,4,5,5,5,5,3,3] #['F','C','A','B','D','E']
bonds_id=[3,4,4,4,4,5,5] #['F-F','C-C','F-C','A-D','D-E','B-D']
angles_id=[0,0,1] # ['A-B-D','D-B-D']
dihedrals_id=[0] # ['E-E-E-E']

#Values used when adding new monomer later in the code
particles_position=np.array([[0,0,0],[0,1,0],[1,1,0],[-1,1,0],[0,1,1],[0,1,-1],[-np.cos(math.pi*1.0/3.0),1+np.sin(math.pi*1.0/3.0),0],[np.cos(math.pi*1.0/3.0),1+np.sin(math.pi*1.0/3.0),0]])
bond_groups=np.array([[0, 1], [1, 2],[1, 3],[1, 4],[1, 5],[1, 6],[1, 7]])
angles_group=np.array([[0,1,6],[0,1,7],[7,1,6]])
dihedrals_group=np.array([2,3,4,5])


#Importing init state
frame = gsd.hoomd.open(name='molecular_init.gsd', mode='r')
frame = frame.__getitem__(len(frame)-1)
number_particles=frame.particles.N

sizex=frame.configuration.box[0]
sizey=frame.configuration.box[1]
sizez=frame.configuration.box[2]
dimensions=frame.configuration.box[0:3]

particles_types=frame.particles.types
bonds_types=frame.bonds.types
angles_types=frame.angles.types
dihedrals_types=frame.dihedrals.types



#############################################
######### CONSTRUCTING POTENTIALS ###########
#############################################



## Apply the harmonic potential on the bonds.
harmonic = hoomd.md.bond.Harmonic()
harmonic.params['F-F'] = dict(k=100, r0=1.0)
harmonic.params['F-C'] = dict(k=100, r0=1.0)
harmonic.params['C-C'] = dict(k=100, r0=1.0)
harmonic.params['A-D'] = dict(k=100, r0=1.0)
harmonic.params['B-D'] = dict(k=100, r0=0.25)
harmonic.params['D-E'] = dict(k=100, r0=0.25)



nl = hoomd.md.nlist.Cell(buffer=0.2,exclusions=('bond',))

lj=hoomd.md.pair.LJ(nlist=nl)
for i in frame.particles.types:
    for j in frame.particles.types:
        
        lj.params[(i, j)] = dict(sigma=1.0, epsilon=1)
        lj.r_cut[(i, j)] = 2**((1.0/6.0))

lj.params[('C', 'A')] = dict(sigma=0, epsilon=0)
lj.r_cut[('C', 'A')] = 0
lj.params[('A', 'C')] = dict(sigma=0, epsilon=0)
lj.r_cut[('A', 'C')] = 0

lj.params[('A', 'B')] = dict(sigma=0, epsilon=0)
lj.r_cut[('A', 'B')] = 0
lj.params[('B', 'A')] = dict(sigma=0, epsilon=0)
lj.r_cut[('B', 'A')] = 0


##### Constructing and applying attractive potential

step=cut_distance/1000.0
rrs=np.arange(0,cut_distance+step,step)

attr_potential=attractive_potential(rrs,depth,cut_distance)
attr_force=attractive_force(rrs,depth,cut_distance)

attractive=hoomd.md.pair.Table(nl, default_r_cut=cut_distance)

for i in frame.particles.types:
    for j in frame.particles.types:

        attractive.params[(i, j)] = dict(r_min=0, U=[0], F=[0])

attractive.params[('C', 'A')] = dict(r_min=0, U=attr_potential, F=attr_force)
attractive.params[('A', 'C')] = dict(r_min=0, U=attr_potential, F=attr_force)
attractive.params[('A', 'B')] = dict(r_min=0, U=attr_potential, F=attr_force)
attractive.params[('B', 'A')] = dict(r_min=0, U=attr_potential, F=attr_force)

## Apply the COSINE potential on the angle B-A-B.
cosinesq = hoomd.md.angle.CosineSquared()
cosinesq.params['A-B-D'] = dict(k=100.0, t0=2.5/3.0*math.pi)
cosinesq.params['D-B-D'] = dict(k=100.0, t0=2.0/3.0*math.pi)

## Apply dihedral potential
opls = hoomd.md.dihedral.OPLS()
opls.params['E-E-E-E'] = dict(k1=0.0, k2=25.0, k3=0.0, k4=0.0)



#############################################
######### SETTING UP SIMULATION #############
#############################################


##### Perform the first simulation.
with gsd.hoomd.open(name='molecular_trajectory.gsd', mode='w') as f:
    f.append(frame)

sim = hoomd.Simulation(device=dev, seed=seedd)
sim.create_state_from_gsd(filename='molecular_init.gsd')
langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=ktt)
integrator = hoomd.md.Integrator(dt=0.005,
                                 methods=[langevin],
                                 forces=[harmonic,lj]) 
gsd_writer = hoomd.write.GSD(filename='molecular_trajectory.gsd',
                             trigger=hoomd.trigger.Periodic(int(1e3)),
                             mode='ab',dynamic=['property','attribute','topology','momentum'])


sim.operations.integrator = integrator
sim.operations.writers.append(gsd_writer)


sim.run(1e3)



################################################
#####    ADDING NEW MOLECULE AT RANDOM   #######
################################################
snapshot = sim.state.get_snapshot()

tt=0
if snapshot.communicator.rank == 0:

    while mol_N<final_N:


        #Getting snapshot and making a frame out of it
        snapshot = sim.state.get_snapshot()
        temp_frame=gsd.hoomd.Frame()


        temp_frame.particles.N=snapshot.particles.N
        temp_num_particles=temp_frame.particles.N
        temp_frame.particles.position=snapshot.particles.position
        temp_frame.particles.types = particles_types
        temp_frame.particles.typeid=snapshot.particles.typeid
        temp_frame.particles.velocity=snapshot.particles.velocity


        temp_frame.configuration.box = snapshot.configuration.box

        temp_frame.bonds.N = snapshot.bonds.N
        temp_num_bonds=temp_frame.bonds.N 
        temp_frame.bonds.types = bonds_types
        temp_frame.bonds.typeid = snapshot.bonds.typeid
        temp_frame.bonds.group = snapshot.bonds.group


        temp_frame.angles.N = snapshot.angles.N
        temp_num_angles=temp_frame.angles.N 
        temp_frame.angles.types = angles_types
        temp_frame.angles.typeid = snapshot.angles.typeid
        temp_frame.angles.group = snapshot.angles.group


        temp_frame.dihedrals.N = snapshot.dihedrals.N
        temp_num_dihedrals=temp_frame.dihedrals.N 
        temp_frame.dihedrals.types = dihedrals_types
        temp_frame.dihedrals.typeid = snapshot.dihedrals.typeid
        temp_frame.dihedrals.group = snapshot.dihedrals.group


     
     
        ##########################################
        ########### ADDING NEW MOLECULE ##########
        ##########################################

        temp_frame.particles.N=temp_frame.particles.N+beads_per_mol

        new_pos_error=1
        
        #Finds a random position to insert the new molecule
        while new_pos_error ==1 :
            random_pos=(np.random.rand(3)-0.5)*sizex
            new_pos=np.remainder(particles_position+random_pos+sizex/2.0,sizex)-sizex/2.0
        
            for i in range(0,beads_per_mol):
                dis_min=np.amin(distance(new_pos[i],temp_frame.particles.position,dimensions))
                
                if dis_min<=insertion_cutoff:
#                    print("Distance between new molecule and existing molecules is smaller than cutoff")
#                    print("Minimum distance found is : "+str(dis_min))
                    new_pos_error=1
                    break
                else:
                    new_pos_error=0
                
                
        new_vel=np.random.normal(0.0,ktt,(beads_per_mol,3))

                
    
        temp_frame.particles.position = np.vstack((temp_frame.particles.position,new_pos))
        temp_frame.particles.typeid=np.concatenate((temp_frame.particles.typeid,np.array(particles_id)))
        temp_frame.particles.velocity=np.vstack((temp_frame.particles.velocity,new_vel))
        
#        print(temp_frame.particles.position)
#        temp_frame.particles.position[-beads_per_mol:]=np.random.normal(0.0,ktt,(beads_per_mol,3))
#        print(temp_frame.particles.position)

    ##    ### Connect particles with bonds.
        temp_frame.bonds.N = temp_frame.bonds.N+bonds_per_mol 
        temp_frame.bonds.typeid=np.concatenate((temp_frame.bonds.typeid,np.array( bonds_id )))
        temp_frame.bonds.group =np.vstack((temp_frame.bonds.group,bond_groups+temp_num_particles))



    ##    ###Create  angle between A-B-C

        temp_frame.angles.N = temp_frame.angles.N+angles_per_mol 

        temp_frame.angles.typeid=np.concatenate((temp_frame.angles.typeid,np.array( angles_id )))

        temp_frame.angles.group =np.vstack((temp_frame.angles.group,angles_group+temp_num_particles))



#    ##    ###Create dihedral 

        temp_frame.dihedrals.N = temp_frame.dihedrals.N+dihedrals_per_mol 
        temp_frame.dihedrals.typeid=np.concatenate((temp_frame.dihedrals.typeid,np.array(dihedrals_id )))
        temp_frame.dihedrals.group =np.vstack((temp_frame.dihedrals.group,dihedrals_group+temp_num_particles))


        with gsd.hoomd.open(name='mol_addition.gsd', mode='a') as f:
            f.append(temp_frame)
            
    
        sim = hoomd.Simulation(device=dev, seed=seedd)
        sim.create_state_from_gsd(filename='mol_addition.gsd')

        integrator = hoomd.md.Integrator(dt=0.005,
                                         methods=[langevin],
                                         forces=[harmonic,cosinesq,opls,lj,attractive]) 
        gsd_writer = hoomd.write.GSD(filename='molecular_trajectory.gsd',
                                     trigger=hoomd.trigger.Periodic(int(1e6)),
                                     mode='ab',dynamic=['property','attribute','topology','momentum'])


        sim.operations.integrator = integrator
        sim.operations.writers.append(gsd_writer)

        added=0

        #Running simulation until added particle has been added
        while (added == 0 ):

            sim.run(every_timestep)
            tt+=every_timestep

            #Getting snapshot and making a frame out of it
            snapshot = sim.state.get_snapshot()
    #            temp_frame=gsd.hoomd.Frame()

    #            temp_frame.particles.typeid=


            positions=snapshot.particles.position
            posA_last=positions[snapshot.particles.typeid==2][-1]
            posB=positions[snapshot.particles.typeid==3]
            posC=positions[snapshot.particles.typeid==1]

            pos_temp=np.vstack((posB,posC))
            dis=distance(posA_last,pos_temp,dimensions)
            dis_argmin=np.argmin(dis)
            dis_min=dis[dis_argmin]

            if dis_min < cut_distance :
               print("Connected at : "+ str(tt))
               added=1
               with open("addition_details.dat", "a") as file:
                   file.write(str(mol_N+1)+"\t"+str(tt)+"\n")
                






        print("Added new molecule")
        print("Number of molecules: "+ str(mol_N+1))
        mol_N=mol_N+1

#########################################################
####    RUNNING LONG SIMULATION TILL EQUILIBRATION   ####
#########################################################

#frame = gsd.hoomd.open(name='molecular_final.gsd', mode='w')
print("Finished addition")

gsd_writer = hoomd.write.GSD(filename='final_configuration.gsd',
                             trigger=hoomd.trigger.Periodic(int(1)),
                             mode='ab',dynamic=['property','attribute','topology','momentum'])


sim.operations.writers.append(gsd_writer)



sim.run(1)
#sim.run(1e6)
#print("Finished script")

