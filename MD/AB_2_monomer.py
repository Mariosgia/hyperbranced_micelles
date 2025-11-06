import numpy as np
import gsd.hoomd
import math


# Parameters for the molecule
beads_per_mol = 8
bonds_per_mol = 7
angles_per_mol = 3
dihedrals_per_mol = 1

# Particle types (IDs), bonds, angles, dihedrals
particles_id = [2, 4, 5, 5, 5, 5, 3, 3]
bonds_id = [3, 4, 4, 4, 4, 5, 5]
angles_id = [0, 0, 1]
dihedrals_id = [0]

sizex=10
sizey=10
sizez=10

# Particle positions
particles_position = np.array([[0, 0, 0],
                               [0, 1, 0],
                               [1, 1, 0],
                               [-1, 1, 0],
                               [0, 1, 1],
                               [0, 1, -1],
                               [-np.cos(math.pi * 1.0 / 3.0), 1 + np.sin(math.pi * 1.0 / 3.0), 0],
                               [np.cos(math.pi * 1.0 / 3.0), 1 + np.sin(math.pi * 1.0 / 3.0), 0]])

# Bond groups
bond_groups = np.array([[0, 1],
                        [1, 2],
                        [1, 3],
                        [1, 4],
                        [1, 5],
                        [1, 6],
                        [1, 7]])

# Angles groups
angles_group = np.array([[0, 1, 6],
                         [0, 1, 7],
                         [7, 1, 6]])

# Dihedrals groups
dihedrals_group = np.array([2, 3, 4, 5])

frame = gsd.hoomd.Frame()

# Place a polymer in the box.
frame.particles.N = beads_per_mol




frame.particles.position=particles_position.tolist()
frame.particles.types = ['F','C','A','B','D','E']
frame.particles.typeid = particles_id
frame.configuration.box = [sizex, sizey, sizez, 0, 0, 0]



### Connect particles with bonds.
frame.bonds.N = bonds_per_mol
frame.bonds.types = ['F-F','F-C','C-C','A-D','D-E','B-D']
frame.bonds.typeid = bonds_id
frame.bonds.group = bond_groups.tolist()

with gsd.hoomd.open(name='ab2_monomer.gsd', mode='x') as f:
    f.append(frame)
