This repository contains all the necessary files for running the analysis and generation of all the results found in the article named:

"Theoretical study on micelle forming Linear-Hyperbranched block copolymers: A comparison between random and well-defined polymer architectures" 

The analysis and production of the figures can be found in the folder 'figures', while the necessary raw data along with the simulations scripts can be found in the "results" folder. 

The simulation scripts for the MD simulations are found in "results/MD_main_files/" while the SCF simulation scripts can be found in "results/SCF_main_files/". 

For the former case additional python scripts that export the topology of the produced LHBC polymers are included as well as some bash files that automate the running of the simulations. 

In the latter the "results/SCF_main_files/Main_Hyper" includes the Makefile that compiles several modules from "results/SCF_main_files/A*" and "results/SCF_main_files/Main_Hyper/main_hyper.f90" to produce "main_hyper.exe" which can be done using "make main_hyper". This is the executable used for the SCF simulations. Also contained in "results/SCF_main_files/Main_Hyper/analysis" which contains "analysis.f90" and is used for further data production of the converged states. This can be compiled using "make analysis". 

SCF simulation guide:

The SCF simulations as mentioned before are run using  "results/SCF_main_files/Main_Hyper/main_hyper.exe" which takes input parameters using the files, "input_general", "input_interactions" , "hyperbranched_data/topology_data.dat" and  "hyperbranched_data/weights.dat".

In "input_general" , parameters such as the size of the simulation box, the number of grid points, mixing parameters, initialization parameters, constraint parameters and ensemble parameters are chosen here.

In "input_interactions" parameters such as the $\chi \bar{N}$ and relative volumes of the monomers ($v_{\alpha}/v^{*}$) and the solvent ($v_{S}/(v^{*}\bar{N})$) are given.

The rest files, "hyperbranched_data/topology_data.dat" and  "hyperbranched_data/weights.dat",  contain information on the topology and proportion of chains $n_{i}/n_{T}$ ( They need to add up to 1).

The "topology_data.dat" file has the following form:

<pre>Number_Molecules, 3 
Mol-1,20
1 , 3 , 0 , 97 , 125 , 7 , 2 , 5 , 3
2 , 4 , 0 , 97 , 131 , 7 , 2 , 5 , 3
..
..
20 , 84 , -1 , 0 , 83 , 0 , 1 , 0 , 0
Mol-2,38
1 ....
..
..
38
Mol-3,10
1
..
..
10 </pre>

The first two lines are comments and are ignored and the third line is given as "Number_Molecules, [Number of molecule types]"
where [Number of molecule types] is the number of molecule types.

For each molecule a line starting with "Mol-[Mol id],[Number of blocks]" needs to be given.The number of such lines needs to match [Number of molecule types]. Below each such line the details for the topology of the polymer molecule are given in [Number of blocks] rows and 9 columns.

The 1st column is the id of the block and must be distinct number from 1 to [Number of blocks].
The 2nd column is the number of propagation steps that are taken for such block and needs to be above 1.
The 3rd column is the type of the block and can either be stem (-1), terminal (0) or internal (1). Their purpose is to define the propagation direction and initial condition properties, thus only one block can serve as stem.
The 4th and 5th columns refer to the source and target nodes refer to the branching points found on the polymer molecule. Each branching point is given a unique number (could be any integer) and nodes of each block are given such that source to target nodes are defined as the "forward" direction.
The 6th column is the generation number and we take the stem to have the lowest generation number (zero in all our case).
The 7th column takes a number which defines its monomer type and the number follows the convention found in "input_interactions".
The 8th and 9th columns are there to save time when computing the propagators of certain blocks in a polymer that obey symmetries by copying the propagator of "similar" blocks. If uncertain about how to use, simply set both to 0, which means the forward and backward propagators are calculated explicitly. If you want to copy the propagators then set the 8th column to the block id of the block that you wish to copy from. In the 9th column then you either choose between 1,2 or 3 for both, only forward or only backward propagators to be copied respectivelly. For example in the case of LDBC only for each generation only the propagators of one of the blocks for each generation needs to be calculated explicitly, while the rest of the block in that generation can be copied. An example of generation 2 LDBC topology file is shown below which demonstrates this point:

<pre># 0.Id 1.Length 2.Type 3.Source node 4.Target node 5.Generation 6. Monomer_type 7.Similar to id 8.Accelerate
#Mol-* 1.No. segments
Number_Molecules, 1
Mol-1, 7 
1 , 84 , -1 , 0 , 1 , 0 , 1 , 0 , 0
2 , 14 , 1 , 1 , 2 , 1 , 2 , 0 , 0
3 , 14 , 1 , 1 , 3 , 1 , 2 , 2 , 1
4 , 14 , 0 , 2 , 4 , 2 , 2 , 0 , 0
5 , 14 , 0 , 2 , 5 , 2 , 2 , 4 , 1
6 , 14 , 0 , 3 , 6 , 2 , 2 , 4 , 1
7 , 14 , 0 , 3 , 7 , 2 , 2 , 4 , 1</pre>

Also in LHBCs the backward propagators of the terminal blocks can be copied from the longest (most propagation steps) backward propagator.

Lastly, note that every propapagation step uses this a single "ds" value and this values is determined based on the weights and the collective number of propagation steps in each polymer.
