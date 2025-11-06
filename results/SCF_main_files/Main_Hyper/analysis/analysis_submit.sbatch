#!/bin/bash
# Access the environment variable
threads=$(grep 'threads' input_general | awk '{print $1}')
echo $threads
# Use a here document to embed the SBATCH directives
sbatch <<EOT
#!/bin/bash
#SBATCH --partition=parallel                  # Queue name 'smp' or 'parallel' on Mogon II
#SBATCH --account=nhr-komet1schmid
#SBATCH -J scf_parallel           # Job name
#SBATCH -o scf.%j.out    # Specify stdout output file (%j expands to jobId)
#SBATCH -n 1                     # Total number of tasks
#SBATCH -N 1                     # Total number of nodes requested 
#SBATCH -c ${threads} 
#SBATCH --mem-per-cpu 4G
#SBATCH --time=0:02:00
#SBATCH --output=log.log


module load numlib/FFTW/3.3.10-GCC-12.3.0
export OMP_NUM_THREADS=${threads}
ulimit -s unlimited

srun -c ${threads} ./analysis.exe 
EOT
