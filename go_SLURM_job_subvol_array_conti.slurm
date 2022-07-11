#!/bin/bash

# Always specify CPU and RAM resources needed, and walltime
# The default is one task per node
# not that --cpus-per-task option will change this default

#SBATCH --job-name="GrainMapping_slurm_s2"                            # Job name
#SBATCH --mail-type=START,END,FAIL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=haixing.fang@grenoble-inp.fr     # Where to send mail	
#SBATCH --partition=nice                             # Run on NICE machine, a limit of 12 h
# #SBATCH --partition=nice-long                        # Run on NICE-Long machinie, no limit of time
# #SBATCH --nodes=12                                 # Number of nodes
# #SBATCH --ntasks=12                                  # Number of processors (CPU)
# #SBATCH --ntasks-per-node=3                          # Number of tasks per node

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=10G                            # Job memory request
#SBATCH --time=12:00:00                              # Time limit hrs:min:sec
# #SBATCH --output=GrainMapping_s2_subvol1.log                # Standard output and error log
#SBATCH --output=GM_slurm_%x.%j.out                         # %j job id; %x job name
#SBATCH --error=GM_slurm_%x.%j.err                          # error message
#SBATCH --array=1-8                                  #1-4 for 1,2,3,4; 1-7:2 for 1,3,5,7; 0-15%4 for 4 maximum at the same time
# echo Running on $HOSTNAME : /home/esrf/haixing0a/Documents/Matlab_script/GrainRecon $(sed -n ${SLURM_ARRAY_TASK_ID}p subvol_slurm.params)
# #/home/esrf/haixing0a/Documents/Matlab_script/GrainRecon $(sed -n ${SLURM_ARRAY_TASK_ID}p subvol_slurm.params)

# Operations

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

export LD_LIBRARY_PATH=/home/esrf/haixing0a/Documents/Matlab_script/dip/Linuxa64/lib:$LD_LIBRARY_PATH
cd /home/esrf/haixing0a/Documents/Matlab_script/GrainRecon
# #srun matlab -nodisplay -r "run_GrainMapping_fun_continue_v3($SLURM_ARRAY_TASK_ID)"
srun matlab -nodisplay -r "run_GrainMapping_subvol_fun_continue($SLURM_ARRAY_TASK_ID)"

