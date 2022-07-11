#!/bin/bash

# Always specify CPU and RAM resources needed, and walltime
# The default is one task per node
# not that --cpus-per-task option will change this default

#SBATCH --job-name="GrainMapping_gpu_slurm"                            # Job name
#SBATCH --mail-type=START,END,FAIL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=haixing.fang@grenoble-inp.fr     # Where to send mail	
#SBATCH --partition=nice                             # Run on NICE machine, limit of 12 h
# #SBATCH --partition=nice-long                        # Run on NICE-long machine, no limit on time
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=28
# #SBATCH --mem-per-cpu=10G                            # Job memory request
# #SBATCH --time=12:00:00                              # Time limit hrs:min:sec

# #use gpu
#SBATCH --nodes=1
#SBATCH --exclusive 
#SBATCH --tasks-per-node=2
#SBATCH --time=12:00:00                       # maximum is 12 h
#SBATCH --gres=gpu:1                          # number of gpu, one node has two gpus
#SBATCH -p gpu

#SBATCH --output=GM_slurm_%x.%j.out                         # %j job id; %x job name
#SBATCH --error=GM_slurm_%x.%j.err                          # error message
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
matlab -nodisplay -r "run_GrainMapping_fullvol_fun(4)"
# matlab -nodisplay -r "run_GrainMapping_fullvol_fun_continue(3)"
# #matlab -nodisplay -r "run_GrainMapping_subvol_fun(2)"    # run a job for a subvolume
# #matlab -nodisplay -r "run_GrainMapping_subvol_fun_continue(1)" # continue for a subvolume
# #matlab -nodisplay -r "test_grain_recon"
