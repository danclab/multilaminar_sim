#!/bin/bash

# ----- #
# Job name
#SBATCH --job-name=snr_sim

# ----- #
# Jobs to run; each element corresponds to a subject.
# Only run n at a time with (%n) at the end of the command!
# Counting the python way!
#SBATCH --array=0-799%150

# ----- #
# Computational resources.
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G

#SBATCH --exclude=ccwslurm0368,ccwslurm0370,ccwslurm0371,ccwslurm0372

# ----- #
# Task time limit (D-HH:MM:SS)
#SBATCH --time=0-12:00:00

#SBATCH --mail-user=jbonaiuto@gmail.com          # Where to send the e-mail
#SBATCH --mail-type=ALL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)

#SBATCH --licenses=sps                # Declaration of storage and/or software resources

# ----- #
# Python activation.
module add Programming_Languages/anaconda/3.11

# Activation of virtual python environment.
conda activate lameg_0.0.5

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

export FONTCONFIG_PATH=/etc/fonts
export FONTCONFIG_FILE=/etc/fonts/fonts.conf
export XDG_CACHE_HOME=/tmp/$USER/$SLURM_JOB_ID/fontcache
mkdir -p "$XDG_CACHE_HOME"
fc-cache -r >/dev/null 2>&1

# ----- #
# Run script.
# Standard output and standard error are NOT redirected to the same file.
python -u /pbs/home/b/bonaiuto/laminar_erf/pipeline_12_snr_interference_simulations.py > /sps/isc/bonaiuto/laminar_erf/output/output_snr_interference_$SLURM_ARRAY_TASK_ID.txt 2> /sps/isc/bonaiuto/laminar_erf/output/error_snr_interference_$SLURM_ARRAY_TASK_ID.txt ${SLURM_ARRAY_TASK_ID}

