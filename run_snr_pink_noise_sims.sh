#!/bin/bash

# ----- #
# Job name
#SBATCH --job-name=snr_sim


# ----- #
# Jobs to run; each element corresponds to a subject.
# Only run n at a time with (%n) at the end of the command!
# Counting the python way!
# SBATCH --array=0-2099%150
# SBATCH --array=0-699%200
#SBATCH --array=0-899%150


# ----- #
# Computational resources.
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G

# Instead of specifying a nodelist which will ask for all the nodes to be
# available for each job, exclude the nodes that are not contained in the
# nodelist. Each job will only occupy one node this way.
# SBATCH --nodelist=node[13-21]
# SBATCH --exclude=node[2-12]


# ----- #
# Task time limit (D-HH:MM:SS)
#SBATCH --time=1-00:00:00


# ----- #
# Output and error filenames.
# Currently skipped and instead used directly when calling the python script.
# --output=fichier_de_sortie${SLURM_ARRAY_TASK_ID}.txt
# --error=sortie_erreur.err

#SBATCH --licenses=sps


# ----- #
# Python activation.
module add Programming_Languages/anaconda/3.11

# Activation of virtual python environment.
conda activate lameg

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
#export JAVA_TOOL_OPTIONS="-Djava.awt.headless=true"

export FONTCONFIG_PATH=/etc/fonts
export FONTCONFIG_FILE=/etc/fonts/fonts.conf
export XDG_CACHE_HOME=/tmp/$USER/$SLURM_JOB_ID/fontcache
mkdir -p "$XDG_CACHE_HOME"
fc-cache -r >/dev/null 2>&1

# ----- #
# Run script.
# Standard output and standard error are NOT redirected to the same file.
python -u /pbs/home/b/bonaiuto/laminar_erf/pipeline_10_snr_pink_noise_simulations.py > /sps/isc/bonaiuto/laminar_erf/output/output_snr_pink_$SLURM_ARRAY_TASK_ID.txt 2> /sps/isc/bonaiuto/laminar_erf/output/error_snr_pink_$SLURM_ARRAY_TASK_ID.txt ${SLURM_ARRAY_TASK_ID}

