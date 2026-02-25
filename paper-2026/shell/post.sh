#!/bin/bash
#SBATCH --gpus=1
#SBATCH --output=post%j.out
#SBATCH --error=post%j.err
#SBATCH --export=ALL

module load cuda/12.9
julia --project=~/run/TurbCircStat\
	~/run/TurbCircStat/paper-2026/postdata.jl\
    -p ~/run -m "w"
