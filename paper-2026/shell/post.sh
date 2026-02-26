#!/bin/bash
#SBATCH --gpus=1
#SBATCH --output=sout/post%j.out
#SBATCH --error=sout/post%j.err
#SBATCH --export=ALL

module load cuda/12.9
julia --project=~/run/TurbCircStat\
	~/run/TurbCircStat/paper-2026/figures/getdata.jl\
    -p ~/run -m "w"
