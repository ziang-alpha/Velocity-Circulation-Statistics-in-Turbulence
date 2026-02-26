#!/bin/bash
#SBATCH --gpus=1
#SBATCH --output=sout/sim%j.out
#SBATCH --error=sout/sim%j.err
#SBATCH --export=ALL

module load cuda/12.9
julia --project=~/run/TurbCircStat\
	~/run/TurbCircStat/src/simulations/ns2.jl\
    --re $1 --ngrid $2 --ndata $3\
	--nstep $4 --gpu --path ~/run --cfl $5