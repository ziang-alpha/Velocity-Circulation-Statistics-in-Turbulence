#!/bin/bash 
BASE_DIR="$HOME/Project/Velocity-Circulation-Statistics-in-Turbulence"
julia --project="$BASE_DIR" "$BASE_DIR/paper-2026/simulation.jl" $1 $2 $3 $4