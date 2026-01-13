#!/bin/bash 
BASE_DIR="$HOME/run/Velocity-Circulation-Statistics-in-Turbulence"
export LD_LIBRARY_PATH="/data/home/sczd921/run/julia/julia-1.11.7/lib:${LD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="/data/home/sczd921/.julia/artifacts/7f39a18d94f87b5135df6731a327b61b8c463af6/lib:${LD_LIBRARY_PATH:-}"
julia --project="$BASE_DIR" "$BASE_DIR/paper-2026/simulation.jl" $1 $2 $3 $4
