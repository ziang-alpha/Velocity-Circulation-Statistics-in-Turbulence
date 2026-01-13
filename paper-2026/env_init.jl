using Pkg
Pkg.activate((@__DIR__) * "/..")
Pkg.update()
Pkg.resolve()

using GeophysicalFlows, Random, CUDA, HDF5, ProgressMeter