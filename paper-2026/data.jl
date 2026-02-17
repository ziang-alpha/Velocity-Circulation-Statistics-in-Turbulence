using HDF5, ArgParse
include((@__DIR__) * "/../src/circulation.jl")
include((@__DIR__) * "/../src/diagnostics.jl")
include((@__DIR__) * "/../src/statistics.jl")

aps = ArgParseSettings()
@add_arg_table! aps begin
    "--path", "-p"
    help = "Path of the raw data."
    arg_type = String
    default = (@__DIR__)
    "--gpu"
    help = "Using GPU acceleration."
    action = :store_true
    "--diags"
    help = "Compute the diagnostic data, including spectra and velocity pdfs."
    action = :store_true
    "--aprule"
    help = "Compute the data of area and perimeter rules."
    action = :store_true
end
args = parse_args(aps)

dev = args["gpu"] ? GPU() : CPU()
filenames = readdir(args["path"] * "/.output/")

rdfiles = map(filenames) do name
    h5open(args["path"] * "/.output/" * name, "r")
end

rawdata = map(rdfiles, filenames) do fid, name
    ζhs = begin
        nframe = maximum(parse.(Int, keys(fid)))
        Iterators.map(1:nframe) do n
            ζh = read(fid, "$(n)")
            device_array(dev)(ζh)
        end
    end
    regex = r"Re(\d+\.\d+)_N(\d+)\.h5"
    m = match(regex, name)
    Re, N = parse(Float64, m[1]), parse(Int, m[2])
    grid = TwoDGrid(dev; nx=N, Lx=2π)
    fh = device_array(dev)(float(grid.Krsq .< (N / 3sqrt(Re) - 5)^2))
    fζhs = Iterators.map(ζhs) do ζh
        ζh .* fh
    end

    (Re=Re, ζhs=ζhs, fζhs=fζhs, grid=grid)
end


args["diags"] && include("diagnostics.jl")
args["aprule"] && include("areaperi_rule.jl")

close.(rdfiles)
