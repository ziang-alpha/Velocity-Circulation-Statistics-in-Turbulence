using HDF5, ArgParse, Measurements
include((@__DIR__) * "/../src/postproc/circulation.jl")
include((@__DIR__) * "/../src/postproc/diagnostics.jl")
include((@__DIR__) * "/../src/postproc/statistics.jl")

aps = ArgParseSettings()
@add_arg_table! aps begin
	"--path", "-p"
	help = "Path of the raw data."
	arg_type = String
	default = (@__DIR__)
	"--mode", "-m"
	help = "Mode of the output .h5 file."
	arg_type = String
	default = "w"
	"--gpu"
	help = "Using GPU acceleration."
	action = :store_true
end
args = parse_args(aps)

dev = args["gpu"] ? GPU() : CPU()
filenames = readdir(args["path"] * "/.output/")

rfiles = map(filenames) do name
	h5open(args["path"] * "/.output/" * name, "r")
end

rawdata = map(rfiles, filenames) do fid, name
	m = match(r"Re(\d+\.\d+)_N(\d+)\.h5", name)
	Re, N = parse(Float64, m[1]), parse(Int, m[2])
	grid = TwoDGrid(dev; nx = N, Lx = 2π)

	ζhs = Iterators.map(fid) do ds
		ζh = read(ds)
		device_array(dev)(ζh)
	end
	
	fh = device_array(dev)(float(grid.Krsq .< (N / 3sqrt(Re) - 5)^2))
	fζhs = Iterators.map(fid) do ds
		ζh = read(ds) .* fh
		device_array(dev)(ζh)
	end

	Re => (ζhs = ζhs, fζhs = fζhs, grid = grid)
end
rawdata = Dict(rawdata)

pfile = h5open(args["path"] * "/postdata.h5", args["mode"])

include("datalist.jl")

close(pfile)
close.(rfiles)


