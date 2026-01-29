fid = h5open(args["path"] * "/diagnostics.h5", "w")

catvecArray(x...) = cat(vec.(Array.(x))..., dims = 2)

@time begin
	@info "Energy spectra."
	group = create_group(fid, "energy_spectra")
	for datum in raw_data
		kr = radialk(datum[:grid])
		# Compute the averaged 2D energy spectrum
		Ehrs = map(datum[:ζhs]) do ζh
			Eh = energy(ζh, datum[:grid])
			radialspectrum(Eh, datum[:grid])
		end
		mean, std = measure(Ehrs)
		group[datum[:name]] = catvecArray(kr, mean, std)
	end
end

@time begin
	@info "Energy fluxes."
	group = create_group(fid, "energy_fluxes")
	for datum in raw_data
		kr = radialk(datum[:grid])
		# Compute the averaged 2D energy spectrum
		Fhrs = map(datum[:ζhs]) do ζh
			Nh = energyInjectionRate(ζh, datum[:grid])
			cumradialspectrum(Nh, datum[:grid])
		end
		mean, std = measure(Fhrs)
		group[datum[:name]] = catvecArray(kr, mean, std)
	end
end

@time begin
	@info "Enstrophy fluxes."
	group = create_group(fid, "enstrophy_fluxes")
	for datum in raw_data
		kr = radialk(datum[:grid])
		# Compute the averaged 2D energy spectrum
		Fhrs = map(datum[:ζhs]) do ζh
			Nh = enstrophyInjectionRate(ζh, datum[:grid])
			cumradialspectrum(Nh, datum[:grid])
		end
		mean, std = measure(Fhrs)
		group[datum[:name]] = catvecArray(kr, mean, std)
	end
end

@time begin
	@info "Velocity PDFs."
	group = create_group(fid, "u_pdfs")
	for datum in raw_data
		grid = datum[:grid]

		u0 = grid.rfftplan \ (im * first(datum[:ζhs]) .* grid.invKrsq .* grid.l)
		umax = maximum(abs.(u0))
		bin = range(-3*umax, 3*umax, 100)
		bc = (bin[1:(end-1)] .+ bin[2:end]) / 2
		
        pdfs = map(datum[:ζhs]) do ζh
			uh = im * ζh .* grid.invKrsq .* grid.l
			u = grid.rfftplan \ uh
			pdf(u, bin)
		end
		mean, std = measure(pdfs)
		group[datum[:name]] = catvecArray(bc, mean, std)
	end
end

@time begin
	@info "Large-scale velocity PDF."
    
	group = create_group(fid, "u_pdfs_filtered")
	for datum in raw_data
		grid = datum[:grid]

        u0 = grid.rfftplan \ (im * first(datum[:ζhs]) .* grid.invKrsq .* grid.l)
		umax = maximum(abs.(u0))
		bin = range(-3*umax, 3*umax, 100)
		bc = (bin[1:(end-1)] .+ bin[2:end]) / 2
		pdfs = map(datum[:ζhs]) do ζh
			uh = im * ζh .* grid.invKrsq .* grid.l
			uh .*= datum[:fh]
			u = grid.rfftplan \ uh
			pdf(u, bin)
		end
		mean, std = measure(pdfs)
		group[datum[:name]] = catvecArray(bc, mean, std)
	end
end

close(fid)
