fid_diag = h5open(args["path"] * "/diagnostics.h5", "w")

@time begin
	@info "Energy spectra."
	group = create_group(fid_diag, "energy spectra")
	for datum in rawdata
		kr = radialk(datum.grid)
		# Compute the averaged 2D energy spectrum
		Ehrs = map(datum.ζhs) do ζh
			Eh = energy(ζh, datum.grid)
			radialspectrum(Eh, datum.grid)
		end
		mean, std = measure(Ehrs)
		group["Re$(datum.Re)/radialk"] = Array(kr)
		group["Re$(datum.Re)/mean"] = Array(mean)
		group["Re$(datum.Re)/std"] = Array(std)		
	end
end

@time begin
	@info "Energy fluxes."
	group = create_group(fid_diag, "energy fluxes")
	for datum in rawdata
		kr = radialk(datum.grid)
		# Compute the averaged 2D energy spectrum
		Fhrs = map(datum.ζhs) do ζh
			Nh = energyInjectionRate(ζh, datum.grid)
			cumradialspectrum(Nh, datum.grid)
		end
		mean, std = measure(Fhrs)
		group["Re$(datum.Re)/radialk"] = Array(kr)
		group["Re$(datum.Re)/mean"] = Array(mean)
		group["Re$(datum.Re)/std"] = Array(std)
	end
end

@time begin
	@info "Enstrophy fluxes."
	group = create_group(fid_diag, "enstrophy fluxes")
	for datum in rawdata
		kr = radialk(datum.grid)
		# Compute the averaged 2D energy spectrum
		Fhrs = map(datum.ζhs) do ζh
			Nh = enstrophyInjectionRate(ζh, datum.grid)
			cumradialspectrum(Nh, datum.grid)
		end
		mean, std = measure(Fhrs)
		group["Re$(datum.Re)/radialk"] = Array(kr)
		group["Re$(datum.Re)/mean"] = Array(mean)
		group["Re$(datum.Re)/std"] = Array(std)
	end
end

@time begin
	@info "Velocity PDFs."
	group = create_group(fid_diag, "velocity pdfs")
	for datum in rawdata
		grid = datum.grid

		u0 = grid.rfftplan \ (im * first(datum.ζhs) .* grid.invKrsq .* grid.l)
		umax = maximum(abs.(u0))
		bin = range(-3*umax, 3*umax, 100)
		bc = (bin[1:(end-1)] .+ bin[2:end]) / 2
		
        pdfs = map(datum.ζhs) do ζh
			uh = im * ζh .* grid.invKrsq .* grid.l
			u = grid.rfftplan \ uh
			pdf(u, bin)
		end
		mean, std = measure(pdfs)
		group["Re$(datum.Re)/bincenter"] = Array(bc)
		group["Re$(datum.Re)/mean"] = Array(mean)
		group["Re$(datum.Re)/std"] = Array(std)
	end
end

@time begin
	@info "Large-scale velocity PDFs."
	group = create_group(fid_diag, "filtered velocity pdfs")
	for datum in rawdata
		grid = datum.grid

		u0 = grid.rfftplan \ (im * first(datum.fζhs) .* grid.invKrsq .* grid.l)
		umax = maximum(abs.(u0))
		bin = range(-3*umax, 3*umax, 100)
		bc = (bin[1:(end-1)] .+ bin[2:end]) / 2
		
        pdfs = map(datum.fζhs) do ζh
			uh = im * ζh .* grid.invKrsq .* grid.l
			u = grid.rfftplan \ uh
			pdf(u, bin)
		end
		mean, std = measure(pdfs)
		group["Re$(datum.Re)/bincenter"] = Array(bc)
		group["Re$(datum.Re)/mean"] = Array(mean)
		group["Re$(datum.Re)/std"] = Array(std)
	end
end

close(fid_diag)
