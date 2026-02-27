using FourierFlows
# Compute the radial wavenumbers
" Radial wavenumbers for a given grid "
radialk(grid) = begin
	fh = device_array(grid)(zeros(grid.nkr, grid.nl))
	kr, _ = FourierFlows.radialspectrum(fh, grid; refinement = 1)
	kr
end
# Compute the radial spectrum
" Radial spectrum from given 2D spectrum "
radialspectrum(fh, grid) = begin
	_, fhr = FourierFlows.radialspectrum(fh, grid; refinement = 1)
	vec(fhr)
end
# Compute the cumulate radial spectrum (e.g. fluxes)
" Cumulative radial spectrum from given 2D spectrum "
cumradialspectrum(fh, grid) = begin
	kr, fhr = FourierFlows.radialspectrum(fh, grid; refinement = 1)
	cumsum(vec(fhr)) .* (kr[2] - kr[1])
end

# Diagnostic Quantities for 2D Navier Stokes equations
" 2D enstrophy spectrum for a single frame of vorcity "
enstrophy(ζh, grid) = abs2.(ζh) / (grid.nx * grid.ny)^2 / 2

" 2D energy spectrum for a single frame of vorcity "
energy(ζh, grid) = abs2.(ζh) .* grid.invKrsq / (grid.nx * grid.ny)^2 / 2

" 2D spectrum of the enstrophy injection rate of a single frame of vorticity "
enstrophyInjectionRate(ζh, grid) = begin
	ζx = grid.rfftplan \ (@. im * grid.kr * ζh)
	ζy = grid.rfftplan \ (@. im * grid.l * ζh)
	u = grid.rfftplan \ (@. im * grid.l * grid.invKrsq * ζh)
	v = grid.rfftplan \ (@. -im * grid.kr * grid.invKrsq * ζh)
	Nh = grid.rfftplan * (@. u * ζx + v * ζy)
	real.(Nh .* conj(ζh)) / (grid.nx * grid.ny)^2
end

" 2D spectrum of the enstrophy injection rate of a single frame of vorticity "
energyInjectionRate(ζh, grid) = begin
	ζx = grid.rfftplan \ (@. im * grid.kr * ζh)
	ζy = grid.rfftplan \ (@. im * grid.l * ζh)
	u = grid.rfftplan \ (@. im * grid.l * grid.invKrsq * ζh)
	v = grid.rfftplan \ (@. -im * grid.kr * grid.invKrsq * ζh)
	Nh = grid.rfftplan * (@. u * ζx + v * ζy)
	real.(Nh .* conj(ζh) .* grid.invKrsq) / (grid.nx * grid.ny)^2
end




