using FourierFlows
"""
	Build a loop by combining a sequence of rectangular loops. 
	'heights', 'widths', 'startxs', 'startys' specify the height, width
	and the location left-bottom corner of each rectangle.
"""
function rectsloop(grid, heights, widths, startxs, startys)
	hs = zeros(Bool, grid.nx, grid.ny)
	for (h, w, x, y) in zip(heights, widths, startxs, startys)
		hs .= hs .|| rectangle(h, w, x, y, grid)
	end
	hs = device_array(grid)(hs)
	return grid.rfftplan * hs

end

function rectangle(height::Int, width::Int, startx::Int, starty::Int, grid)
	x = Bool.(vcat(zeros(startx-1), ones(width), zeros(grid.nx+1-startx-width)))
	y = Bool.(vcat(zeros(starty-1), ones(height), zeros(grid.ny+1-starty-height)))
	return x .&& transpose(y)
end

"""
	Compute the velocity circulation
	using the convolution of the vorcity field 'ζ' and the loop Heaviside 'hs'.
	'ζh' and 'hsh' are their Fourier transforms respectively.
"""
function getΓ(ζh, hsh, grid)
	Γh = device_array(grid)(ζh .* hsh) * (grid.dx * grid.dy)
	return grid.rfftplan \ Γh
end

"""
	Functional Π^α[C] defined as,
	Π^α[C] = ∬ₛ∬ₛ (-⧊)^(α/2)δ(x-x')dσdσ'
	which proportional to the area for α = 0, and to perimeter for α = 2
"""
periarea(hsh, α, grid) = FourierFlows.parsevalsum(abs2.(hsh) .* grid.Krsq .^ (α / 2), grid)


