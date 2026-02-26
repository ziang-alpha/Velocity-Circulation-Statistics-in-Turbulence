using CairoMakie, HDF5, LinearAlgebra
using FourierFlows, CUDA
include((@__DIR__) * "/../src/circulation.jl")

mytheme = merge(
	Theme(
		Axis = (
			xgridvisible = false,
			xminorticksvisible = true,
			xminorticksize = 2,
			xticksize = 3,
			xticklabelsize = 8,
			xlabelsize = 8,
			xlabelpadding = 2,
			xtickalign = 1,
			xminortickalign = 1,
			ygridvisible = false,
			yminorticksvisible = true,
			yminorticksize = 2,
			yticksize = 3,
			yticklabelsize = 8,
			ylabelsize = 8,
			ylabelpadding = 2,
			ytickalign = 1,
			yminortickalign = 1,
		),
		Lines = (
			joinstyle = :round,
			linewidth = 1,
		),
		Scatter = (
			markersize = 5,
		),
		Legend = (
			labelsize = 8,
			padding = 6,
			patchsize = (10, 5),
			colgap = 7,
			patchlabelgap = 3,
			framevisible = false,
		),
		Errorbars = (
			linewidth = 0.5,
			whiskerwidth = 2,
		),
		Text = (
			fontsize = 8,
		),
	),
	theme_latexfonts(),
)

## Global constants
pfile = h5open((@__DIR__) * "/postdata.h5", "r")
cm = 96/2.54
pt = 4/3
colormaps = cgrad(:roma, 8; categorical = true)
Reαβs = map(keys(pfile["wavenumber-energy_spectra"])) do str
	Re = parse(Float64, match(r"Re(\d+\.\d+)", str)[1])
	data = pfile["wavenumber-energy_spectra"][str]
	x, y = read(data["x"]), read(data["y"])
	αβ = hcat(y[1:300], x[1:300] .^ 2 .* y[1:300])\(2π .* x[1:300])
	(Re, αβ[1], αβ[2])
end

## Figure 1
with_theme(mytheme) do
	fig = Figure(size = (8.5cm, 8.2cm))

	# Layout
	begin
		ga = GridLayout(fig[1, 1])
		axa = Axis(ga[1, 1];
			xlabel = L"k",
			ylabel = L"E(k)",
			xscale = log10,
			yscale = log10,
			limits = (1, 512, 1e-6, 1),
			xticks = 2 .^ (0:2:8),
			xminorticks = IntervalsBetween(3),
			yminorticks = IntervalsBetween(9),
		)

		gb = GridLayout(fig[1, 2])
		axb = Axis(gb[1, 1];
			xlabel = L"u/\sqrt{⟨u^2⟩}",
			ylabel = "PDF",
			yscale = log10,
			limits = (-7, 7, 1e-10, 1),
			yminorticksvisible = false,
		)

		gc = GridLayout(fig[2, 1])
		axct = Axis(gc[1, 1];
			ylabel = L"\Pi_E/ε",
			xscale = log10,
			limits = (1, 512, -1, 0.5),
			xticklabelsvisible = false,
			xticks = 2 .^ (0:2:8),
			xminorticks = IntervalsBetween(3),
		)
		axcb = Axis(gc[2, 1];
			xlabel = L"k",
			ylabel = L"\Pi_\Omega/η",
			xscale = log10,
			limits = (1, 512, -0.5, 1),
			xticks = 2 .^ (0:2:8),
			xminorticks = IntervalsBetween(3),
		)
		rowgap!(gc, 8)

		gd = GridLayout(fig[2, 2])
		axd = Axis(gd[1, 1];
			xlabel = L"Re",
			ylabel = L"l_\text{eq}k_f",
			yscale = log10,
			limits = (3.2, 4.1, 1, 10^3),
			yminorticksvisible = false,
		)

		for (gl, lb) in zip([ga, gb, gc, gd], ["(a)", "(b)", "(c)", "(d)"])
			Label(gl[1, 1, TopLeft()], lb;
				fontsize = 10,
				font = :bold,
				padding = (0, 22, -10, 0),
			)
		end

		rowgap!(fig.layout, 10)
		colgap!(fig.layout, 10)
	end

	for (data, (Re, α, β), c) in zip(pfile["wavenumber-energy_spectra"], Reαβs, colormaps)
		x, y = read(data["x"]), read(data["y"])
		lines!(axa, x, y; color = c)
		if α > 0 && β > 0
			xref = x[4:300]
			eqSpec(x) = @. 2π * x / (α + β * x^2)
			lines!(axa, xref, eqSpec(xref); linestyle = :dashdot, color = :black)
			scatter!(axa,sqrt(α/β),eqSpec(sqrt(α/β)); color = :black)
		end
	end
	lines!(axa, 2:64, 2e-6 .* (2:64); linestyle = :dash, color = :red)
	lines!(axa, 2:64, 1e-1 .* (2:64) .^ (-1); linestyle = :dash, color = :blue)
	text!(axa, 16, 1e-5; text = L"k", color = :red)
	text!(axa, 16, 1e-2; text = L"k^{-1}", color = :blue)


	for (data, c) in zip(pfile["velocity-pdf:filtered"], colormaps)
		x, y, Δy = read(data["x"]), read(data["y"]), read(data["Δy"])
		std = sqrt(sum(@. y*x^2) * (x[2]-x[1]))
		scatter!(axb, x ./ std, y .* std; color = c)
		errorbars!(axb, x ./ std, y .* std, Δy .* std; color = c)
	end
	xref = -7:0.05:7
	yref = @. exp(-xref^2/2)/sqrt(2π)
	lines!(axb, xref, yref; linestyle = :dashdot, color = :black)

	for (efdata, wfdata, c) in zip(pfile["wavenumber-normalized_energy_fluxes"], pfile["wavenumber-enstrophy_fluxes"], colormaps)
		x1, y1 = read(efdata["x"]), read(efdata["y"])
		x2, y2 = read(wfdata["x"]), read(wfdata["y"])
		lines!(axct, x1, y1; color = c)
		lines!(axcb, x2, y2; color = c)
	end

	x = [Reαβ[1] for Reαβ in Reαβs if Reαβ[2] > 0 && Reαβ[3] > 0]
	y = [sqrt(Reαβ[3]/Reαβ[2])*1024 / 3sqrt(Reαβ[1]) for Reαβ in Reαβs if Reαβ[2] > 0 && Reαβ[3] > 0]
	lines!(axd, x, y; color = :black, alpha = 0.2)
	scatter!(axd, x, y; color = 1:7, colormap = colormaps)

	Legend(gd[1, 1],
		[LineElement(color = c, linestyle = nothing) for c in colormaps],
		[match(r"Re(\d+\.\d+)", str)[1] for str in keys(pfile["wavenumber-energy_spectra"])];
		height = Relative(0.45),
		width = Relative(0.65),
		halign = 0,
		valign = 1,
		nbanks = 4,
		orientation = :horizontal,
	)
	save((@__DIR__) * "/figures/diag.pdf", fig)
end

## Figure 2


with_theme(mytheme) do
	fig = Figure(size = (8.5cm, 8.2cm))

	# Layout
	begin
		ga = GridLayout(fig[1, 1])
		axa = Axis(ga[1, 1];
			xlabel = L"\sqrt{l_1l_2}/2\pi l_\text{eq}",
			ylabel = L"\sqrt{⟨Γ_C^2⟩/⟨Γ_□^2⟩}",
			xscale = log10,
			limits = (1e-2, 1e2, 1, 1.3),
		)

		gb = GridLayout(fig[1, 2])
		axb = Axis(gb[1, 1];
			xlabel = L"l_1/l_2",
			ylabel = L"(⟨Γ_C^p⟩/⟨Γ_□^p⟩)^{1/p}",
			xscale = log10,
			limits = (1e-2, 1, 1e-10, 1),
		)

		gc = GridLayout(fig[2, 1])
		axc = Axis(gc[1, 1];
			xlabel = L"(l_1 + l_2)/4\pi l_\text{eq}",
			ylabel = L"\sqrt{⟨Γ_C^2⟩/⟨Γ_□^2⟩}",
			xscale = log10,
			limits = (1e-2, 1e2, 0.6, 1.1),
		)

		gd = GridLayout(fig[2, 2])
		axd = Axis(gd[1, 1];
			xlabel = L"l_1/l_2",
			ylabel = L"(⟨Γ_C^p⟩/⟨Γ_□^p⟩)^{1/p}",
			xscale = log10,
			limits = (3.2, 4.1, 1, 10.0^1.5),
		)

		for (gl, lb) in zip([ga, gb, gc, gd], ["(a)", "(b)", "(c)", "(d)"])
			Label(gl[1, 1, TopLeft()], lb;
				fontsize = 10,
				font = :bold,
				padding = (0, 22, -10, 0),
			)
		end

		rowgap!(fig.layout, 10)
		colgap!(fig.layout, 10)
	end

	# for (data, (Re, α, β), c) in zip(pfile["loop_size-var_ratio:equal area"], Reαβs, colormaps)
	# 	if α > 0 && β > 0
	# 		leq = 2π * sqrt(β/α)
	# 		x, y = read(data["x"]) ./ leq, read(data["y"])
	# 		lines!(axa, x, y; color = c)
	# 	end
	# end

	for (data, (Re, α, β), c) in zip(pfile["loop_size-var_ratio:equal perimeter"], Reαβs, colormaps)
		if α > 0 && β > 0
			leq = 2π * sqrt(β/α)
			x, y = read(data["x"]) ./ leq, read(data["y"])
			lines!(axc, x, y; color = c)
		end
	end
	
	aaa(α,β) = begin
		ngrid = 1024
		grid = TwoDGrid(GPU(); nx = ngrid, Lx = 2π)
		loopsizes = (16:8:floor(Int, grid.nx/16))
		squares = [rectsloop(grid, l, l, 0, 0) for l in loopsizes]
		rects = [rectsloop(grid, 2 * l, l / 2, 0,0) for l in loopsizes]
		Eh = @. grid.Krsq/(α+β*grid.Krsq)
		Eh[grid.Krsq .> (ngrid/3)^2] .= 0
		CUDA.@allowscalar Eh[1,1] = 0
		varSquares = map(squares) do hsh
			FourierFlows.parsevalsum(Eh .* abs2.(hsh), grid)
		end

		varRects = map(rects) do hsh
			FourierFlows.parsevalsum(Eh .* abs2.(hsh), grid)
		end
		lines!(axa, loopsizes./ (2π *sqrt(β/α)), varRects ./ varSquares; color = :black)
	end
	autolimits!(axa)
	save((@__DIR__) * "/fig2.pdf", fig)
end
