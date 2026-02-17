using Measurements

fid_aprule = h5open(args["path"] * "/areaperirules.h5", "w")

@time begin
    @info "Equal area variance ratio."
    group = create_group(fid_aprule, "eqarea var ratio")

    for datum in rawdata
        grid = datum.grid
        dl = 2π / grid.nx
        loopsizes = (10:10:floor(Int, grid.nx / 2)) .* dl
        squares = [rectsloop(grid, l, l, -π, -π) for l in loopsizes]
        rects = [rectsloop(grid, 2 * l, l / 2, -π, -π) for l in loopsizes]

        varSquares = map(squares) do hsh
            vars = map(datum.fζhs) do fζh
                Γ = getΓ(fζh, hsh, grid)
                moment(Γ, 2)
            end
            mean, std = measure(vars)
            mean ± std
        end

        varRects = map(rects) do hsh
            vars = map(datum.fζhs) do fζh
                Γ = getΓ(fζh, hsh, grid)
                moment(Γ, 2)
            end
            mean, std = measure(vars)
            mean ± std
        end

        varRatios = varRects ./ varSquares
        group["Re$(datum.Re)/loop sizes"] = Array(loopsizes)
        group["Re$(datum.Re)/values"] = Array(Measurements.value.(varRatios))
        group["Re$(datum.Re)/uncertainties"] = Array(Measurements.uncertainty.(varRatios))
    end
end

@time begin
    @info "Equal perimeter variance ratio."
    group = create_group(fid_aprule, "eqperi var ratio")

    for datum in rawdata
        grid = datum.grid
        dl = 2π / grid.nx
        loopsizes = (10:10:floor(Int, grid.nx / 2)) .* dl
        squares = [rectsloop(grid, l, l, -π, -π) for l in loopsizes]
        rects = [rectsloop(grid, 8l / 5, 2l / 5, -π, -π) for l in loopsizes]

        varSquares = map(squares) do hsh
            vars = map(datum.fζhs) do fζh
                Γ = getΓ(fζh, hsh, grid)
                moment(Γ, 2)
            end
            mean, std = measure(vars)
            mean ± std
        end

        varRects = map(rects) do hsh
            vars = map(datum.fζhs) do fζh
                Γ = getΓ(fζh, hsh, grid)
                moment(Γ, 2)
            end
            mean, std = measure(vars)
            mean ± std
        end

        varRatios = varRects ./ varSquares
        group["Re$(datum.Re)/loop sizes"] = Array(loopsizes)
        group["Re$(datum.Re)/values"] = Array(Measurements.value.(varRatios))
        group["Re$(datum.Re)/uncertainties"] = Array(Measurements.uncertainty.(varRatios))
    end
end

@time begin
    datum = rawdata[1]
    @info "Fixed area moments ratio with Re $(datum.Re)."
    group = create_group(fid_aprule, "eqarea mom ratio Re$(datum.Re)")
    grid = datum.grid

    dl = 2π / grid.nx
    fixarea = 8100dl^2
    loopsizes = [15; 20; 30; 40; 45; 60; 90] .* dl
    rects = [rectsloop(grid, l, fixarea / l, -π, -π) for l in loopsizes]

    momentsRatio = map(2:10) do order
        momRects = map(rects) do hsh
            moms = map(datum.fζhs) do fζh
                Γ = getΓ(fζh, hsh, grid)
                moment(Γ, order)
            end
            mean, std = measure(moms)
            mean ± std
        end
        momRatios = momRects ./ last(momRects)
        group["Re$(datum.Re) Order$(order)/aspect ratios"] = Array(loopsizes.^2 ./ fixarea)
        group["Re$(datum.Re) Order$(order)/values"] = Array(Measurements.value.(momRatios))
        group["Re$(datum.Re) Order$(order)/uncertainties"] = Array(Measurements.uncertainty.(momRatios))
    end

end

@time begin
    datum = rawdata[3]
    @info "Fixed perimeter moments ratio with Re $(datum.Re)."
    group = create_group(fid_aprule, "eqperi mom ratio Re$(datum.Re)")
    grid = datum.grid

    dl = 2π / grid.nx
    fixperi = 180dl
    loopsizes = (10:10:90) .* dl
    rects = [rectsloop(grid, l, fixperi -  l, -π, -π) for l in loopsizes]

    momentsRatio = map(2:10) do order
        momRects = map(rects) do hsh
            moms = map(datum.fζhs) do fζh
                Γ = getΓ(fζh, hsh, grid)
                moment(Γ, order)
            end
            mean, std = measure(moms)
            mean ± std
        end
        momRatios = momRects ./ last(momRects)
        group["Re$(datum.Re) Order$(order)/aspect ratios"] = Array(loopsizes ./ (fixperi .- loopsizes))
        group["Re$(datum.Re) Order$(order)/values"] = Array(Measurements.value.(momRatios))
        group["Re$(datum.Re) Order$(order)/uncertainties"] = Array(Measurements.uncertainty.(momRatios))
    end

end

close(fid_aprule)
