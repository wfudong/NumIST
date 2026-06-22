export fourier_wavenumbers,
       etdrk4_coefficients,
       solve_kdv_etdrk4,
       interpolate_periodic,
       sampled_invariants,
       invariant_table

using FFTW
using Statistics

"""
    fourier_wavenumbers(n, period)

Return FFT-ordered angular wavenumbers for an even periodic grid.
"""
function fourier_wavenumbers(n::Int, period::Real)
    iseven(n) || throw(ArgumentError("n must be even"))
    period > 0 || throw(ArgumentError("period must be positive"))
    return (2π / period) .* vcat(0:n÷2, -n÷2+1:-1)
end

"""
    etdrk4_coefficients(k, dt; m=32)

Precompute Kassam-Trefethen ETDRK4 coefficients for the KdV equation
`q_t + 6 q q_x + q_xxx = 0`.
"""
function etdrk4_coefficients(k, dt; m::Int=32)
    m > 0 || throw(ArgumentError("m must be positive"))
    L = 1im .* k.^3
    E = exp.(dt .* L)
    E2 = exp.(dt .* L ./ 2)
    r = exp.(1im * π .* (((1:m) .- 0.5) ./ m))
    LR = dt .* L .+ transpose(r)
    Q = dt .* vec(mean((exp.(LR ./ 2) .- 1) ./ LR; dims=2))
    f1 = dt .* vec(mean((-4 .- LR .+ exp.(LR) .* (4 .- 3 .* LR .+ LR.^2)) ./ LR.^3; dims=2))
    f2 = dt .* vec(mean((2 .+ LR .+ exp.(LR) .* (-2 .+ LR)) ./ LR.^3; dims=2))
    f3 = dt .* vec(mean((-4 .- 3 .* LR .- LR.^2 .+ exp.(LR) .* (4 .- LR)) ./ LR.^3; dims=2))
    return L, E, E2, Q, f1, f2, f3
end

"""
    solve_kdv_etdrk4(q0, xgrid, times; dt, dealias=true)

Solve the periodic KdV equation `q_t + 6 q q_x + q_xxx = 0` on a uniform grid.
The first requested time must be zero, and every requested time must be an
integer multiple of `dt`.
"""
function solve_kdv_etdrk4(q0, xgrid, times; dt::Real, dealias::Bool=true)
    n = length(xgrid)
    n == length(q0) || throw(ArgumentError("q0 and xgrid must have equal length"))
    n >= 8 && iseven(n) || throw(ArgumentError("xgrid length must be even and at least 8"))
    issorted(times) || throw(ArgumentError("times must be sorted"))
    first(times) == 0 || throw(ArgumentError("the first requested time must be 0"))
    dt > 0 || throw(ArgumentError("dt must be positive"))

    dx = xgrid[2] - xgrid[1]
    period = xgrid[end] - xgrid[1] + dx
    k = fourier_wavenumbers(n, period)
    _, E, E2, Q, f1, f2, f3 = etdrk4_coefficients(k, dt)
    mask = dealias ? abs.(k) .<= (2 / 3) * maximum(abs.(k)) : trues(n)

    function nonlinear(vhat)
        u = real.(ifft(vhat))
        nhat = 3im .* k .* fft(u.^2)
        nhat .*= mask
        return nhat
    end

    v = fft(q0)
    out = Dict{Float64, Vector{Float64}}(0.0 => Float64.(q0))
    t = 0.0

    for target in times[2:end]
        target_float = Float64(target)
        while t < target_float - 0.5dt
            Nv = nonlinear(v)
            a = E2 .* v .+ Q .* Nv
            Na = nonlinear(a)
            b = E2 .* v .+ Q .* Na
            Nb = nonlinear(b)
            c = E2 .* a .+ Q .* (2 .* Nb .- Nv)
            Nc = nonlinear(c)
            v = E .* v .+ f1 .* Nv .+ 2 .* f2 .* (Na .+ Nb) .+ f3 .* Nc
            v .*= mask
            t += dt
        end
        abs(t - target_float) < 10sqrt(eps(Float64)) ||
            throw(ArgumentError("target time $target is not an integer multiple of dt=$dt"))
        out[target_float] = real.(ifft(v))
    end

    return out
end

"""
    interpolate_periodic(xgrid, values, xplot)

Linearly interpolate periodic sampled data from `xgrid` to `xplot`.
"""
function interpolate_periodic(xgrid, values, xplot)
    length(xgrid) == length(values) || throw(ArgumentError("xgrid and values must have equal length"))
    dx = xgrid[2] - xgrid[1]
    xmin = xgrid[1]
    n = length(xgrid)
    period = n * dx
    out = similar(xplot, Float64)
    for i in eachindex(xplot)
        ξ = mod(xplot[i] - xmin, period) / dx
        j0 = floor(Int, ξ)
        θ = ξ - j0
        j1 = mod(j0, n) + 1
        j2 = mod(j0 + 1, n) + 1
        out[i] = (1 - θ) * values[j1] + θ * values[j2]
    end
    return out
end

"""
    sampled_invariants(xgrid, u)

Return sampled mass and squared L2 norm using the periodic grid spacing.
"""
function sampled_invariants(xgrid, u)
    length(xgrid) == length(u) || throw(ArgumentError("xgrid and u must have equal length"))
    dx = xgrid[2] - xgrid[1]
    return (; mass=dx * sum(u), l2=dx * sum(abs2, u))
end

"""
    invariant_table(xgrid, solution)

Build a drift table for a `Dict` returned by `solve_kdv_etdrk4`.
"""
function invariant_table(xgrid, solution::AbstractDict)
    rows = []
    for t in sort(collect(keys(solution)))
        invs = sampled_invariants(xgrid, solution[t])
        push!(rows, (; time=t, invs.mass, invs.l2))
    end
    m0, e0 = rows[1].mass, rows[1].l2
    return [(; row.time,
              row.mass,
              mass_drift=abs(row.mass - m0),
              row.l2,
              l2_drift=abs(row.l2 - e0)) for row in rows]
end
