using FFTW
using LinearAlgebra
using Printf
using Statistics

try
    using Plots
catch err
    @error "Plots.jl is required for saving figures in this script" exception=(err, catch_backtrace())
    rethrow()
end

gr()

sech(x) = inv(cosh(x))

function besselj0_quad(x; ntheta::Int=512)
    h = π / ntheta
    s = 0.5 * (cos(x * sin(0.0)) + cos(x * sin(π)))
    for j in 1:ntheta-1
        θ = j * h
        s += cos(x * sin(θ))
    end
    return h * s / π
end

function fourier_wavenumbers(n::Int, period::Real)
    @assert iseven(n)
    return (2π / period) .* vcat(0:n÷2, -n÷2+1:-1)
end

function etdrk4_coefficients(k, dt; m::Int=32)
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

function solve_kdv_etdrk4(q0, xgrid, times; dt::Real, dealias::Bool=true)
    n = length(xgrid)
    @assert n >= 8 && iseven(n)
    @assert issorted(times)
    @assert first(times) == 0

    period = xgrid[end] - xgrid[1] + (xgrid[2] - xgrid[1])
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
    out = Dict{Float64, Vector{Float64}}(0.0 => copy(q0))
    t = 0.0

    for target in times[2:end]
        while t < target - 0.5dt
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
        @assert abs(t - target) < 10sqrt(eps(Float64)) "target time is not an integer multiple of dt"
        out[Float64(target)] = real.(ifft(v))
    end

    return out
end

function interpolate_periodic(xgrid, values, xplot)
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

function paper_subplot!(plt, x, y, row; xlims, ylims, labeltext, initial_text=false)
    plot!(plt[row], x, y;
          color=:black,
          linewidth=1.6,
          legend=false,
          xlims=xlims,
          ylims=ylims,
          framestyle=:box,
          grid=false,
          xlabel="x",
          ylabel="",
          xticks=collect(ceil(Int, xlims[1] / 2) * 2:2:floor(Int, xlims[2] / 2) * 2),
          yticks=collect(ceil(Int, ylims[1]):floor(Int, ylims[2])),
          tickfont=font(9),
          guidefont=font(10))
    annotate!(plt[row], xlims[1] + 0.13 * (xlims[2] - xlims[1]),
              ylims[2] - 0.15 * (ylims[2] - ylims[1]),
              text(labeltext, 10, :black))
    if initial_text
        annotate!(plt[row], xlims[1] + 0.48 * (xlims[2] - xlims[1]),
                  ylims[2] - 0.16 * (ylims[2] - ylims[1]),
                  text("Initial Data", 10, :black))
    end
    return plt
end

function save_both(plt, basename)
    mkpath(dirname(basename))
    savefig(plt, basename * ".png")
    savefig(plt, basename * ".pdf")
end

function figure1(; outdir="figures")
    xmin, xmax = -40.0, 40.0
    n = 4096
    xgrid = collect(range(xmin, xmax; length=n + 1))[1:end-1]
    q0 = xgrid .* exp.(-xgrid.^2)
    times = [0.0, 0.08, 0.16]
    sol = solve_kdv_etdrk4(q0, xgrid, times; dt=2e-4)
    xplot = collect(range(-5.0, 7.0; length=900))

    plt = plot(layout=(3, 1), size=(560, 650), margin=3Plots.mm)
    labels = ["t=0", "t=0.08", "t=0.16"]
    for (i, t) in enumerate(times)
        y = interpolate_periodic(xgrid, sol[t], xplot)
        paper_subplot!(plt, xplot, y, i; xlims=(-5, 7), ylims=(-1, 1),
                       labeltext=labels[i], initial_text=(i == 1))
    end
    save_both(plt, joinpath(outdir, "figure1_example91"))
    return (; xgrid, sol)
end

function figure2(; outdir="figures")
    c = π
    times = [0.0, 0.5, 1.0]
    xplot = collect(range(-5.0, 7.0; length=900))
    plt = plot(layout=(3, 1), size=(560, 650), margin=3Plots.mm)
    labels = ["t=0", "t=0.5", "t=1"]
    for (i, t) in enumerate(times)
        y = @. -c / 2 * sech(sqrt(c) * (xplot - c * t) / 2)^2
        paper_subplot!(plt, xplot, y, i; xlims=(-5, 7), ylims=(-2, 0),
                       labeltext=labels[i], initial_text=(i == 1))
    end
    save_both(plt, joinpath(outdir, "figure2_example92"))
    return nothing
end

function example93_initial(x)
    if x < 0
        return exp(x) * cos(4x)
    end
    return exp(-x) * besselj0_quad(2x)
end

function figure3(; outdir="figures")
    xplot = collect(range(-8.0, 8.0; length=1200))
    y = example93_initial.(xplot)
    plt = plot(xplot, y;
               color=:black,
               linewidth=1.6,
               legend=false,
               xlims=(-8, 8),
               ylims=(-0.5, 1.05),
               xlabel="",
               ylabel="q(x)",
               framestyle=:box,
               grid=false,
               size=(560, 280),
               xticks=-8:2:8,
               yticks=[-0.5, 0, 0.5, 1.0],
               tickfont=font(9),
               guidefont=font(10),
               margin=3Plots.mm)
    save_both(plt, joinpath(outdir, "figure3_example93_initial"))
    return y
end

function figure4(; outdir="figures")
    xmin, xmax = -48.0, 48.0
    n = 4096
    xgrid = collect(range(xmin, xmax; length=n + 1))[1:end-1]
    q0 = example93_initial.(xgrid)
    times = [0.0, 0.015, 0.03]
    sol = solve_kdv_etdrk4(q0, xgrid, times; dt=5e-5)
    xplot = collect(range(-7.0, 7.0; length=900))

    plt = plot(layout=(3, 1), size=(560, 650), margin=3Plots.mm)
    labels = ["t=0", "t=0.015", "t=0.03"]
    for (i, t) in enumerate(times)
        y = interpolate_periodic(xgrid, sol[t], xplot)
        paper_subplot!(plt, xplot, y, i; xlims=(-7, 7), ylims=(-1, 1),
                       labeltext=labels[i], initial_text=(i == 1))
    end
    save_both(plt, joinpath(outdir, "figure4_example93"))
    return (; xgrid, sol)
end

function invariant_table(example)
    xgrid = example.xgrid
    sol = example.sol
    dx = xgrid[2] - xgrid[1]
    rows = []
    for t in sort(collect(keys(sol)))
        u = sol[t]
        push!(rows, (time=t, mass=dx * sum(u), l2=dx * sum(abs2, u)))
    end
    m0, e0 = rows[1].mass, rows[1].l2
    return [(; row.time,
              row.mass,
              mass_drift=abs(row.mass - m0),
              row.l2,
              l2_drift=abs(row.l2 - e0)) for row in rows]
end

function write_validation_summary(outdir, ex1, ex4)
    summary_path = joinpath(outdir, "validation_summary.txt")
    open(summary_path, "w") do io
        println(io, "Validation summary for 152.pdf numerical examples")
        println(io)
        println(io, "Figure 1 / Example 9.1")
        println(io, "  Initial data: q(x) = x exp(-x^2)")
        println(io, "  Paper time slices: t = 0, 0.08, 0.16")
        println(io, "  Paper reports lambda1 approx -0.0138384593995, alpha- approx 0.2055954681199, alpha+ approx 0.0416040800785.")
        for row in invariant_table(ex1)
            @printf(io, "  t=%0.5f mass=% .12e mass_drift=% .3e L2=% .12e L2_drift=% .3e\n",
                    row.time, row.mass, row.mass_drift, row.l2, row.l2_drift)
        end
        println(io)
        println(io, "Figure 2 / Example 9.2")
        println(io, "  Initial data: q(x) = -c/2 sech^2(sqrt(c)x/2), c = pi")
        @printf(io, "  Exact lambda1 = -pi/4 = %.14f\n", -π / 4)
        @printf(io, "  Exact alpha+ = alpha- = sqrt(pi) = %.14f\n", sqrt(π))
        println(io, "  Paper time slices: t = 0, 0.5, 1")
        println(io, "  Figure generated from the exact solitary-wave formula reported in the paper.")
        println(io)
        println(io, "Figure 3 / Example 9.3")
        println(io, "  Initial data: exp(x) cos(4x) for x < 0, exp(-x) J0(2x) for x > 0")
        println(io, "  J0 is evaluated by its integral representation J0(x) = pi^-1 int_0^pi cos(x sin(theta)) dtheta.")
        println(io)
        println(io, "Figure 4 / Example 9.3")
        println(io, "  Paper time slices: t = 0, 0.015, 0.03")
        for row in invariant_table(ex4)
            @printf(io, "  t=%0.5f mass=% .12e mass_drift=% .3e L2=% .12e L2_drift=% .3e\n",
                    row.time, row.mass, row.mass_drift, row.l2, row.l2_drift)
        end
    end
    return summary_path
end

function run_all_examples(; outdir="figures")
    mkpath(outdir)
    ex1 = figure1(outdir=outdir)
    figure2(outdir=outdir)
    figure3(outdir=outdir)
    ex4 = figure4(outdir=outdir)
    summary_path = write_validation_summary(outdir, ex1, ex4)
    return (; ex1, ex4, summary_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_all_examples()
end
