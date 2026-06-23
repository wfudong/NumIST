include(joinpath(@__DIR__, "..", "src", "NumISTNext.jl"))
using .NumISTNext
using Plots
using Printf

# This script still contains the old PDE time-step diagnostics for Examples 9.1
# and 9.3. The paper's Section 9 finite-section inverse-scattering algorithm is
# implemented in `src/ist/finite_section.jl`; use `reproduce_example92_ist.jl`
# for the currently migrated paper-algorithm reproduction.

gr()

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

function figure1(; outdir)
    xmin, xmax = -40.0, 40.0
    n = 4096
    xgrid = collect(range(xmin, xmax; length=n + 1))[1:end-1]
    q0 = paper_example91_initial.(xgrid)
    times = section93_time_slices().example91
    sol = solve_kdv_etdrk4(q0, xgrid, times; dt=2e-4)
    xplot = collect(range(-5.0, 7.0; length=900))

    plt = plot(layout=(5, 1), size=(560, 650), margin=3Plots.mm)
    labels = ["t=0", "t=0.28", "t=0.4", "t=0.5", "t=0.6"]
    for (i, t) in enumerate(times)
        y = interpolate_periodic(xgrid, sol[t], xplot)
        paper_subplot!(plt, xplot, y, i; xlims=(-5, 7), ylims=(-1, 1),
                       labeltext=labels[i], initial_text=(i == 1))
    end
    save_both(plt, joinpath(outdir, "figure1_example91"))
    return (; xgrid, sol)
end

function figure2(; outdir)
    times = section93_time_slices().example92
    xplot = collect(range(-5.0, 7.0; length=900))
    plt = plot(layout=(3, 1), size=(560, 650), margin=3Plots.mm)
    labels = ["t=0", "t=0.5", "t=1"]
    for (i, t) in enumerate(times)
        y = paper_example92_exact.(xplot, t)
        paper_subplot!(plt, xplot, y, i; xlims=(-5, 7), ylims=(-2, 0),
                       labeltext=labels[i], initial_text=(i == 1))
    end
    save_both(plt, joinpath(outdir, "figure2_example92"))
    return nothing
end

function figure3(; outdir)
    xplot = collect(range(-8.0, 8.0; length=1200))
    y = paper_example93_initial.(xplot)
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

function figure4(; outdir)
    xmin, xmax = -48.0, 48.0
    n = 4096
    xgrid = collect(range(xmin, xmax; length=n + 1))[1:end-1]
    q0 = paper_example93_initial.(xgrid)
    times = section93_time_slices().example93
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

function write_validation_summary(outdir, ex1, ex4)
    summary_path = joinpath(outdir, "validation_summary.txt")
    open(summary_path, "w") do io
        println(io, "Validation summary for 152.pdf numerical examples")
        println(io)
        println(io, "Figure 1 / Example 9.1")
        println(io, "  Initial data: q(x) = x exp(-x^2)")
        for row in invariant_table(ex1.xgrid, ex1.sol)
            @printf(io, "  t=%0.5f mass=% .12e mass_drift=% .3e L2=% .12e L2_drift=% .3e\n",
                    row.time, row.mass, row.mass_drift, row.l2, row.l2_drift)
        end
        println(io)
        println(io, "Figure 2 / Example 9.2")
        println(io, "  Initial data: q(x) = -c/2 sech^2(sqrt(c)x/2), c = pi")
        @printf(io, "  Exact lambda1 = -pi/4 = %.14f\n", -π / 4)
        @printf(io, "  Exact alpha+ = alpha- = sqrt(pi) = %.14f\n", sqrt(π))
        println(io)
        println(io, "Figure 3 / Example 9.3")
        println(io, "  Initial data: exp(x) cos(4x) for x < 0, exp(-x) J0(2x) for x > 0")
        println(io)
        println(io, "Figure 4 / Example 9.3")
        for row in invariant_table(ex4.xgrid, ex4.sol)
            @printf(io, "  t=%0.5f mass=% .12e mass_drift=% .3e L2=% .12e L2_drift=% .3e\n",
                    row.time, row.mass, row.mass_drift, row.l2, row.l2_drift)
        end
    end
    return summary_path
end

function run_all_examples(; outdir=joinpath(@__DIR__, "..", "output", "figures"))
    mkpath(outdir)
    ex1 = figure1(outdir=outdir)
    figure2(outdir=outdir)
    figure3(outdir=outdir)
    ex4 = figure4(outdir=outdir)
    summary_path = write_validation_summary(outdir, ex1, ex4)
    return (; ex1, ex4, summary_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    results = run_all_examples()
    println("Wrote ", results.summary_path)
end
