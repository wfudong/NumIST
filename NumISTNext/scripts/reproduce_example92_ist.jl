include(joinpath(@__DIR__, "..", "src", "NumISTNext.jl"))
using .NumISTNext
using Plots
using Printf

gr()

function save_both(plt, basename)
    mkpath(dirname(basename))
    savefig(plt, basename * ".png")
    savefig(plt, basename * ".pdf")
end

function figure_example92_ist(; outdir=joinpath(@__DIR__, "..", "output", "figures_ist"))
    data = paper_example92_scattering_data(c=π)
    x = collect(range(-5.0, 7.0; length=900))
    times = section93_time_slices().example92

    plt = plot(layout=(3, 1), size=(560, 650), margin=3Plots.mm)
    max_errors = Float64[]
    for (i, t) in enumerate(times)
        result = kdv_ist_solution(data.right, data.left, x, t; Ns=8, Ni=256)
        exact = paper_example92_exact.(x, t)
        push!(max_errors, maximum(abs.(result.q[5:end-4] .- exact[5:end-4])))

        plot!(plt[i], x, result.q;
              color=:black,
              linewidth=1.6,
              legend=false,
              xlims=(-5, 7),
              ylims=(-2, 0),
              framestyle=:box,
              grid=false,
              xlabel="x",
              ylabel="",
              xticks=-4:2:6,
              yticks=-2:1:0,
              tickfont=font(9),
              guidefont=font(10))
        annotate!(plt[i], -3.4, -0.3, text("t=$t", 10, :black))
        i == 1 && annotate!(plt[i], 0.6, -0.32, text("Initial Data", 10, :black))
    end

    save_both(plt, joinpath(outdir, "figure2_example92_finite_section_ist"))

    summary_path = joinpath(outdir, "validation_summary_example92_ist.txt")
    open(summary_path, "w") do io
        println(io, "Finite-section IST validation for 152.pdf Example 9.2")
        println(io, "Algorithm: evolve scattering data by Section 4, solve (9.2)/(9.3), recover q by (2.19)/(2.20).")
        for (t, err) in zip(times, max_errors)
            @printf(io, "  t=%0.3f max_abs_error=% .6e\n", t, err)
        end
    end

    return (; figure=joinpath(outdir, "figure2_example92_finite_section_ist"), summary_path, max_errors)
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = figure_example92_ist()
    println("Wrote ", result.figure, ".png/.pdf")
    println("Wrote ", result.summary_path)
end
