export plot_potential, plot_comparison

using Plots
using LaTeXStrings

"""
    plot_potential(x, q; kwargs...)

Plot a sampled KdV potential. Numerical routines do not depend on plotting.
"""
function plot_potential(x, q; label=L"q(x)", kwargs...)
    return plot(
        x,
        q;
        label=label,
        xlabel=L"x",
        ylabel=L"q(x)",
        linewidth=2.5,
        framestyle=:box,
        kwargs...,
    )
end

"""
    plot_comparison(x, q, qref; kwargs...)

Plot computed and reference potentials on the same axes.
"""
function plot_comparison(x, q, qref; label=L"Crum", reflabel=L"exact", kwargs...)
    plt = plot_potential(x, q; label=label, kwargs...)
    plot!(plt, x, qref; label=reflabel, linewidth=2, linestyle=:dash)
    return plt
end
