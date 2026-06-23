export ZeroReflection,
       zero_reflection,
       ScatteringData,
       paper_example92_scattering_data,
       evolve_scattering,
       reflection_integral,
       finite_section_coefficients,
       finite_difference_1_2,
       recover_potential,
       kdv_ist_solution

"""
    ZeroReflection()

Callable zero reflection coefficient. Use this for reflectionless scattering
data so the finite-section assembly skips quadrature work exactly.
"""
struct ZeroReflection end
(::ZeroReflection)(rho) = zero(ComplexF64)

const zero_reflection = ZeroReflection()

"""
    ScatteringData(side, tau, alpha, reflection=zero_reflection)

Scattering data for the KdV inverse-scattering step in Section 9 of `152.pdf`.

- `side = :right` represents `J+` and uses (6.4)-(6.5), then (9.2).
- `side = :left` represents `J-` and uses (6.6)-(6.7), then (9.3).
- `tau[k] > 0` corresponds to `lambda_k = -tau[k]^2`.
- `alpha[k]` is the norming constant for the chosen side.
- `reflection(rho)` evaluates `s+(rho)` or `s-(rho)` on the real axis.
"""
struct ScatteringData{F}
    side::Symbol
    tau::Vector{Float64}
    alpha::Vector{Float64}
    reflection::F
end

function ScatteringData(side::Symbol, tau, alpha, reflection=zero_reflection)
    side in (:right, :left) || throw(ArgumentError("side must be :right or :left"))
    tauv = Float64.(collect(tau))
    alphav = Float64.(collect(alpha))
    length(tauv) == length(alphav) || throw(ArgumentError("tau and alpha must have the same length"))
    all(>(0), tauv) || throw(ArgumentError("all tau values must be positive"))
    all(isfinite, tauv) || throw(ArgumentError("tau values must be finite"))
    all(isfinite, alphav) || throw(ArgumentError("alpha values must be finite"))
    return ScatteringData(side, tauv, alphav, reflection)
end

"""
    paper_example92_scattering_data(; c=pi)

Return exact reflectionless scattering data for paper Example 9.2:
`q(x) = -c/2 sech(sqrt(c) x / 2)^2`.
"""
function paper_example92_scattering_data(; c=π)
    c > 0 || throw(ArgumentError("c must be positive"))
    tau = [sqrt(c) / 2]
    alpha = [sqrt(c)]
    return (
        right=ScatteringData(:right, tau, alpha, zero_reflection),
        left=ScatteringData(:left, tau, alpha, zero_reflection),
    )
end

"""
    evolve_scattering(data, t)

Apply the KdV scattering-data evolution law from Section 4 of the paper:
`s±(rho,t)=s±(rho) exp(±8im*rho^3*t)` and
`alpha±_k(t)=alpha±_k exp(±8*tau_k^3*t)`.
"""
function evolve_scattering(data::ScatteringData, t::Real)
    sgn = data.side == :right ? 1.0 : -1.0
    evolved_alpha = data.alpha .* exp.(sgn .* 8 .* data.tau.^3 .* t)
    evolved_reflection = data.reflection isa ZeroReflection ? zero_reflection :
        (rho -> data.reflection(rho) * exp(sgn * 8im * rho^3 * t))
    return ScatteringData(data.side, data.tau, evolved_alpha, evolved_reflection)
end

# (-1)^k occurs repeatedly in the Laguerre-basis finite-section systems.
@inline _signpow(k::Integer) = isodd(k) ? -1.0 : 1.0

"""
    reflection_integral(reflection, x, m, n; Ni=10_000)

Compute the reflection integral appearing in (6.4)-(6.7) using the unit-circle
change of variables (9.4). The midpoint trapezoidal rule avoids the endpoint
`z=-1`, which corresponds to `|rho|=infinity`.
"""
function reflection_integral(reflection, x::Real, m::Int, n::Int; Ni::Int=10_000)
    n - m >= 1 || throw(ArgumentError("expected n - m >= 1"))
    Ni > 0 || throw(ArgumentError("Ni must be positive"))
    reflection isa ZeroReflection && return zero(ComplexF64)

    total = zero(ComplexF64)
    dtheta = 2π / Ni
    for j in 0:Ni-1
        theta = -π + (j + 0.5) * dtheta
        z = cis(theta)
        rho = real(im * (1 - z) / (2 * (1 + z)))
        total += reflection(rho) *
                 exp((z - 1) / (z + 1) * x) *
                 z^(m + 1) *
                 (z + 1)^(n - m - 2)
    end
    return dtheta * total
end

function _discrete_sum(data::ScatteringData, x::Real, power::Int, denom_shift::Int)
    expo_sign = data.side == :right ? -1.0 : 1.0
    value = zero(ComplexF64)
    for (tau, alpha) in zip(data.tau, data.alpha)
        # Bound-state contribution in (6.4)-(6.7), later assembled into
        # the finite systems (9.2) and (9.3).
        value += alpha * exp(expo_sign * 2 * tau * x) *
                 (0.5 - tau)^power / (0.5 + tau)^(power + denom_shift)
    end
    return value
end

function _reflection_term(data::ScatteringData, x::Real, m::Int, n::Int; Ni::Int)
    # The left Jost system is written with the reflected coordinate in the
    # paper, hence the sign flip before applying the unit-circle quadrature.
    xeff = data.side == :right ? x : -x
    return reflection_integral(data.reflection, xeff, m, n; Ni=Ni) / (2π)
end

"""
    finite_section_coefficients(data, x; Ns=5, Ni=10_000)

Solve the truncated finite-section system (9.2) or (9.3) at one spatial point.
Returns `[a0, a1, ..., a_Ns]` for `J+` data and `[b0, b1, ..., b_Ns]` for
`J-` data.
"""
function finite_section_coefficients(data::ScatteringData, x::Real; Ns::Int=5, Ni::Int=10_000)
    Ns >= 0 || throw(ArgumentError("Ns must be nonnegative"))
    dim = Ns + 1
    matrix = Matrix{ComplexF64}(I, dim, dim)
    rhs = zeros(ComplexF64, dim)

    refl_matrix = [_reflection_term(data, x, p, p + 2; Ni=Ni) for p in 0:2Ns]
    refl_rhs = [_reflection_term(data, x, m, m + 1; Ni=Ni) for m in 0:Ns]

    for m in 0:Ns
        rhs[m + 1] = _signpow(m + 1) * (_discrete_sum(data, x, m, 1) + refl_rhs[m + 1])
        for n in 0:Ns
            p = m + n
            matrix[m + 1, n + 1] += _signpow(p) * (_discrete_sum(data, x, p, 2) + refl_matrix[p + 1])
        end
    end

    return matrix \ rhs
end

"""
    finite_difference_1_2(x, y)

Return first and second derivatives on a uniform grid. The paper differentiates
a spline fit to `a0` or `b0`; this lightweight implementation uses a five-point
centered stencil in the interior and lower-order boundary stencils.
"""
function finite_difference_1_2(x::AbstractVector, y::AbstractVector)
    n = length(x)
    n == length(y) || throw(ArgumentError("x and y must have equal length"))
    n >= 5 || throw(ArgumentError("need at least five grid points"))
    h = x[2] - x[1]
    maximum(abs.((x[2:end] .- x[1:end-1]) .- h)) < 100eps(Float64) * max(1, abs(h)) ||
        throw(ArgumentError("x must be uniformly spaced"))

    dy = similar(ComplexF64.(y))
    d2y = similar(ComplexF64.(y))

    dy[1] = (-3y[1] + 4y[2] - y[3]) / (2h)
    dy[end] = (3y[end] - 4y[end-1] + y[end-2]) / (2h)
    d2y[1] = (2y[1] - 5y[2] + 4y[3] - y[4]) / h^2
    d2y[end] = (2y[end] - 5y[end-1] + 4y[end-2] - y[end-3]) / h^2

    for i in 2:n-1
        dy[i] = (y[i+1] - y[i-1]) / (2h)
        d2y[i] = (y[i+1] - 2y[i] + y[i-1]) / h^2
    end

    for i in 3:n-2
        dy[i] = (-y[i+2] + 8y[i+1] - 8y[i-1] + y[i-2]) / (12h)
        d2y[i] = (-y[i+2] + 16y[i+1] - 30y[i] + 16y[i-1] - y[i-2]) / (12h^2)
    end

    return dy, d2y
end

function _recover_one_side(data::ScatteringData, xgrid; Ns::Int, Ni::Int)
    c0 = [finite_section_coefficients(data, x; Ns=Ns, Ni=Ni)[1] for x in xgrid]
    dc0, d2c0 = finite_difference_1_2(xgrid, c0)
    # Equations (2.19) and (2.20): recover q from a0 for J+ or b0 for J-.
    q = data.side == :right ? (d2c0 .- dc0) ./ (1 .+ c0) :
                               (d2c0 .+ dc0) ./ (1 .+ c0)
    return real.(q), c0
end

function _local_stitch_gap(qright, qleft, xgrid; switch::Real=0.0, window::Real=0.25)
    idx = findall(x -> abs(x - switch) <= window, xgrid)
    isempty(idx) && return 0.0
    return maximum(abs.(qright[idx] .- qleft[idx]))
end

"""
    recover_potential(data, xgrid; Ns=5, Ni=10_000)

Recover `q(x)` from one side of scattering data via (9.2) or (9.3), followed
by (2.19) or (2.20). Returns `(q, c0)`.
"""
function recover_potential(data::ScatteringData, xgrid; Ns::Int=5, Ni::Int=10_000)
    return _recover_one_side(data, xgrid; Ns=Ns, Ni=Ni)
end

"""
    recover_potential(right_data, left_data, xgrid; Ns=5, Ni=10_000, switch=0)

Recover using both sides of the scattering data. Following Section 9.3, this
uses system (9.2), i.e. `J+`/right data, for `x < switch`, and system (9.3),
i.e. `J-`/left data, for `x >= switch`. The result is returned as a named
tuple containing the combined potential, both branch reconstructions, and
`stitch_gap`, a diagnostic measuring the local branch mismatch near `switch`.
"""
function recover_potential(right_data::ScatteringData, left_data::ScatteringData, xgrid;
                           Ns::Int=5, Ni::Int=10_000, switch::Real=0.0)
    right_data.side == :right || throw(ArgumentError("right_data must have side=:right"))
    left_data.side == :left || throw(ArgumentError("left_data must have side=:left"))

    qright, a0 = _recover_one_side(right_data, xgrid; Ns=Ns, Ni=Ni)
    qleft, b0 = _recover_one_side(left_data, xgrid; Ns=Ns, Ni=Ni)
    q = similar(qright)
    for i in eachindex(xgrid)
        q[i] = xgrid[i] < switch ? qright[i] : qleft[i]
    end

    return (; q, qright, qleft, a0, b0,
            stitch_gap=_local_stitch_gap(qright, qleft, xgrid; switch=switch))
end

"""
    kdv_ist_solution(right0, left0, xgrid, t; Ns=5, Ni=10_000, switch=0)

Compute `u(x,t)` for the paper's KdV convention `u_t - 6u*u_x + u_xxx = 0`
by evolving scattering data as in Section 4 and solving the finite-section
inverse problem as in Section 9.2.
"""
function kdv_ist_solution(right0::ScatteringData, left0::ScatteringData, xgrid, t;
                          Ns::Int=5, Ni::Int=10_000, switch::Real=0.0)
    return recover_potential(
        evolve_scattering(right0, t),
        evolve_scattering(left0, t),
        xgrid;
        Ns=Ns,
        Ni=Ni,
        switch=switch,
    )
end
