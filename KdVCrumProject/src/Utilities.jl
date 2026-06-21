export sech_stable,
       exact_one_soliton,
       evolve_beta,
       shifted_norming_constants,
       decimal_digits_to_bits,
       trapezoid,
       relative_l2_error,
       conserved_errors,
       alternating_norming_constants,
       uniform_kappa,
       soliton_gas_kappa,
       clip_reflectionless!

"""
    sech_stable(x)

Stable real hyperbolic secant.
"""
function sech_stable(x::Real)
    e = exp(-abs(x))
    return 2 * e / (1 + e^2)
end

"""
    exact_one_soliton(x, kappa; beta=1)

Exact one-soliton reflectionless KdV potential generated from `q0=0`:

`q(x) = 2*kappa^2*sech(kappa*x - 0.5*log(beta))^2`.
"""
function exact_one_soliton(x::AbstractVector, kappa::Real; beta::Real=1, shift::Real=0)
    beta > 0 || throw(ArgumentError("the one-soliton beta must be positive"))
    beta_eff = beta * exp(2 * kappa * shift)
    center_phase = 0.5 * log(beta_eff)
    return [2 * kappa^2 * sech_stable(kappa * xi - center_phase)^2 for xi in x]
end

"""
    shifted_norming_constants(kappa, beta, shifts; precision=nothing, digits=nothing)

Apply per-eigenvalue spatial/phase shifts to target norming constants.

With the conventions in this project, a one-soliton with norming constant
`beta > 0` is centered at `log(beta)/(2*kappa)`. Therefore a requested shift
`s_n` is encoded by

`beta_n -> beta_n * exp(2*kappa_n*s_n)`.

The sign of each `beta_n` is preserved. Use `digits` or `precision` to compute
the exponentials in `BigFloat` arithmetic.
"""
function shifted_norming_constants(kappa::AbstractVector, beta::AbstractVector, shifts::AbstractVector; precision=nothing, digits=nothing)
    length(kappa) == length(beta) == length(shifts) ||
        throw(ArgumentError("kappa, beta, and shifts must have equal length"))

    bits = _precision_bits(precision, digits)
    if bits === nothing
        return beta .* exp.(2 .* kappa .* shifts)
    end

    return setprecision(BigFloat, bits) do
        kb = BigFloat.(kappa)
        bb = BigFloat.(beta)
        sb = BigFloat.(shifts)
        bb .* exp.(2 .* kb .* sb)
    end
end

"""
    evolve_beta(kappa, beta0, t)

Time-evolve target norming constants for the KdV convention used in the
paper, `q_t + 6 q q_x + q_xxx = 0`.

The eigenvalues `kappa` are fixed and
`beta_n(t) = beta_n(0) * exp(8*kappa_n^3*t)`.
"""
function _evolve_beta_impl(kappa::AbstractVector, beta0::AbstractVector, t::Real)
    length(kappa) == length(beta0) || throw(ArgumentError("kappa and beta0 must have equal length"))
    return beta0 .* exp.(8 .* (kappa .^ 3) .* t)
end

function evolve_beta(kappa::AbstractVector, beta0::AbstractVector, t::Real; precision=nothing, digits=nothing)
    bits = _precision_bits(precision, digits)
    bits === nothing && return _evolve_beta_impl(kappa, beta0, t)

    return setprecision(BigFloat, bits) do
        kb = BigFloat.(kappa)
        bb = BigFloat.(beta0)
        tb = BigFloat(t)
        _evolve_beta_impl(kb, bb, tb)
    end
end

"""
    decimal_digits_to_bits(digits)

Convert a requested decimal digit count to a BigFloat precision in bits, with
a small guard margin for intermediate arithmetic.
"""
decimal_digits_to_bits(digits::Integer) = ceil(Int, digits * log2(10)) + 16

function _precision_bits(precision, digits)
    precision === nothing || digits === nothing ||
        throw(ArgumentError("use either precision=<bits> or digits=<decimal digits>, not both"))

    if digits !== nothing
        digits isa Integer && digits > 0 || throw(ArgumentError("digits must be a positive integer"))
        return decimal_digits_to_bits(digits)
    elseif precision !== nothing
        precision isa Integer && precision > 0 || throw(ArgumentError("precision must be a positive integer bit count"))
        return precision
    else
        return nothing
    end
end

"""
    trapezoid(x, y)

Composite trapezoidal integration on a sampled grid.
"""
function trapezoid(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || throw(ArgumentError("x and y must have equal length"))
    length(x) >= 2 || return zero(eltype(y))

    s = zero(promote_type(eltype(x), eltype(y), Float64))
    @inbounds for i in 1:(length(x)-1)
        s += 0.5 * (x[i + 1] - x[i]) * (y[i + 1] + y[i])
    end
    return s
end

"""
    relative_l2_error(y, yref)

Relative discrete 2-norm error.
"""
function relative_l2_error(y::AbstractVector, yref::AbstractVector)
    length(y) == length(yref) || throw(ArgumentError("vectors must have equal length"))
    return norm(y .- yref) / norm(yref)
end

"""
    conserved_errors(x, q, kappa; q0=zeros)

Return `(E1, E2)` from Prins-Wahls Eq. (48), useful when the exact target
potential is unknown.
"""
function conserved_errors(x::AbstractVector, q::AbstractVector, kappa::AbstractVector; q0=nothing)
    q0v = q0 === nothing ? zeros(eltype(q), length(q)) : q0
    length(q0v) == length(q) || throw(ArgumentError("q0 and q must have equal length"))

    i1 = trapezoid(x, q .- q0v)
    i2 = trapezoid(x, q .^ 2 .- q0v .^ 2)
    denom1 = 4 * sum(kappa)
    denom2 = 16 * sum(kappa .^ 3)

    E1 = abs(i1 / denom1 - 1)
    E2 = abs(3 * i2 / denom2 - 1)
    return E1, E2
end

"""
    alternating_norming_constants(N)

Norming constants `beta_n=(-1)^(N-n)`, a useful regular family for symmetric
reflectionless examples.
"""
function alternating_norming_constants(N::Integer)
    N >= 0 || throw(ArgumentError("N must be nonnegative"))
    return [isodd(N - n) ? -1.0 : 1.0 for n in 1:N]
end

"""
    uniform_kappa(N; start=0.5, spacing=0.35)

Generate a strictly increasing positive `kappa` vector.
"""
function uniform_kappa(N::Integer; start::Real=0.5, spacing::Real=0.35)
    N >= 0 || throw(ArgumentError("N must be nonnegative"))
    start > 0 || throw(ArgumentError("start must be positive"))
    spacing > 0 || throw(ArgumentError("spacing must be positive"))
    return [float(start + (n - 1) * spacing) for n in 1:N]
end

"""
    soliton_gas_kappa(N, lambda1)

Generate the ordered soliton-gas eigenvalue magnitudes

`kappa_j = lambda1/N * sqrt(2*j*N - j^2)`, for `j = 1, ..., N`.

The output is positive, strictly increasing, and has `kappa_N = lambda1`.
"""
function soliton_gas_kappa(N::Integer, lambda1::Real)
    N > 0 || throw(ArgumentError("N must be positive"))
    lambda1 > 0 || throw(ArgumentError("lambda1 must be positive"))
    isfinite(lambda1) || throw(ArgumentError("lambda1 must be finite"))

    scale = float(lambda1) / N
    T = typeof(scale)
    return [scale * sqrt(T(2 * j * N - j^2)) for j in 1:N]
end

"""
    clip_reflectionless!(q; rtol=1000*eps())

For reflectionless potentials from a zero background, remove tiny negative
roundoff excursions without hiding large instability.
"""
function clip_reflectionless!(q::AbstractVector{T}; rtol=1000 * eps(T)) where {T<:AbstractFloat}
    scale = max(one(T), maximum(abs, q))
    tol = rtol * scale
    @inbounds for i in eachindex(q)
        if q[i] < 0 && abs(q[i]) <= tol
            q[i] = zero(T)
        end
    end
    return q
end
