export crum1_step, crum1_matrices, apply_crum1_dressing

const _Pair = Tuple{T,T} where {T<:Real}

"""
    crum1_step(q, x, kappa, theta; rescale=true)

One Darboux/Crum step using the stable branch formulas in Prins-Wahls
Eqs. (32)-(35), evaluated for the potential update Eq. (30).

The two exponential branches are intentionally kept separate.
"""
function crum1_step(q::Real, x::Real, kappa::Real, theta::Tuple{<:Real,<:Real}; rescale::Bool=true)
    theta1, theta2 = rescale ? stable_rescale_components(theta[1], theta[2]) : theta

    if x <= 0
        e = exp(2 * kappa * x)
        w = theta1 * e + theta2
        return -q + 8 * kappa^2 * theta1 * theta2 * e / (w * w)
    else
        e = exp(-2 * kappa * x)
        w = theta1 + theta2 * e
        return -q + 8 * kappa^2 * theta1 * theta2 * e / (w * w)
    end
end

"""
    crum1_step(q, x, kappa, theta_seed; rescale=true)

Vectorized sampled 1-fold step.
"""
function crum1_step(q::AbstractVector, x::AbstractVector, kappa::Real, theta::SeedTrajectory; rescale::Bool=true)
    length(q) == length(x) == length(theta) || throw(ArgumentError("q, x, and theta must have equal length"))
    out = similar(q)
    @inbounds for i in eachindex(x)
        out[i] = crum1_step(q[i], x[i], kappa, component(theta, i); rescale=rescale)
    end
    return out
end

"""
    crum1_matrices(x, kappa, theta, z=0; rescale=true)

Return `(Mminus1, M0, M1)` from Eqs. (32)-(34).

The argument `z` denotes `j*zeta`. In this project the update of a remaining
bound-state trajectory uses `z=-kappa_target`, while the potential update uses
`z=0`.
"""
function crum1_matrices(
    x::Real,
    kappa::Real,
    theta::Tuple{<:Real,<:Real},
    z=zero(kappa);
    rescale::Bool=true,
)
    theta1, theta2 = rescale ? stable_rescale_components(theta[1], theta[2]) : theta

    if x <= 0
        e = exp(2 * kappa * x)
        w = theta1 * e + theta2
        pref = 2 * kappa^2 * theta1 * theta2 / (w * w)
        m0 = kappa * (theta2 - theta1 * e) / w
        e12 = exp(2 * (z + kappa) * x)
        e21 = exp(-2 * (z - kappa) * x)
        z0 = zero(pref + m0 + e12 + e21)

        Mminus1 = [pref * e pref * e12; -pref * e21 -pref * e]
        M0 = [m0 z0; z0 m0]
        M1 = [-one(z0) z0; z0 one(z0)]
        return Mminus1, M0, M1
    else
        e = exp(-2 * kappa * x)
        w = theta1 + theta2 * e
        pref = 2 * kappa^2 * theta1 * theta2 / (w * w)
        m0 = kappa * (theta2 * e - theta1) / w
        e12 = exp(2 * (z - kappa) * x)
        e21 = exp(-2 * (z + kappa) * x)
        z0 = zero(pref + m0 + e12 + e21)

        Mminus1 = [pref * e pref * e12; -pref * e21 -pref * e]
        M0 = [m0 z0; z0 m0]
        M1 = [-one(z0) z0; z0 one(z0)]
        return Mminus1, M0, M1
    end
end

@inline function _matvec2(M, v1, v2)
    return M[1, 1] * v1 + M[1, 2] * v2,
           M[2, 1] * v1 + M[2, 2] * v2
end

"""
    apply_crum1_dressing(x, kappa_added, theta_added, kappa_target, theta_target)

Map one remaining seed trajectory through a 1-fold Crum step via Eq. (31).
For a bound state at `zeta=j*kappa_target`, `j*zeta=-kappa_target`.
"""
function apply_crum1_dressing(
    x::Real,
    kappa_added::Real,
    theta_added::Tuple{<:Real,<:Real},
    kappa_target::Real,
    theta_target::Tuple{<:Real,<:Real};
    rescale::Bool=true,
)
    z = -kappa_target
    Mminus1, M0, M1 = crum1_matrices(x, kappa_added, theta_added, z; rescale=rescale)
    r = -kappa_target

    v1, v2 = theta_target
    a1, a2 = _matvec2(Mminus1, v1, v2)
    b1, b2 = _matvec2(M0, v1, v2)
    c1, c2 = _matvec2(M1, v1, v2)

    return inv(r) * a1 + b1 + r * c1,
           inv(r) * a2 + b2 + r * c2
end
