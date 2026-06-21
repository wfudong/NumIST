export crum2_step, crum2_matrices, apply_crum2_dressing

@inline _det2(a, b, c, d) = a * d - b * c

"""
    crum2_step(q, x, kappa1, theta1, kappa2, theta2; rescale=true)

Twofold Crum step using the stable branch formulas in Prins-Wahls
Eqs. (36)-(43), evaluated for the potential update Eq. (30).

In Algorithm 1, `kappa1` is the largest remaining eigenvalue and `kappa2`
is the second largest. The formula itself is order-aware through the
Wronskian and the factor `kappa2^2-kappa1^2`.
"""
function crum2_step(
    q::Real,
    x::Real,
    kappa1::Real,
    theta1::Tuple{<:Real,<:Real},
    kappa2::Real,
    theta2::Tuple{<:Real,<:Real};
    rescale::Bool=true,
)
    a1, b1 = rescale ? stable_rescale_components(theta1[1], theta1[2]) : theta1
    a2, b2 = rescale ? stable_rescale_components(theta2[1], theta2[2]) : theta2
    delta = kappa2^2 - kappa1^2

    if x <= 0
        e1 = exp(2 * kappa1 * x)
        e2 = exp(2 * kappa2 * x)
        A1 = a1 * e1 + b1
        B1 = kappa1 * (a1 * e1 - b1)
        A2 = a2 * e2 + b2
        B2 = kappa2 * (a2 * e2 - b2)
        w = _det2(A1, B1, A2, B2)

        U1 = kappa1^2 * a1 * b1 * exp(2 * kappa1 * x)
        U2 = kappa2^2 * a2 * b2 * exp(2 * kappa2 * x)
        s0 = _det2(U1, A1^2, U2, A2^2)
        return q + 8 * delta * s0 / (w * w)
    else
        e1 = exp(-2 * kappa1 * x)
        e2 = exp(-2 * kappa2 * x)
        A1 = a1 + b1 * e1
        B1 = kappa1 * (a1 - b1 * e1)
        A2 = a2 + b2 * e2
        B2 = kappa2 * (a2 - b2 * e2)
        w = _det2(A1, B1, A2, B2)

        U1 = kappa1^2 * a1 * b1 * exp(-2 * kappa1 * x)
        U2 = kappa2^2 * a2 * b2 * exp(-2 * kappa2 * x)
        s0 = _det2(U1, A1^2, U2, A2^2)
        return q + 8 * delta * s0 / (w * w)
    end
end

"""
    crum2_step(q, x, kappa1, theta1, kappa2, theta2; rescale=true)

Vectorized sampled 2-fold step.
"""
function crum2_step(
    q::AbstractVector,
    x::AbstractVector,
    kappa1::Real,
    theta1::SeedTrajectory,
    kappa2::Real,
    theta2::SeedTrajectory;
    rescale::Bool=true,
)
    length(q) == length(x) == length(theta1) == length(theta2) ||
        throw(ArgumentError("q, x, theta1, and theta2 must have equal length"))
    out = similar(q)
    @inbounds for i in eachindex(x)
        out[i] = crum2_step(q[i], x[i], kappa1, component(theta1, i), kappa2, component(theta2, i); rescale=rescale)
    end
    return out
end

function _crum2_matrices_left(x, kappa1, a1, b1, kappa2, a2, b2, z)
    e1 = exp(2 * kappa1 * x)
    e2 = exp(2 * kappa2 * x)
    A1 = a1 * e1 + b1
    B1 = kappa1 * (a1 * e1 - b1)
    A2 = a2 * e2 + b2
    B2 = kappa2 * (a2 * e2 - b2)
    w = _det2(A1, B1, A2, B2)
    delta = kappa2^2 - kappa1^2

    m(zval) = begin
        C1 = a1^2 * e1^2 - b1^2
        D1 = kappa1 * a1 * b1 * exp(2 * (zval + kappa1) * x)
        C2 = a2^2 * e2^2 - b2^2
        D2 = kappa2 * a2 * b2 * exp(2 * (zval + kappa2) * x)
        _det2(C1, D1, C2, D2)
    end

    s(zval) = begin
        U1 = kappa1^2 * a1 * b1 * exp(2 * (zval + kappa1) * x)
        U2 = kappa2^2 * a2 * b2 * exp(2 * (zval + kappa2) * x)
        _det2(U1, A1^2, U2, A2^2)
    end

    p = 0.25 * delta * A1^2 * A2^2
    m0 = m(zero(z))
    mz = m(z)
    mnegz = m(-z)
    sz = s(z)
    snegz = s(-z)

    prefm1 = 2 * kappa1 * kappa2 * delta / (w * w)
    prefm0 = 2 * delta / (w * w)
    scalar_m1 = delta * A1 * A2 / w
    half_trace = 0.5 * (kappa1^2 + kappa2^2)
    z0 = zero(prefm1 + prefm0 + scalar_m1 + p + mz + sz)

    Mminus1 = prefm1 .* [m0 mz; -mnegz -m0]
    M0 = [prefm0 * p - half_trace prefm0 * sz;
          prefm0 * snegz prefm0 * p - half_trace]
    M1 = [scalar_m1 z0; z0 -scalar_m1]
    M2 = [one(z0) z0; z0 one(z0)]
    return Mminus1, M0, M1, M2
end

function _crum2_matrices_right(x, kappa1, a1, b1, kappa2, a2, b2, z)
    e1 = exp(-2 * kappa1 * x)
    e2 = exp(-2 * kappa2 * x)
    A1 = a1 + b1 * e1
    B1 = kappa1 * (a1 - b1 * e1)
    A2 = a2 + b2 * e2
    B2 = kappa2 * (a2 - b2 * e2)
    w = _det2(A1, B1, A2, B2)
    delta = kappa2^2 - kappa1^2

    m(zval) = begin
        C1 = a1^2 - b1^2 * e1^2
        D1 = kappa1 * a1 * b1 * exp(2 * (zval - kappa1) * x)
        C2 = a2^2 - b2^2 * e2^2
        D2 = kappa2 * a2 * b2 * exp(2 * (zval - kappa2) * x)
        _det2(C1, D1, C2, D2)
    end

    s(zval) = begin
        U1 = kappa1^2 * a1 * b1 * exp(2 * (zval - kappa1) * x)
        U2 = kappa2^2 * a2 * b2 * exp(2 * (zval - kappa2) * x)
        _det2(U1, A1^2, U2, A2^2)
    end

    p = 0.25 * delta * A1^2 * A2^2
    m0 = m(zero(z))
    mz = m(z)
    mnegz = m(-z)
    sz = s(z)
    snegz = s(-z)

    prefm1 = 2 * kappa1 * kappa2 * delta / (w * w)
    prefm0 = 2 * delta / (w * w)
    scalar_m1 = delta * A1 * A2 / w
    half_trace = 0.5 * (kappa1^2 + kappa2^2)
    z0 = zero(prefm1 + prefm0 + scalar_m1 + p + mz + sz)

    Mminus1 = prefm1 .* [m0 mz; -mnegz -m0]
    M0 = [prefm0 * p - half_trace prefm0 * sz;
          prefm0 * snegz prefm0 * p - half_trace]
    M1 = [scalar_m1 z0; z0 -scalar_m1]
    M2 = [one(z0) z0; z0 one(z0)]
    return Mminus1, M0, M1, M2
end

"""
    crum2_matrices(x, kappa1, theta1, kappa2, theta2, z=0; rescale=true)

Return `(Mminus1, M0, M1, M2)` from Eqs. (36)-(43).

The argument `z` denotes `j*zeta`. The two `x <= 0` and `x > 0` branches
are deliberately implemented separately.
"""
function crum2_matrices(
    x::Real,
    kappa1::Real,
    theta1::Tuple{<:Real,<:Real},
    kappa2::Real,
    theta2::Tuple{<:Real,<:Real},
    z=zero(kappa1);
    rescale::Bool=true,
)
    a1, b1 = rescale ? stable_rescale_components(theta1[1], theta1[2]) : theta1
    a2, b2 = rescale ? stable_rescale_components(theta2[1], theta2[2]) : theta2

    if x <= 0
        return _crum2_matrices_left(x, kappa1, a1, b1, kappa2, a2, b2, z)
    else
        return _crum2_matrices_right(x, kappa1, a1, b1, kappa2, a2, b2, z)
    end
end

"""
    apply_crum2_dressing(x, kappa1, theta1, kappa2, theta2, kappa_target, theta_target)

Map one remaining seed trajectory through a 2-fold Crum step via Eq. (31).
"""
function apply_crum2_dressing(
    x::Real,
    kappa1::Real,
    theta1::Tuple{<:Real,<:Real},
    kappa2::Real,
    theta2::Tuple{<:Real,<:Real},
    kappa_target::Real,
    theta_target::Tuple{<:Real,<:Real};
    rescale::Bool=true,
)
    r = -kappa_target
    Mminus1, M0, M1, M2 = crum2_matrices(x, kappa1, theta1, kappa2, theta2, r; rescale=rescale)

    v1, v2 = theta_target
    a1, a2 = _matvec2(Mminus1, v1, v2)
    b1, b2 = _matvec2(M0, v1, v2)
    c1, c2 = _matvec2(M1, v1, v2)
    d1, d2 = _matvec2(M2, v1, v2)

    return inv(r) * a1 + b1 + r * c1 + r^2 * d1,
           inv(r) * a2 + b2 + r * c2 + r^2 * d2
end
