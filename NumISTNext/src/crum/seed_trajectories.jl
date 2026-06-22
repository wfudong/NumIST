export component,
       setcomponent!,
       initial_seed_trajectories,
       stable_scale_exponent,
       stable_rescale_components,
       stable_rescale_seed_at!

"""
    component(seed, i)

Return `(theta1, theta2)` at grid index `i`.
"""
component(seed::SeedTrajectory, i::Integer) = (seed.c1[i], seed.c2[i])

"""
    setcomponent!(seed, i, theta1, theta2)

Set the two VE components at grid index `i`.
"""
function setcomponent!(seed::SeedTrajectory, i::Integer, theta1, theta2)
    seed.c1[i] = theta1
    seed.c2[i] = theta2
    return seed
end

"""
    initial_seed_trajectories(problem)

Construct the VE seed trajectories for the zero background potential.

For `q0(x)=0`, the VE Jost solutions at `j*kappa` are
`phi_E=[1,0]` and `psi_E=[0,1]`, so

`theta_E = phi_E - (-1)^N beta psi_E`.
"""
function initial_seed_trajectories(problem::CrumProblem{T}) where {T}
    N = length(problem.kappa)
    parity = isodd(N) ? -one(T) : one(T)
    beta = effective_beta(problem)
    nx = length(problem.x)

    seeds = Vector{SeedTrajectory{T}}(undef, N)
    for n in 1:N
        c1 = ones(T, nx)
        c2 = fill(-parity * beta[n], nx)
        seeds[n] = SeedTrajectory{T}(c1, c2)
    end
    return seeds
end

"""
    stable_scale_exponent(theta1, theta2)

Return the integer exponent `e` used in Prins-Wahls Eq. (27), so that
`theta <- 2^(-e) theta`.

If one component is exactly zero, use the nonzero component alone. If both
components vanish, return zero and leave the trajectory unchanged.
"""
function stable_scale_exponent(theta1::Real, theta2::Real)
    a = abs(theta1)
    b = abs(theta2)

    if isfinite(a) && isfinite(b) && a > 0 && b > 0
        return round(Int, 0.5 * log2(a) + 0.5 * log2(b))
    elseif isfinite(a) && a > 0
        return round(Int, log2(a))
    elseif isfinite(b) && b > 0
        return round(Int, log2(b))
    else
        return 0
    end
end

"""
    stable_rescale_components(theta1, theta2)

Apply the power-of-two rescaling from Eq. (27). Power-of-two scaling avoids
unnecessary rounding in IEEE floating point arithmetic.
"""
function stable_rescale_components(theta1::T, theta2::T) where {T<:Real}
    e = stable_scale_exponent(theta1, theta2)
    return ldexp(theta1, -e), ldexp(theta2, -e)
end

function stable_rescale_components(theta::Tuple{T,T}) where {T<:Real}
    return stable_rescale_components(theta[1], theta[2])
end

"""
    stable_rescale_seed_at!(seed, i)

Rescale one sampled VE trajectory value in place.
"""
function stable_rescale_seed_at!(seed::SeedTrajectory, i::Integer)
    theta1, theta2 = stable_rescale_components(component(seed, i))
    return setcomponent!(seed, i, theta1, theta2)
end
