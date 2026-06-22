export CrumProblem,
       SeedTrajectory,
       CrumStepRecord,
       CrumTransformResult,
       effective_beta,
       validate_spectrum,
       has_zero_background

"""
    SeedTrajectory(c1, c2)

Vector-valued VE seed trajectory sampled on the `x` grid. The fields `c1`
and `c2` hold the first and second VE components.
"""
struct SeedTrajectory{T<:Real}
    c1::Vector{T}
    c2::Vector{T}
end

function SeedTrajectory(c1::AbstractVector, c2::AbstractVector)
    length(c1) == length(c2) || throw(ArgumentError("seed components must have the same length"))
    T = promote_type(Float64, eltype(c1), eltype(c2))
    return SeedTrajectory{T}(T.(collect(c1)), T.(collect(c2)))
end

Base.length(seed::SeedTrajectory) = length(seed.c1)
Base.copy(seed::SeedTrajectory) = SeedTrajectory(copy(seed.c1), copy(seed.c2))

"""
    CrumProblem(x, kappa, beta; q0=nothing)

Complete sampled Crum-transform input:

- `x`: spatial grid.
- `kappa`: positive eigenvalue magnitudes, strictly increasing.
- `beta`: real nonzero target norming constants.
- `shifts`: optional per-eigenvalue phase/spatial shifts. A shift `s_n`
  multiplies `beta_n` by `exp(2*kappa_n*s_n)`.
- `q0`: sampled background potential. Defaults to zero.
"""
struct CrumProblem{T<:AbstractFloat}
    x::Vector{T}
    kappa::Vector{T}
    beta::Vector{T}
    shifts::Vector{T}
    q0::Vector{T}
end

function CrumProblem(x, kappa, beta; shifts=nothing, q0=nothing, number_type=nothing)
    xraw = collect(x)
    kraw = collect(kappa)
    braw = collect(beta)
    sraw = shifts === nothing ? nothing : collect(shifts)
    qraw = q0 === nothing ? nothing : collect(q0)

    T = if number_type === nothing
        promote_type(
            Float64,
            eltype(xraw),
            eltype(kraw),
            eltype(braw),
            sraw === nothing ? Float64 : eltype(sraw),
            qraw === nothing ? Float64 : eltype(qraw),
        )
    else
        number_type
    end
    T <: AbstractFloat || throw(ArgumentError("number_type must be an AbstractFloat type"))

    xv = T.(xraw)
    kv = T.(kraw)
    bv = T.(braw)
    sv = sraw === nothing ? zeros(T, length(kv)) : T.(sraw)
    qv = qraw === nothing ? zeros(T, length(xv)) : T.(qraw)

    length(xv) == length(qv) || throw(ArgumentError("x and q0 must have the same length"))
    length(kv) == length(sv) || throw(ArgumentError("kappa and shifts must have equal length"))
    all(isfinite, xv) || throw(ArgumentError("x grid must be finite"))
    all(isfinite, sv) || throw(ArgumentError("shifts must be finite"))
    all(isfinite, qv) || throw(ArgumentError("q0 must be finite"))
    validate_spectrum(kv, bv)

    return CrumProblem{T}(xv, kv, bv, sv, qv)
end

"""
    effective_beta(problem)

Return norming constants after applying per-eigenvalue shifts:
`beta_n * exp(2*kappa_n*shift_n)`.
"""
effective_beta(problem::CrumProblem) = shifted_norming_constants(problem.kappa, problem.beta, problem.shifts)

"""
    validate_spectrum(kappa, beta)

Check the ordering and basic assumptions used by Algorithm 1 in Prins-Wahls:
`0 < kappa[1] < ... < kappa[N]`, and real nonzero norming constants.
"""
function validate_spectrum(kappa::AbstractVector, beta::AbstractVector)
    length(kappa) == length(beta) || throw(ArgumentError("kappa and beta must have equal length"))
    all(isfinite, kappa) || throw(ArgumentError("kappa values must be finite"))
    all(isfinite, beta) || throw(ArgumentError("beta values must be finite"))
    all(>(0), kappa) || throw(ArgumentError("all kappa values must be positive"))
    all(!iszero, beta) || throw(ArgumentError("all beta values must be nonzero"))
    if length(kappa) > 1
        all(diff(kappa) .> 0) || throw(ArgumentError("kappa must be strictly increasing"))
    end
    return true
end

"""
    has_zero_background(q0; atol=0)

Return true when the sampled background is zero up to `atol`.
"""
has_zero_background(q0::AbstractVector; atol=0) = all(abs.(q0) .<= atol)

"""
Metadata for one paired or onefold Crum step.
"""
struct CrumStepRecord{T<:AbstractFloat}
    fold::Int
    kappa::Vector{T}
    min_q::T
    max_q::T
end

"""
Result wrapper returned when `return_history=true`.
"""
struct CrumTransformResult{T<:AbstractFloat}
    x::Vector{T}
    kappa::Vector{T}
    beta::Vector{T}
    effective_beta::Vector{T}
    shifts::Vector{T}
    q0::Vector{T}
    q::Vector{T}
    history::Vector{CrumStepRecord{T}}
end
