export besselj0_quad,
       paper_example91_initial,
       paper_example92_exact,
       paper_example93_initial,
       section93_time_slices

_sech(x) = inv(cosh(x))

"""
    besselj0_quad(x; ntheta=512)

Evaluate `J0(x)` by the integral representation used in the old notebook and
script: `pi^-1 * integral_0^pi cos(x sin(theta)) dtheta`.
"""
function besselj0_quad(x; ntheta::Int=512)
    ntheta > 0 || throw(ArgumentError("ntheta must be positive"))
    h = π / ntheta
    s = 0.5 * (cos(x * sin(0.0)) + cos(x * sin(π)))
    for j in 1:ntheta-1
        θ = j * h
        s += cos(x * sin(θ))
    end
    return h * s / π
end

paper_example91_initial(x) = x * exp(-x^2)

function paper_example92_exact(x, t; c=π)
    return -c / 2 * _sech(sqrt(c) * (x - c * t) / 2)^2
end

function paper_example93_initial(x)
    if x < 0
        return exp(x) * cos(4x)
    end
    return exp(-x) * besselj0_quad(2x)
end

function section93_time_slices()
    return (;
        example91=[0.0, 0.28, 0.4, 0.5, 0.6],
        example92=[0.0, 0.5, 1.0],
        example93=[0.0, 0.015, 0.03],
    )
end
