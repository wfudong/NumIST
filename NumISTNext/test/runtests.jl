using Test
using NumISTNext

@testset "Crum one-soliton" begin
    x = collect(range(-12.0, 12.0; length=2001))
    kappa = [1.4]
    beta = [1.0]

    q = crum_transform(x, kappa, beta)
    q_exact = exact_one_soliton(x, kappa[1]; beta=beta[1])

    @test relative_l2_error(q, q_exact) < 1e-10
    @test validate_spectrum(kappa, beta)
end

@testset "KdV ETDRK4 zero solution" begin
    x = collect(range(-π, π; length=65))[1:end-1]
    q0 = zeros(length(x))
    sol = solve_kdv_etdrk4(q0, x, [0.0, 0.01, 0.02]; dt=0.01)

    @test all(iszero, sol[0.0])
    @test maximum(abs, sol[0.01]) < 1e-14
    @test maximum(abs, sol[0.02]) < 1e-14
end

@testset "Paper example helpers" begin
    @test paper_example91_initial(0.0) == 0.0
    @test paper_example92_exact(0.0, 0.0; c=π) ≈ -π / 2
    @test paper_example93_initial(-1.0) ≈ exp(-1.0) * cos(-4.0)
    @test section93_time_slices().example93 == [0.0, 0.015, 0.03]
end

@testset "Section 9 finite-section IST" begin
    data = paper_example92_scattering_data(c=π)
    x = collect(range(-5.0, 7.0; length=501))
    interior = 5:(length(x)-4)

    for t in (0.0, 1.0)
        result = kdv_ist_solution(data.right, data.left, x, t; Ns=8, Ni=64)
        exact = paper_example92_exact.(x, t)
        @test maximum(abs.(result.q[interior] .- exact[interior])) < 2e-3
    end
end
