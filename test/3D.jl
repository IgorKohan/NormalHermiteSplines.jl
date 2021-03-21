@testset "Test 3D" begin

    # function get_3D_model(p::Vector{Float64})
    #     r = p[1] + p[2] + p[3]
    #     return r
    # end
    #
    # function get_3D_model_grad()
    #     return [1.0; 1.0; 1.0]
    # end

    p = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0] # function nodes
    u = [1.0; 1.0; 1.0; 0.0]                                # function values in nodes
    n_1 = size(p, 2)
    dp = Matrix{Float64}(undef, 3, 3*n_1)
    es = Matrix{Float64}(undef, 3, 3*n_1)
    du = Vector{Float64}(undef, 3*n_1)
    grad = [1.0; 1.0; 1.0]

    k = 0
    for i = 1:n_1
        k += 1
        dp[1,k] = p[1,i]
        dp[2,k] = p[2,i]
        dp[3,k] = p[3,i]
        du[k] = grad[1]
        es[1,k] = 1.0
        es[2,k] = 0.0
        es[3,k] = 0.0
        k += 1
        dp[1,k] = p[1,i]
        dp[2,k] = p[2,i]
        dp[3,k] = p[3,i]
        du[k] = grad[2]
        es[1,k] = 0.0
        es[2,k] = 1.0
        es[3,k] = 0.0
        k += 1
        dp[1,k] = p[1,i]
        dp[2,k] = p[2,i]
        dp[3,k] = p[3,i]
        du[k] = grad[3]
        es[1,k] = 0.0
        es[2,k] = 0.0
        es[3,k] = 1.0
    end
    ####

    @testset "Test 3D-RK_H0 kernel" begin
        s = interpolate(p, u, RK_H0(0.00001))
        σ1 = evaluate_one(s, [1.0; 0.0; 0.01])
        @test isapprox(σ1, u[1], atol = 1e-3)
        σ = evaluate(s, [1.0 0.5; 0.0 0.5; 0.01 0.5])
        @test isapprox(σ[1], u[1], atol = 1e-3)

        s = prepare(p, RK_H0(0.00001)) # prepare spline
        s = construct(s, u) # construct spline
        σ1 = evaluate_one(s, [1.0; 0.0; 0.01])
        @test isapprox(σ1, u[1], atol = 1e-3)
        σ = evaluate(s, [1.0 0.5; 0.0 0.5; 0.01 0.5])
        @test isapprox(σ[1], u[1], atol = 1e-3)

        cond = estimate_cond(s)
        @test cond ≈ 1.0e6

        iq = estimate_accuracy(s)
        @test iq ≈ 15.0

        eps = get_epsilon(s)
        @test eps ≈ 1.0e-5

        est_eps = estimate_epsilon(p)
        @test isapprox(est_eps, 1.036, atol = 1e-3)
    end

    @testset "Test 3D-RK_H1 kernel" begin
        s = interpolate(p, u, dp, es, du, RK_H1(0.1))
        σ1 = evaluate_one(s, [1.0; 0.0; 0.01])
        @test isapprox(σ1, u[1], atol = 1e-2)
        σ = evaluate(s, [1.0 0.5; 0.0 0.5; 0.01 0.5])
        @test isapprox(σ[1], u[1], atol = 1e-2)

        s = prepare(p, dp, es, RK_H1(0.1)) # prepare spline
        s = construct(s, u, du) # construct spline
        σ1 = evaluate_one(s, [1.0; 0.0; 0.01])
        @test isapprox(σ1, u[1], atol = 1e-2)
        σ = evaluate(s, [1.0 0.5; 0.0 0.5; 0.01 0.5])
        @test isapprox(σ[1], u[1], atol = 1e-2)

        g = evaluate_gradient(s, [0.01; 0.2; 0.3])
        @test all(isapprox.(g, grad, atol = 0.01))

        cond = estimate_cond(s)
        @test cond ≈ 1.0e5

        iq = estimate_accuracy(s) + 1.0
        @test iq ≈ 14.0 atol = 1.0

        eps = get_epsilon(s)
        @test eps ≈ 0.1

        est_eps = estimate_epsilon(p)
        @test isapprox(est_eps, 1.036, atol = 1e-3)

        est_eps = estimate_epsilon(p, dp)
        @test isapprox(est_eps, 1.415, atol = 1e-3)
    end

    @testset "Test 3D-RK_H2 kernel" begin
        s = interpolate(p, u, dp, es, du, RK_H2(0.1))
        σ1 = evaluate_one(s, [1.0; 0.0; 0.01])
        @test isapprox(σ1, u[1], atol = 2e-2)
        σ = evaluate(s, [1.0 0.5; 0.0 0.5; 0.01 0.5])
        @test isapprox(σ[1], u[1], atol = 2e-2)

        s = prepare(p, dp, es, RK_H2(0.1)) # prepare spline
        s = construct(s, u, du) # construct spline
        σ1 = evaluate_one(s, [1.0; 0.0; 0.01])
        @test isapprox(σ1, u[1], atol = 2e-2)
        σ = evaluate(s, [1.0 0.5; 0.0 0.5; 0.01 0.5])
        @test isapprox(σ[1], u[1], atol =2e-2)

        g = evaluate_gradient(s, [0.01; 0.2; 0.3])
        @test all(isapprox.(g, grad, atol = 0.01))

        cond = estimate_cond(s)
        @test cond ≈ 1.0e8

        iq = estimate_accuracy(s) + 1.0
        @test iq ≈ 13.0 atol = 1.0

        eps = get_epsilon(s)
        @test eps ≈ 0.1

        est_eps = estimate_epsilon(p)
        @test isapprox(est_eps, 1.036, atol = 1e-3)

        est_eps = estimate_epsilon(p, dp)
        @test isapprox(est_eps, 1.415, atol = 1e-3)
    end

end
