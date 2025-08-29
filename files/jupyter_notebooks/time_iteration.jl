using Interpolations, LinearAlgebra, Optim, Plots, Roots
using Parameters


const β = 0.96
const α = 0.4

const kmin = 1e-3;
const kmax = 100.0;
const nk = 1_001;
const kgrid = range(kmin, kmax, length=nk)

f(k) = k^α
u(c) = c ≥ 0.0 ? log(c) : -1e10
u_c(c) = c ≥ 0.0 ? 1/c : 0.0
f_k(k) = α * k^(α - 1)
u_c_inv(c) = c ≥ 0.0 ? 1/c : -1e10

function coleman_operator_roots(polguess_V::Vector{Float64})
    polguess_F = CubicSplineInterpolation(kgrid, polguess_V)
    Kp_V = zeros(nk)
    for i in eachindex(kgrid)
        k = kgrid[i]
        function euler_eq(k′)
            LHS = u_c(f(k) - k′)
            c′ = f(k′) - polguess_F(k′)
            RHS = β*f_k(k′)*u_c(c′)
            return LHS - RHS
        end
        Kp_V[i] = find_zero(euler_eq, (kmin, f(k)))
    end
    return Kp_V
end

function solve_model_coleman()
    polguess_V = zeros(nk)
    iiter = 0
    diff = 1
    max_iter = 600
    tol = 1e-8

    @time while iiter < max_iter 
        iiter += 1
        polnew_V = coleman_operator_roots(polguess_V)

        diff = norm(polnew_V - polguess_V, Inf)

        println("Iteration $iiter, Error: $diff")

        if diff < tol
            println("Converged after $iiter iterations")
            break
        end

        polguess_V = polnew_V
    end
    return polguess_V
end

pol_V = solve_model_coleman()

plot(kgrid, pol_V, label="Coleman Time Iteration", 
    xlabel="Capital", ylabel="Policy Function", 
    title="Policy Function from Coleman Time Iteration")

pol_analytical(k) = α*β*k^α

plot!(kgrid, pol_analytical.(kgrid), label="Analytical Policy Function", 
    linestyle=:dash, color=:red)

pol_F = CubicSplineInterpolation(kgrid, pol_V)

function compute_euler_residual(k::Float64)
    k′ = pol_F(k)
    c = f(k) - k′
    LHS = u_c(c)
    k′′ = pol_F(k′)
    c′ = f(k′) - k′′
    RHS = β * f_k(k′) * u_c(c′)
    return log10(abs(LHS/RHS - 1))
end

kfine_V = range(kmin, kmax, length=5_001)

errors_V = compute_euler_residual.(kfine_V)

#histogram(errors_V)
using Statistics

mean(errors_V)

function compute_euler_residual_improved(k::Float64)
    k′ = pol_F(k)
    c = f(k) - k′
    k′′ = pol_F(k′)
    c′ = f(k′) - k′′
    RHS = β * f_k(k′) * u_c(c′)
    return log10(abs(1 - (u_c_inv(RHS)/c)))
end

errors_V_improved = compute_euler_residual_improved.(kfine_V)

#histogram(errors_V)
using Statistics

mean(errors_V)




function coleman_operator_optim(polguess_V::Vector{Float64})
    polguess_F = CubicSplineInterpolation(kgrid, polguess_V)
    Kp_V = zeros(nk)
    for i in eachindex(kgrid)
        k = kgrid[i]
        function euler_eq(k′)
            LHS = u_c(f(k) - k′)
            c′ = f(k′) - polguess_F(k′)
            RHS = β*f_k(k′)*u_c(c′)
            return LHS - RHS
        end
        results = optimize(k′ -> euler_eq(k′)^2, kmin, f(k))
        Kp_V[i] = results.minimizer
    end
    return Kp_V
end

