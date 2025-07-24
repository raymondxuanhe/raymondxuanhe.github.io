#Code to replicate QuantEcon's VFI code to solve the optimal growth model
#Original: https://julia.quantecon.org/dynamic_programming/optgrowth.html
#Last updated: June 19 2025

using LinearAlgebra, QuantEcon, NLsolve, Interpolations, Optim
using Distributions, Expectations
using BenchmarkTools, Parameters
using LaTeXStrings, Plots
using Base.Threads
using Optim: maximum, maximizer

@with_kw mutable struct params
    β::Float64 = 0.96    #HH discount factor
    γ::Float64 = 1.0    #HH risk aversion
    u::Function = (c, γ = γ) -> isone(γ) ? log(c) : (c^(1-γ)-1)/(1-γ)
    α::Float64 = 0.4 #Capital production elasticity
    f::Function = (k) -> k^α #Production function
    δ::Float64 = 1.0 #Depreciation rate
    kmin::Float64 = 1e-3 #Minimum capital
    kmax::Float64 = 90.0 #Maximum capital
    nk::Int = 1_001 #Number of gridpoints for capital
    k::LinRange = range(kmin, kmax, length = nk) #Grid for capital stock
end

p = params()

function T(w; tol = 1e-10)
    w_func = LinearInterpolation(p.k, w, extrapolation_bc = Line())
    Tw = similar(w)
    sigma = similar(w)
    for (i, k_val) in enumerate(p.k)
        # solve maximization for each point in k, using k itself as initial condition.
        results  = maximize(c -> p.u(c, 1) + p.β*w_func.(p.f(k_val) - c), tol, k_val)
        Tw[i]    = maximum(results)
        sigma[i] = maximizer(results)
    end
    return (; w = Tw, sigma) # returns named tuple of results
end

function v_star(k)
    c1 = log(1 - p.α * p.β) / (1 - p.β)
    c2 = (p.α * log(p.α * p.β)) / (1 - p.α)
    c3 = 1 / (1 - p.β)
    c4 = 1 / (1 - p.α * p.β)
    return c1 + c2 * (c3 - c4) + c4 * log(k^p.α)
end

w_star = v_star.(p.k)  # evaluate closed form value along grid of k's

w = T(w_star).w # evaluate operator, access Tw results

function solve_optgrowth(initial_w; iterations = 500, m = 3)
    results = fixedpoint(w -> T(w).w, 
        initial_w; 
        iterations, 
        xtol = 1e-10, m, 
        autodiff = :forward)       # Anderson iteration
    v_star = results.zero          #Fixed point
    sigma = T(results.zero).sigma  #Consumption policy function
    return (; v_star, sigma, results)
end

initial_w = 5 * log.((p.k).^p.α); #initial guess

@time sol = solve_optgrowth(initial_w); #Finding the fixed point

v_star_approx = sol.v_star; #Extracting the solution from the code
println("Converged in $(sol.results.iterations) to an ||residuals||_∞ of $(sol.results.residual_norm)")

v_star_approx


plot(p.k, v_star_approx)
