

using Interpolations, LinearAlgebra, Optim, Plots, Roots
using Parameters

@with_kw struct Params{F1, F2, F3, F4}
    #Economic parameters
    β::Float64  # Discount factor
    α::Float64   # Capital share in production

    #Grid
    kmin::Float64
    kmax::Float64
    nk::Int64
    kgrid::AbstractVector{<:Real}

    #Numerical parameters
    tol::Float64
    max_iter::Int64
    
    #Functions
    u::F1
    f::F2
    u_c::F3
    f_k::F4
end

function make_params(; β = 0.96, α = 0.4, kmin = 1e-3, kmax = 100.0, nk = 1_001, tol = 1e-8, max_iter = 600)
    kgrid = range(kmin, kmax, length=nk)
    u(c) = c ≥ 0.0 ? log(c) : -1e10
    f(k) = k^α
    u_c(c) = c ≥ 0.0 ? 1/c : -1e10
    f_k(k) = α * k^(α - 1)
    return Params(β = β, α = α, kmin = kmin, kmax = kmax, nk = nk, kgrid = kgrid, tol = tol, max_iter = max_iter, u = u, f = f, u_c = u_c, f_k = f_k)
end

p = make_params()

function coleman_operator(polV_old::Vector{Float64}; p::Params)
    @unpack kgrid, nk, kmin, f, β, f_k, u_c = p
    polF_old = CubicSplineInterpolation(kgrid, polV_old)
    polV_new = zeros(nk)
    for i in eachindex(kgrid)
        k = kgrid[i]
        function euler_eq(k′)
            LHS = u_c(f(k) - k′)
            c′ = f(k′) - polF_old(k′)
            RHS = β*f_k(k′)*u_c(c′)
            return LHS - RHS
        end
        polV_new[i] = find_zero(euler_eq, (kmin, f(k)))
    end
    return polV_new
end

function solve_model_coleman(p::Params)
    @unpack tol, max_iter, nk = p

    polV_old = zeros(nk)
    iiter = 0
    diff = 1

    @time while iiter < max_iter 
        iiter += 1
        polV_new = coleman_operator(polV_old, p = p)

        diff = norm(polV_new - polV_old, Inf)
        println("Iteration $iiter, Error: $diff")

        polV_old = copy(polV_new)

        if diff < tol
            println("Converged after $iiter iterations")
            break
        end
    end
    return polV_old
end

polV = solve_model_coleman(p);

