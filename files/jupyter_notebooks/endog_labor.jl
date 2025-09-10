

using Interpolations, LinearAlgebra, Optim, Plots, Roots
using Parameters

@with_kw struct Params{F1, F2, F3, F4, F5, F6}
    #Economic parameters
    β::Float64  # Discount factor
    α::Float64  # Capital share in production
    σ::Float64  # Risk aversion
    ψ::Float64  # Disutility of labor parameter

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
    u_l::F4
    f_k::F5
    f_l::F6
end

function make_params(; β = 0.96, α = 0.4, σ = 2.0, ψ = 0.9, kmin = 1e-3, kmax = 100.0, nk = 1_001, tol = 1e-8, max_iter = 600)
    kgrid = range(kmin, kmax, length=nk)
    u(c,l) = c ≥ 0.0 ? (c^(1-σ))/(1-σ) - (l^(1-ψ))/(1-ψ) : -1e10
    u_c(c) = c ≥ 0.0 ? c^-σ : -1e10
    u_l(l) = -l^ψ
    f(k,l) = (k^α * l^(1-α))
    f_k(k,l) = α * (k/l)^(α - 1)
    f_l(k,l) = (1 - α) * (k/l)^α
    return Params(β = β, α = α, σ = σ, ψ = ψ, 
                kmin = kmin, kmax = kmax, nk = nk, kgrid = kgrid, 
                tol = tol, max_iter = max_iter,
                u = u, f = f, u_c = u_c, u_l = u_l, f_k = f_k, f_l = f_l)
end

p = make_params()

polK_V = zeros(p.nk)
polL_V = zeros(p.nk)

function bellman_operator(valueV_old::Vector{Float64}, p::Params)
    @unpack kgrid, nk, kmin, kmax, f, β, f_k, u, u_c, u_l = p
    valueF_old = CubicSplineInterpolation(kgrid, valueV_old)
    valueV_new = zeros(nk)
    for i in eachindex(kgrid)
        k = kgrid[i]
        println("Solving for state $i out of $nk")
        function objective(x)
            k′, l = x
            if k′ > kmax || k′ < kmin
                return 1e10
            end

            c = f(k,l) - k′

            if c < 0
                return 1e10
            end

            return -(u(c,l) + β * valueF_old(k′)) 
        end
        result = optimize(objective, [0.5, 0.5], NelderMead())
        valueV_new[i] = -result.minimum
        polK_V[i] = clamp(result.minimizer[1], kmin, kmax)
        polL_V[i] = result.minimizer[2]
    end
    return valueV_new
end


valueV_old = zeros(p.nk)

for i in eachindex(p.kgrid)
    k = p.kgrid[i]

end


test = bellman_operator(valueV_old, p)

