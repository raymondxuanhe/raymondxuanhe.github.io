#Julia code to solve ECO387C PSet 4 Question 2
#Author: Raymond He
#Last updated: October 20, 2025

using Interpolations, LinearAlgebra, NonlinearSolve, Parameters

@with_kw struct Params{F1, F2}
    #Economic parameters
    β::Float64  # Discount factor
    α::Float64   # Capital share in production
    δ::Float64   # Depreciation rate
    γ::Float64 

    #Grid
    kmin::Float64
    kmax::Float64
    nk::Int64
    nz::Int64 #Number of gridpoints for productivity (2 since it's a 2 state Markov process)
    kgrid::AbstractVector{<:Real}
    zL::Float64
    zH::Float64
    zgrid::AbstractVector{<:Real}
    pHH::Float64 #Probability of STAYING in high state
    pLL::Float64 #Probability of STAYING in low state
    pHL::Float64 #Probability of TRANSITIONING from high to low state
    pLH::Float64 #Probability of TRANSITIONING from low to high state
    Ptrans::Matrix{Float64}

    #Numerical parameters
    tol::Float64
    
    #Functions
    u::F1
    f::F2
end

function make_params(; β = 0.96, α = 0.3, δ = 0.06, γ = 1.0, kmin = 1e-3, kmax = 100.0, nz = 2,
                        nk = 1_001, zL = 0.95, zH = 1.05, pHH = 0.9, pLL = 0.4, tol = 1e-6)

    pHL = 1.0 - pHH
    pLH = 1.0 - pLL
    kgrid = range(kmin, kmax, length = nk)
    zgrid = [zL, zH]
    Ptrans = [pLL pLH; pHL pHH]
    u(c,n) = c ≥ 0.0 ? log(c) - (n^(1 + γ))/(1+γ) : -1e10
    f(k,n) = k^α * n^(1 - α)
    return Params(β = β, α = α, δ = δ, γ = γ, kmin = kmin, kmax = kmax, nk = nk, nz = nz,
                zL = zL, zH = zH, pHH = pHH, pLL = pLL, pHL = pHL, pLH = pLH, tol = tol,
                Ptrans = Ptrans, kgrid = kgrid, zgrid = zgrid, u = u, f = f)
end

@with_kw mutable struct Solution
    vA::Array{Float64, 2} #value function
    polCA::Array{Float64, 2} #policy function for consumption
    polKA::Array{Float64,2} #policy function for capital
    polNA::Array{Float64,2} #policy function for labor
end

p = make_params(zL = 0.85, zH = 1.15)

function solve_period_T_problem(p::Params)
    """
    Initializes all the desired equilibrium objects to correspond to their period T values.
    """
    @unpack nk, nz, zgrid, kgrid = p
    @unpack α, γ, δ, u = p
    polKA = zeros(nk,nz)
    polCA = zeros(nk,nz)
    polNA = ones(nk,nz)
    vA = zeros(nk,nz)

    for iz in 1:nz
        z = zgrid[iz]
        for ik in 1:nk
            k = kgrid[ik]

            function objective(du, u, p)
                k = p[1]
                z = p[2]
                n = u[1]
                numerator = (1 - α) * z * k^α
                denominator = z * (k^α) * n^(1-α) + (1 - δ)*k
                du[1] = n^(γ + α) - (numerator / denominator)
            end

            prob = NonlinearProblem(objective, [0.5], [k,z])
            result = solve(prob)
            if result.retcode != :Success
                println("Solver did not converge at state = ($k, $z)")
            end
            n = result.u[1]
            polNA[ik, iz] = n
            polCA[ik, iz] = z * k^α * n^(1 - α) + (1 - δ)*k
            vA[ik, iz] = u(polCA[ik, iz], polNA[ik, iz])
        end
    end
    
    return Solution(vA = vA, polCA = polCA, 
                    polKA = polKA, polNA = polNA)
end

sol = solve_period_T_problem(p)

function time_iteration_step(polCF, polNF, p::Params)
    @unpack kgrid, zgrid, nk, nz, Ptrans = p
    @unpack β, α, δ, γ, u = p

    polCA_new = zeros(nk,nz)
    polKA_new = zeros(nk,nz)
    polNA_new = zeros(nk,nz)
    vA_new = zeros(nk,nz)

    for iz in 1:nz
        z = zgrid[iz]
        for ik in 1:nk
            k = kgrid[ik]
            function objective(du, u, p)
                k = p[1]
                z = p[2]
                c = u[1]

                n = ((1 - α) * z * k^α/c)^(1/ (γ + α))
                k′ = z * k^α * n^(1 - α) + (1 - δ) * k - c
                
                LHS = 1/c

                RHS = 0.0
                prob = Ptrans[iz, :]
                
                for j in 1:nz
                    z′ = zgrid[j]
                    n′ = polNF(k′, z′)
                    c′ = polCF(k′, z′)
                    RHS += prob[j] * (z′ * α * (k′/n′)^(α - 1) + 1 - δ)/c′
                end
                RHS *= β 
                du[1] = LHS - RHS
            end

            prob = NonlinearProblem(objective, [0.05 * k], [k, z])
            result = solve(prob)

            if result.retcode != :Success
                println("Solver did not converge at state = ($k, $z)")
            end
            c = result.u[1]
            polCA_new[ik, iz] = c
            n = ((1 - α) * z * k^α/c)^(1/ (γ + α))
            polKA_new[ik, iz] = z * k^α * n^(1 - α) + (1 - δ) * k - c
            polNA_new[ik, iz] = n
            vA_new[ik, iz] = u(c, n)
        end
    end
    return polCA_new, polKA_new, polNA_new, vA_new
end

function solve_model!(sol::Solution, p::Params, max_iter::Int64)
    @unpack kgrid, zgrid, tol = p

    iiter = 0

    while iiter < max_iter
        iiter += 1
        
        polNF = LinearInterpolation((kgrid, zgrid), sol.polNA)
        polCF = LinearInterpolation((kgrid, zgrid), sol.polCA)

        polCA_new, polKA_new, polNA_new, vA_new = time_iteration_step(polCF, polNF, p)

        diffV = norm(sol.vA - vA_new, Inf)
        diffC = norm(sol.polCA - polCA_new, Inf)
        diffK = norm(sol.polKA - polKA_new, Inf)
        diffN = norm(sol.polNA - polNA_new, Inf)

        println("Iteration $iiter: diffV = $diffV, diffC = $diffC, diffK = $diffK, diffN = $diffN")

        if diffV < tol
            println("Converged after $iiter iterations.")
            break
        end

        sol.polCA .= polCA_new
        sol.polKA .= polKA_new
        sol.polNA .= polNA_new
        sol.vA .= vA_new
    end
    return sol
end

solve_model!(sol, p, 1000)

using Plots

@unpack kgrid, zgrid = p

polCF_low = LinearInterpolation(kgrid, sol.polCA[:, 1])
polCF_high = LinearInterpolation(kgrid, sol.polCA[:, 2])
polNF_low = LinearInterpolation(kgrid, sol.polNA[:, 1])
polNF_high = LinearInterpolation(kgrid, sol.polNA[:, 2])
polKF_low = LinearInterpolation(kgrid, sol.polKA[:, 1])
polKF_high = LinearInterpolation(kgrid, sol.polKA[:, 2])

default(fontfamily = "Computer Modern")

plot(kgrid, polCF_low.(kgrid), title = "Policy function for consumption", xlabel = "Capital - k", label = "Low")
plot!(kgrid, polCF_high.(kgrid), label = "High")

plot(kgrid, polNF_low.(kgrid), title = "Policy function for labor", xlabel = "Capital - k", label = "Low")
plot!(kgrid, polNF_high.(kgrid), label = "High")

plot(kgrid, polKF_low.(kgrid), title = "Policy function for k′", xlabel = "Capital - k", label = "Low")
plot!(kgrid, polKF_high.(kgrid), label = "High")

