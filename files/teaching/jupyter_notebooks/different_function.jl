

using Interpolations, LinearAlgebra, Roots
using Parameters
using NonlinearSolve

const β = 0.96
const α = 0.8

const kmin = 1e-3;
const kmax = 100.0;
const nk = 1_001;
const kgrid = range(kmin, kmax, length=nk)

f(k) = k^α
u(c) = c ≥ 0.0 ? log(c) : -1e10
u_c(c) = c ≥ 0.0 ? 1/c : -1e10
f_k(k) = α * k^(α - 1)
u_c_inv(c) = c ≥ 0.0 ? 1/c : -1e10





function K(g; p)
    #Construct a linear interpolation approximation of the consumption policy function g
    g_func = LinearInterpolation(p.k_grid, g, extrapolation_bc = Line())
    function h!(H, c, k)
        H[1] = p.dudc(c[1]) .- p.β*(p.dudc.(g_func(p.f(k) .- c[1]))*p.dfdk(p.f(k) .- c[1]))
    end
    #Solving for the updated consumption value
    Kg = Vector{Vector{Float64}}(undef, p.nk)
    for (i, k) in collect(enumerate(p.k_grid))
        Kg[i] = nlsolve((H,x) -> h!(H, x, k), [0.02*k]).zero  # find the c that makes the EE 0 in those specific grid-points!
    end
    return map(v -> v[1], Kg) #Converts the vector of 1-element vectors, into just a vector
end




function coleman_operator(cpolV_old::Vector{Float64})
    cpolF_old = CubicSplineInterpolation(kgrid, cpolV_old)
    cpolV_new = zeros(nk)
    for i in eachindex(kgrid)
        k = kgrid[i]
        function euler_eq(du, u, p)
            c = u[1]
            LHS = u_c(c)
            k′ = f(k) - c

            if k′ > f(k)
                du[1] = 1e10
            end

            c′ = cpolF_old(k′)
            RHS = β * f_k(k′) * u_c(c′)
            du[1] = LHS - RHS
        end
        u0 = [cpolF_old(k)]
        prob = NonlinearProblem(euler_eq, u0)
        cpolV_new[i] = solve(prob, NewtonRaphson())[1]
    end
    return cpolV_new
end

cpolV_old = f.(kgrid)

coleman_operator(cpolV_old)

function solve_model_coleman()
    cpolV_old = f.(kgrid)
    iiter = 0
    diff = 1
    max_iter = 600
    tol = 1e-8

    @time while iiter < max_iter 
        iiter += 1
        cpolV_new = coleman_operator_roots(cpolV_old)

        diff = norm(cpolV_new - cpolV_old, Inf)

        println("Iteration $iiter, Error: $diff")

        cpolV_old = cpolV_new

        if diff < tol
            println("Converged after $iiter iterations")
            break
        end
    end
    return cpolV_old
end

cpol_V = solve_model_coleman()

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



