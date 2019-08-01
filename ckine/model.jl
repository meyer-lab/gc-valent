using DifferentialEquations
using Sundials
using DiffEqSensitivity
    
include("reaction.jl")

function runCkine(tps, params, sensi, IL2case)
    if IL2case
        _, _, _, trafP = IL2param(params)
        f = IL2Deriv
    else
        _, _, _, trafP = fullParam(params)
        f = fullDeriv
    end

    u0 = solveAutocrine(trafP)

    if sensi
        prob = ODELocalSensitivityProblem(f, u0, (0.0, maximum(tps)), params)
    else
        prob = ODEProblem(f, u0, (0.0, maximum(tps)), params)
    end

    sol = solve(prob, Vern9())

    if sensi
        return extract_local_sensitivities(sol, tps)
    else
        return sol(tps), nothing
    end
end
