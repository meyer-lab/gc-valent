using OrdinaryDiffEq
using DiffEqSensitivity

include("reaction.jl")

function runCkine(tps, params, IL2case)
    if IL2case
        _, _, _, trafP = IL2param(params)
        f = IL2Deriv
    else
        _, _, _, trafP = fullParam(params)
        f = fullDeriv
    end

    u0 = solveAutocrine(trafP)

    prob = ODEProblem(f, u0, (0.0, maximum(tps)), params)

    sol = solve(prob, Rodas5())

    adjsol = adjoint_sensitivities(sol, Rodas5(), adjG, nothing, adjGd)

    return sol(tps)
end
