using gcSolver

@Base.ccallable Int64 function runCkineC(tps::Vector{Float64}, params::Vector{Float64}, paramsLen::Int64, outt::Ptr{Float64})::Int64
    
    vecc = zeros(paramsLen)
    
    unsafe_copyto!(Ptr(vecc), Ptr(params), paramsLen)
    
    output = runCkine(tps, vecc)
    
    unsafe_copyto!(outt, Ptr(output), gcSolver.Nspecies)
    
    return 0
end