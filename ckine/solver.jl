using gcSolver


@Base.ccallable function runCkineC(tpsIn::Ptr{Cdouble}, nTps::Csize_t, outt::Ptr{Cdouble}, params::Ptr{Cdouble}, paramsLen::Csize_t)::Cint
    vecc = unsafe_wrap(Vector{Float64}, params, paramsLen)
    tps = unsafe_wrap(Vector{Float64}, tpsIn, nTps)
    
    output = runCkine(tps, vecc)
    
    unsafe_copyto!(outt, Ref(output), gcSolver.Nspecies)
    
    return 0
end