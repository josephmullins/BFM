using DataFrames, DataFramesMeta, CSV, GLM, Optim

include("stage1.jl")
include("stage2.jl")
include("stage3.jl")
include("stage4.jl")
include("stage5.jl")
include("input_output.jl")

# this function combines each of the stages, written above.
function estimate_model(x4_0,x5_0,θ,θk,V,F,data,panel; R = 10, num_iter = 1000, show_trace = true)
    θ = stage1(θ, F, P)
    θ = stage2(θ, K)
    θ = stage3(θ, K, F.τgrid, θ.cprobs)
    θ = stage4(x4_0,θ,V,F,data,panel; R, num_iter, show_trace)
    dat = prep_sim_data(data, panel; R)
    mod = (;θ,values=V,F)
    θk = stage5(x5_0, θk, mod, dat, K; num_iter, show_trace)
    return θ,θk
end
