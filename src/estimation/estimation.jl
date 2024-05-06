using DataFrames, DataFramesMeta, CSV, FixedEffectModels, Optim

include("stage1.jl")
include("stage2.jl")
include("stage3.jl")
include("stage4.jl")
include("stage5.jl")
include("input_output.jl")

# this function combines each of the stages, written above.
function estimate_model(x4_0,x5_0,θ,θk,V,F,data,panel,kid_data; R = 10, num_iter = 1000, show_trace = true, seed4 = 1234, seed5 = 20240220)
    θ = stage1(θ, F, panel)
    θ = stage2(θ, kid_data)
    θ = stage3(θ, kid_data, F.τgrid, θ.cprobs)
    θ = stage4(x4_0,θ,V,F,data,panel; R, num_iter, show_trace, seed = seed4)
    dat = prep_sim_data(data, panel; R)
    mod = (;θ,values=V,F)
    θk = stage5(x5_0, θk, mod, dat, kid_data; num_iter=1000, show_trace, seed = seed5)
    return θ,θk
end

include("bootstrap.jl")