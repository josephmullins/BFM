include("../src/model/model.jl")
include("../src/estimation/estimation.jl")
using DelimitedFiles

M = CSV.read("data/MarriageFile.csv",DataFrame,missingstring="NA")
P = CSV.read("data/MarriagePanel.csv",DataFrame,missingstring="NA")
K = CSV.read("data/KidPanel.csv",DataFrame,missingstring="NA")
cprobs = CSV.read("data/CustodyMomentsSimple.csv",DataFrame).p[:]

F = FixedParams()
θ = Params(F)
V = values(F);
θ = (;θ...,cprobs)

# this works, but we need to set up the sequence.
θ = stage1(θ,F,P)
θ = stage2(θ,K)
θ = stage3(θ,K,F.τgrid,cprobs)

L = convert(Vector{Int64},P.L)
moms0 = data_moms(M,P)

dat = prep_sim_data(M,P;R = 10)

x4 = readdlm("output/stage4")[:]
θ = update(x4,θ,F)

mod = (;θ,values=V,F);

solve_all!(mod)
sim_data = data_gen(mod,dat);

kid_data = prep_child_data(sim_data,dat,cprobs);

S = zeros(length(kid_data.ΩK))
θk = (;γ_ψ0 = -10,
    γ_ψ = [0.,0.,0.5,0.],
    δW = fill(0.1,18),δH = fill(0.1,18),δk = 0.9,
    Γa = zeros(19))

θk =  updateθk(x0,θk,θ)

# the data moments:
kmoms0 = kidmoms_data(K)

x0 = [-20.,0.75,0.75,0.05,0.00,2.5,log(0.2),log(0.7),log(20.)]

res = optimize(x->obj_stage5(S,updateθk(x,θk,θ),θ,F,kid_data,kmoms0),x0,Optim.Options(show_trace = true))
θk =  updateθk(res.minimizer,θk,θ)

writedlm("output/stage5",res.minimizer)

x1,x2,x3,x4,x5 = stack_ests(θ,θk)
writedlm("output/est_stage1",x1)
writedlm("output/est_stage2",x2)
writedlm("output/est_stage3",x3)
writedlm("output/est_stage4",x4)
writedlm("output/est_stage5",x5)