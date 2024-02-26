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
θk = prod_pars()

x4_0 = readdlm("output/Xsave_round2")[:,1]
x5_0 = [-20.,0.75,0.75,0.05,0.00,2.5,log(0.2),log(0.7),log(20.)]

θ,θk = estimate_model(x4_0, x5_0, θ, θk, V, F, M, P ; R = 10, num_iter = 1000, show_trace = true)
x1,x2,x3,x4,x5 = stack_ests(θ,θk)
writedlm("output/est_stage1",x1)
writedlm("output/est_stage2",x2)
writedlm("output/est_stage3",x3)
writedlm("output/est_stage4",x4)
writedlm("output/est_stage5",x5)
# next: run this again to make sure it all works.

# then: how shall we write up the results? in a jupyter notebook? or just go for the tables?