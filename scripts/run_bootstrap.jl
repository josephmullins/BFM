# load the estimates
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

x1 = readdlm("output/est_stage1")[:]
x2 = readdlm("output/est_stage2")[:]
x3 = readdlm("output/est_stage3")[:]
x4 = readdlm("output/est_stage4")[:]
x5 = readdlm("output/est_stage5")[:]

θ,θk = update_all((x1,x2,x3,x4,x5),θ,θk,F);

X1b,X2b,X3b,X4b,X5b = bootstrap_pars(θ,θk,V,F,M,P,K,50 ; R = 10, num_iter = 100, show_trace = false)

writedlm("output/boot_stage1",X1b)
writedlm("output/boot_stage2",X2b)
writedlm("output/boot_stage3",X3b)
writedlm("output/boot_stage4",X4b)
writedlm("output/boot_stage5",X5b)
