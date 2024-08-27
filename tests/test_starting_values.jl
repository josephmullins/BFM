include("../src/model/model.jl")
include("../src/estimation/estimation.jl")
using DelimitedFiles

M = CSV.read("data/MarriageFile.csv",DataFrame,missingstring="NA")
P = CSV.read("data/MarriagePanel.csv",DataFrame,missingstring="NA")
K = CSV.read("data/KidPanelv2.csv",DataFrame,missingstring="NA")
cprobs = CSV.read("data/CustodyMomentsSimple.csv",DataFrame).p[:]

F = FixedParams()
θ = Params(F)
V = values(F);
θ = (;θ...,cprobs)
θk = prod_pars()

#x4_0 = readdlm("output/Xsave_round2_old")[:,1] #<- can't do better than this.

θ = stage1(θ, F, P)
θ = stage2(θ, K)
θ = stage3(θ, K, F.τgrid, θ.cprobs)

X4 = readdlm("output/Xsave_round2")

dat = prep_sim_data(M,P;R = 10)
moms0 = data_moms(M,P)

for i in 1:10
    println("Doing guess number $i")
    res1 = optimize(x->ssq(update(x,θ,F),V,F,moms0,dat),X4[:,i],Optim.Options(iterations=100,show_trace=false))
    @show res1.minimum
end