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

x4_0 = readdlm("output/Xsave_round2_old")[:,1] #<- can't do better than this.
#x5_0 = [-20.,0.75,0.75,0.05,0.00,2.5,log(0.2),log(0.7),log(20.)]
x5_0 = [-15.,0.75,0.75,0.05,0.00,2.5,log(0.2),log(0.7),log(12.)] #<- slightly better guess

θ = stage1(θ, F, P)
θ = stage2(θ, K)
θ = stage3(θ, K, F.τgrid, θ.cprobs)

X4 = readdlm("output/Xsave_round1")

dat = prep_sim_data(M,P;R = 10)
moms0 = data_moms(M,P)

seed4 = 1234

@show ssq(update(x4_0,θ,F),V,F,moms0,dat)

for i in axes(X4,2)
    @show ssq(update(X4[:,i],θ,F),V,F,moms0,dat)
end


#θ = stage4(x4_0,θ,V,F,data,panel; R, num_iter, show_trace, seed = seed4)


# θ,θk = estimate_model(x4_0, x5_0, θ, θk, V, F, M, P, K ; R = 10, num_iter = 1000, show_trace = true)
