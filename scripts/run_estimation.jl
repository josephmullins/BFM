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

# θ = stage1(θ, F, P)
# θ = stage2(θ, K)
# θ = stage3(θ, K, F.τgrid, θ.cprobs)
# θ = stage4(x4_0,θ,V,F,M,P; R = 10, num_iter = 5)
# dat = prep_sim_data(M, P; R = 10)
# mod = (;θ,values=V,F)
# sim_data = data_gen(mod,dat);
# kid_data = prep_child_data(sim_data,dat,cprobs);
# S = zeros(length(kid_data.ΩK))

# kmoms0 = kidmoms_data(K)

# kmoms1 = kid_moments(S,θk,θ,F,kid_data)

# θk = stage5(x5_0, θk, mod, dat, K; num_iter = 5)



θ,θk = estimate_model(x4_0, x5_0, θ, θk, V, F, M, P ; R = 10, num_iter = 5, show_trace = true)