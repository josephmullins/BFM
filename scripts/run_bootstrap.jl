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

X1b,X2b,X3b,X4b,X5b = bootstrap_pars(θ,θk,V,F,M,P,K,50 ; R = 10, num_iter = 1, show_trace = false)

break
X1c,X2c,X3c,X4c,X5c = bootstrap_pars(θ,θk,V,F,M,P,K,5 ; R = 10, num_iter = 1, show_trace = false)



writedlm("output/boot_stage1",X1b)
writedlm("output/boot_stage2",X2b)
writedlm("output/boot_stage3",X3b)
writedlm("output/boot_stage4",X4b)
writedlm("output/boot_stage5",X5b)

# something random is happening?
seed0 = 1010
seed1 = 2020
seed2 = 3030
x4_0 = get_x(θ)
x5_0 = get_xk(θk)
x1,x2,x3,x4,x5 = stack_ests(θ,θk)

b = 5
Random.seed!(seed0+b)
Mb,Pb,Kb = draw_boot_sample(M,P,K)
seed4 = seed1+b
seed5 = seed2+b
θ,θk = estimate_model(x4_0, x5_0, θ, θk, V, F, Mb, Pb, Kb ; R = 10, num_iter=1, show_trace=false, seed4, seed5)
x1b,x2b,x3b,x4b,x5b = stack_ests(θ,θk)


# let's unpack a little more:
θ = stage1(θ, F, Pb)
θ = stage2(θ, Kb)
θ = stage3(θ, Kb, F.τgrid, θ.cprobs)
θ = stage4(x4_0,θ,V,F,Mb,Pb; R = 10, num_iter = 1, show_trace = true, seed = seed4)
dat = prep_sim_data(Mb, Pb; R)
mod = (;θ,values=V,F)
θk = stage5(x5_0, θk, mod, dat, Kb; num_iter = 1, show_trace = true, seed = seed5)

solve_all!(mod)
sim_data = data_gen(mod,dat);

kid_data = prep_child_data(sim_data,dat,cprobs);

S = zeros(length(kid_data.ΩK))

# the data moments:
kmoms0 = kidmoms_data(K)
predict_k!(S,θk,θ,F,kid_data ; seed = seed5)