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

x0 = readdlm("output/stage5")[:]
θk =  updateθk(x0,θk,θ)

x1,x2,x3,x4,x5 = stack_ests(θ,θk);
θb,θkb = update_all((x1,x2,x3,x4,x5),θ,θk,F);

ma = get_moments(θ,V,F,dat)
mb = get_moments(θb,V,F,dat)
[ma mb]

kmomsa = kid_moments(S,θk,θ,F,kid_data);
kmomsb = kid_moments(S,θkb,θb,F,kid_data);
[kmomsa[1] kmomsb[1]]
[kmomsa[2] kmomsb[2]]
[kmomsa[3] kmomsb[3]]
