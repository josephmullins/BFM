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

mod = (;θ,values=V,F)

dat = prep_sim_data(M,P;R = 10)
moms0 = data_moms(M,P)

X0 = readdlm("output/Xsave_round2")

# this produces the best outcome from all the guesses
res1 = optimize(x->ssq(update(x,θ,F),V,F,moms0,dat),X0[:,1],Optim.Options(iterations=1000,show_trace=true))

θ = update(res1.minimizer,θ,F)
m1 = get_moments(θ,V,F,dat)
display([moms0 m1])

