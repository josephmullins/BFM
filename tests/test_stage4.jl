include("../src/model/model.jl")
include("../src/estimation/estimation.jl")

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

# model_moms(TD,TF,D,L_sim)
# have to re-pack this model each time. I don't like that. Change the arguments?
# mod = (;F,θ,values=V);
# solve_all!(mod)
# model_moms(TD,TF,D,L_sim) 
moms0 = data_moms(M,P)
ssq(θ,V,F,moms0,M,P.L)
# NEXT: write an update function for the parameters.

x0 = get_x(θ)
optimize(x->ssq(update(x,θ,F),V,F,moms0,M,P.L),x0,Optim.Options(iterations=10))