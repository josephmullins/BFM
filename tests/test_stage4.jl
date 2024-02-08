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

#x0 = get_x(θ)
#res1 = optimize(x->ssq(update(x,θ,F),V,F,moms0,M,P.L),x0,Optim.Options(iterations=100))

blocks = [1:2,3:4,5:6,7:11]
for b in blocks
    xcurrent = get_x(θ)
    x0 = xcurrent[b]
    res = optimize(x->ssq(update(x,b,θ,F),V,F,moms0,M,P.L),x0,Optim.Options(iterations=30))
    xcurrent[b] .= res.minimizer
    θ = update(xcurrent,θ,F)
end

m1 = get_moments(θ,V,F,M,P.L)
display([moms0 m1])

