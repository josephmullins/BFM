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
L = convert(Vector{Int64},P.L)
M[!,:age_mar] = M.YMAR .- M.YBIRTH_M
age_mar = convert(Vector{Int64},M.YMAR .- M.YBIRTH_M)
tlength = convert(Vector{Int64},M.tlength)
moms0 = data_moms(P).value

dat = repeat(M,10)
legal = repeat(L,10)

θ = (;θ...,
    σ_L = 2., α_l = 1.,
    α_F = -2.,σ_F = 10.,
    α_ω = [0., 5.],σ_ω = 3., # σ_ω = 2.
    π_ω = 0.7, Π_ω = transmat_ω(0.7,F.N_ω),
    α_νH0 = [3.,-0.1],α_νH1 = [7.,0.2])


#wght = [1.,1.,1.,1.,0.,0.,1.,0.,0.,1.,50.,0.,0.,1.]
x0 = get_x(θ)

res1 = optimize(x->ssq(update(x,θ,F),V,F,moms0,dat,legal),x0,Optim.Options(iterations=1000,show_trace=true))
res1 = optimize(x->ssq(update(x,θ,F),V,F,moms0,dat,legal),res1.minimizer,Optim.Options(iterations=1000,show_trace=true))
θ = update(res1.minimizer,θ,F)
m1 = get_moments(θ,V,F,M,L)
#display([moms0 m1])

using Plots
M1 = reshape(m1,19,5)
M0 = reshape(moms0,19,5)
plot!(22:40,M1,layout = (2,3))
plot!(22:40,M0,layout = (2,3))
