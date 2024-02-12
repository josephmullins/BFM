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
moms0 = data_moms(M,P)

# θ = (;θ...,
#     σ_L = 2., α_l = 1.,
#     α_F = 2.,σ_F = 10.,
#     α_ω = [1., 5.],σ_ω = 4., # σ_ω = 2.
#     π_ω = 0.9, Π_ω = transmat_ω(0.9,F.N_ω),
#     α_νH0 = [3.,0.],α_νH1 = [4.,0.1])

θ = (;θ...,
    σ_L = 2., α_l = 1.,
    α_F = -2.,σ_F = 10.,
    α_ω = [0., 5.],σ_ω = 4., # σ_ω = 2.
    π_ω = 0.7, Π_ω = transmat_ω(0.7,F.N_ω),
    α_νH0 = [3.,-0.1],α_νH1 = [7.,0.2])


wght = [1.,1.,1.,1.,0.,0.,1.,0.,0.,1.,10.,0.,0.,1.]
x0 = get_x(θ)
res1 = optimize(x->ssq(update(x,θ,F),wght,V,F,moms0,M,L),x0,Optim.Options(iterations=200))

#V2 = values(F);
#     return (;θ..., α_l=exp(x[1]), σ_L=exp(x[2]),α_F=x[3], σ_F=exp(x[4]), α_ω=[x[5],exp(x[6])], α_νH0=x[7:8], α_νH1=x[9:10], σ_ω=exp(x[11]), π_ω, Π_ω)
# block = [3,4,5,6,7,9,11,12]
# moms1 = get_moments(θ,V,F,M,L)
# [moms0 moms1]
#ssq(θ,V,F,moms0,M,L)

#wght = 1 ./ (moms0 .* (1 .- moms0))
θ = update(res1.minimizer,θ,F)

# now go after the unweighted version:
res2 = optimize(x->ssq(update(x,θ,F),V,F,moms0,M,L),res1.minimizer,Optim.Options(iterations=1000))
#res2 = optimize(x->ssq(update(x,θ,F),V,F,moms0,M,L),res2.minimizer,Optim.Options(iterations=200))

θ = update(res2.minimizer,θ,F)

m1 = get_moments(θ,V,F,M,L)
display([moms0 m1])

break

idx = LinearIndices((F.N_ϵ,F.N_ω,F.N_d,2,2,2));

x1 = res1.minimizer;
# doesn't seem to work
θ2 = (;θ...,);
θ2.α_ω[1] = 0.
θ2.α_νH0[1] = 3.88
θ2 = (;θ2...,α_F = -40. )
θ2 = (;θ2..., π_ω = 0.9, Π_ω = transmat_ω(0.9,F.N_ω) )

m1 = get_moments(θ2,V,F,M,L)
display([moms0 m1])



# α_ω,σ_ω,α_F,σ_F,α_νH[0],α_νH1[1], π_ω, (too many parameters)
# divorce bfore kids, divorce eventually, fertility eventually


# # doesn't work that well, try the simultaneous version?
# blocks = [1:2,3:4,5:6,7:11]
# for b in blocks
#     xcurrent = get_x(θ)
#     x0 = xcurrent[b]
#     res = optimize(x->ssq(update(x,b,θ,F),V,F,moms0,M,L),x0,Optim.Options(iterations=50))
#     xcurrent[b] .= res.minimizer
#     θ = update(xcurrent,θ,F)
# end


