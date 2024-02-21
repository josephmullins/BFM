
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

L = convert(Vector{Int64},P.L)
moms0 = data_moms(M,P)

dat = prep_sim_data(M,P;R = 10)

θ = (;θ...,
    σ_L = 2., α_l = 1.,
    α_νd = 0.,σ_F = 10.,
    α_ω = [0., 5.],σ_ω = 3., # σ_ω = 2.
    π_ω = 0.7, Π_ω = transmat_ω(0.7,F.N_ω),
    α_νH0 = [3.,-0.1],α_νH1 = [7.,0.2])

x0 = get_x(θ)
moms1 = get_moments(θ,V,F,dat)

mod = (;θ,values=V,F)

solve_all!(mod)
sim_data = data_gen(mod,dat);

kid_data = prep_child_data(sim_data,dat);

