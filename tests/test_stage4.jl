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
mod = (;θ,values=V,F)

data_gen(mod,dat);
@time data_gen(mod,dat)

break


wght = [[1.,0.,1.,0.];ones(10)]

res1 = optimize(x->ssq(update(x,θ,F),wght,V,F,moms0,dat,legal),x0,Optim.Options(iterations=200,show_trace=true))
θ = update(res1.minimizer,θ,F)
m1 = get_moments(θ,V,F,dat,legal)
display([moms0 m1])

(;Π_ϵ,Π_ω) = θ
(;N_ϵ,N_t,N_ω,A_bar,N_d,T_f) = F
# = MersenneTwister(seed)
#Random.seed!(seed)
# interpolate choice probabilities
policies = interpolate_probs(V,F);
idx = LinearIndices((N_ϵ,N_ω,N_d,2,2,2))
πϵ = [Categorical(Π_ϵ[ϵi,:]) for ϵi=1:N_ϵ]
πϵ0 = eigvecs(I(N_ϵ) .- Π_ϵ')[:,1]
πϵ0 ./= sum(πϵ0)
πϵ_init = Categorical(πϵ0)

πω = [Categorical(Π_ω[ωi,:]) for ωi=1:N_ω]
πω0 = Categorical(fill(1/N_ω,N_ω))
distributions = (πϵ,πϵ_init,πω,πω0);

N = size(dat,1)
TD = zeros(Int64,N)
TF = zeros(Int64,N)
NT = size(legal,1)
L_sim = zeros(NT)
Ω_sim = zeros(NT)
D = zeros(Bool,NT)
nt = 1

d = (;edW = dat.edW,edH = dat.edH,tlength = dat.tlength,YBIRTH_M = dat.YBIRTH_M,YBIRTH_F = dat.YBIRTH_F,YMAR = dat.YMAR,exp0 = dat.exp0);

@code_warntype sim_panel!(TD,TF,L_sim,1,nt,F,policies,distributions,idx,d,legal)

@code_warntype data_gen(mod,d,legal)

break

# blocks = [
#     [5,6,7,9,11,12],
#     [3,4,5,7],
#     [3,5,7,8,9,10,11],
#     [1,2]
# ]
# θ = estimate_blocks(θ,blocks,wght,V,F,moms0,M,L)
# θ = estimate_blocks(θ,blocks,wght,V,F,moms0,M,L)


x0 = get_x(θ)

wght = ones(14); wght[11] = 50.
res1 = optimize(x->ssq(update(x,θ,F),wght,V,F,moms0,M,L),x0,Optim.Options(iterations=200,show_trace=true))
res1 = optimize(x->ssq(update(x,θ,F),wght,V,F,moms0,M,L),res1.minimizer,Optim.Options(iterations=1000,show_trace=true))
θ = update(res1.minimizer,θ,F)
m1 = get_moments(θ,V,F,M,L)
#display([moms0 m1])

res2 = optimize(x->ssq(update(x,θ,F),V,F,moms0,M,L),x0,Optim.Options(iterations=200,show_trace=true))
res2 = optimize(x->ssq(update(x,θ,F),V,F,moms0,M,L),res2.minimizer,Optim.Options(iterations=1000,show_trace=true))
θ = update(res2.minimizer,θ,F)
m2 = get_moments(θ,V,F,M,L)

mstring = ["FT - married","PT - married","FT - divorced","PT - divorced"]
mstring = [mstring;["divorce < $x years" for x in (5,10,15)]]
mstring = [mstring;["birth < $x years" for x in (2,4,6)]]
mstring = [mstring;["time to divorce - time to birth < $x years" for x in (0,5,10,15)]]
display([mstring moms0 m1 m2])

