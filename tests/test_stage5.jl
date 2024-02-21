
include("../src/model/model.jl")
include("../src/estimation/estimation.jl")
using DelimitedFiles

M = CSV.read("data/MarriageFile.csv",DataFrame,missingstring="NA")
P = CSV.read("data/MarriagePanel.csv",DataFrame,missingstring="NA")
K = CSV.read("data/KidPanel.csv",DataFrame,missingstring="NA")
cprobs = CSV.read("data/CustodyMomentsSimple.csv",DataFrame).p[:]

K[!,:dgroup] .= 1 .+ (K.YDIV.<9000) .+ (K.TSD.>=0) #<- THIS NEEDS TO GO WITH THE DATA PREP!!!


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

mod = (;θ,values=V,F)

solve_all!(mod)
sim_data = data_gen(mod,dat);

kid_data = prep_child_data(sim_data,dat,cprobs);

S = zeros(length(kid_data.ΩK))
θk = (;γ_ψ0 = -5,
    γ_ψ = [0.,0.,0.2,0.],
    δW = fill(0.1,18),δH = fill(0.1,18),δk = 0.9,
    Γa = zeros(19))

x0 = [-5.,0.,0.,0.2,0.,1.,log(0.1),log(0.1)]
θk =  updateθk(x0,θk,θ)

# the data moments:
k = kidmoms_data(K)
kmoms0 = k.S

# we can do the same thing for the data easily enough
obj_stage5(S,θk,θ,F,kid_data,kmoms0)

res = optimize(x->obj_stage5(S,updateθk(x,θk,θ),θ,F,kid_data,kmoms0),x0,Optim.Options(show_trace = true))

θk = updateθk(res.minimizer,θk,θ)
kmoms1 = kid_moments(S,θk,θ,F,kid_data)

# now we just need the update function and we're ready to go
k[!,:Smod] .= kmoms1

using StatsPlots
@df stack(k) plot(
    :AGE,
    :value,
    group = (:dgroup,:variable),
    linecolor = :dgroup
)

p = plot()
K0 = reshape(kmoms0,15,3)
K1 = reshape(kmoms1,15,3)
clrs = ["blue","red"]
for i in 2:3
    plot!(p,K0[:,i] .- K0[:,1],linecolor=clrs[i-1])
    plot!(p,K1[:,i] .- K1[:,1],linecolor=clrs[i-1],linestyle=:dash)
end
p


