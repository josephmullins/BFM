
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

mod = (;θ,values=V,F);

solve_all!(mod)
sim_data = data_gen(mod,dat);
TD = sim_data[1];
L_sim = sim_data[3];
Ω_sim = sim_data[4];
tsd = Int64[]
for n in eachindex(TD)
    tlength = dat.tlength[n]
    for t in 1:tlength
        push!(tsd,t-1-TD[n])
    end
end

pt = [mean(L_sim[tsd.==x].==1) for x in -6:10]
ft = [mean(L_sim[tsd.==x].==2) for x in -6:10]
L = [mean(L_sim[tsd.==x].>0) for x in -6:10]
plot(-6:10,ft)
plot!(-6:10,pt)
plot!(-6:10,L)

ft = [mean(L_sim[Ω_sim.==x].==2) for x in 1:5]
pt = [mean(L_sim[Ω_sim.==x].==1) for x in 1:5]
plot(1:5,ft)

plot(-6:0,[mean(Ω_sim[tsd.==x].==g) for x in -6:0, g in 1:5])
plot(-6:0,[mean(Ω_sim[tsd.==x]) for x in -6:0])
plot(-10:10,[mean(L_sim[tsd.==x].==g) for x in -10:10, g in 1:2])

idx = LinearIndices((N_ϵ,N_ω,N_d,2,2,2))
pL1,pL2,pL3,pL4,pL5,pD1,pD2,pD3,pF = interpolate_probs(V,F);


kid_data = prep_child_data(sim_data,dat,cprobs);

S = zeros(length(kid_data.ΩK))
θk = (;γ_ψ0 = -10,
    γ_ψ = [0.,0.,0.5,0.],
    δW = fill(0.1,18),δH = fill(0.1,18),δk = 0.9,
    Γa = zeros(19))

x0 = [-5.,0.,0.,0.2,0.,1.,log(0.4)]
θk =  updateθk(x0,θk,θ)

# the data moments:
kmoms0 = kidmoms_data(K)
kmoms1 = kid_moments(S,θk,θ,F,kid_data)

wght = ones(41)
wght[40:41] .= 1000.
# we can do the same thing for the data easily enough
obj_stage5(S,θk,θ,F,kid_data,kmoms0,wght)

res = optimize(x->obj_stage5(S,updateθk(x,θk,θ),θ,F,kid_data,kmoms0,wght),x0,Optim.Options(show_trace = true))





# ---- a cheeky plot to show what's what
θk = updateθk(res.minimizer,θk,θ)
kmoms1 = kid_moments(S,θk,θ,F,kid_data)

using StatsPlots
p = plot()
K0 = reshape(kmoms0[1:end-2],13,3)
K1 = reshape(kmoms1[1:end-2],13,3)
clrs = ["blue","red"]
for i in 2:3
    plot!(p,K0[:,i] .- K0[:,1],linecolor=clrs[i-1])
    plot!(p,K1[:,i] .- K1[:,1],linecolor=clrs[i-1],linestyle=:dash)
end
p


x0 = [-7.5,0.3,0.,0.05,0.00,2.5,log(0.2)]
θk =  updateθk(x0,θk,θ)
kmoms1 = kid_moments(S,θk,θ,F,kid_data)

p = plot()
K0 = reshape(kmoms0[1:end-2],13,3)
K1 = reshape(kmoms1[1:end-2],13,3)
clrs = ["blue","red","purple"]
for i in 1:3
    plot!(p,K0[:,i],linecolor=clrs[i])
    plot!(p,K1[:,i],linecolor=clrs[i],linestyle=:dash)
end
p


m = @chain K begin
    @subset :AGE.<=15
    @subset :AGE.>=3
    #@transform :AGE = round.(:AGE ./ 4)
    groupby([:dgroup,:AGE])
    @combine :S = mean(skipmissing(:AP_std))
end
plot(reshape(m.S,13,3))

m = @chain K begin
    @subset :AGE.<=15
    @subset :AGE.>=3
    @select :AGE :dgroup :AP_raw
    dropmissing()
    groupby(:AGE)
    @transform :S = (:AP_raw .- mean(:AP_raw)) ./ std(:AP_raw)
    @transform :AGE = round.(:AGE ./ 4)
    groupby([:dgroup,:AGE])
    @combine :S = mean(skipmissing(:S))
end
plot(reshape(m.S,4,3))
plot(reshape(m.S,13,3))

m = @chain K begin
    @subset :AGE.<=15
    @subset :AGE.>=3
    @select :AGE :dgroup :LW_raw
    @transform :AGE = round.(:AGE ./ 4)
    groupby([:dgroup,:AGE])
    @combine :S = mean(skipmissing(:S))
end
plot(reshape(m.S,4,3))