# load everything
# load the estimates
include("../src/model/model.jl")
include("../src/estimation/estimation.jl")
include("../src/counterfactuals.jl")
using DelimitedFiles

M = CSV.read("data/MarriageFile.csv",DataFrame,missingstring="NA")
P = CSV.read("data/MarriagePanel.csv",DataFrame,missingstring="NA")
K = CSV.read("data/KidPanel.csv",DataFrame,missingstring="NA")
cprobs = CSV.read("data/CustodyMomentsSimple.csv",DataFrame).p[:]

F = FixedParams()
θ = Params(F)
V = values(F);
θ = (;θ...,cprobs)
θk = prod_pars()
θse = Params(F)
θkse = prod_pars()

x1 = readdlm("output/est_stage1")[:]
x2 = readdlm("output/est_stage2")[:]
x3 = readdlm("output/est_stage3")[:]
x4 = readdlm("output/est_stage4")[:]
x5 = readdlm("output/est_stage5")[:]

θ,θk = update_all((x1,x2,x3,x4,x5),θ,θk,F);

X1b = readdlm("output/boot_stage1")
X2b = readdlm("output/boot_stage2")
X3b = readdlm("output/boot_stage3")
X4b = readdlm("output/boot_stage4")
X5b = readdlm("output/boot_stage5")


θk = (;θk...,ρ = θ.ρ)

# ----- Simulate Data for model fit
dat = prep_sim_data(M,P;R = 10)
model = (;θ,values=V,F)
solve_all!(model)
sim_data = data_gen(model,dat)
kid_data = prep_child_data(sim_data,dat,cprobs);

# get baseline stats
stats_baseline = counterfactual_statistics(kid_data,dat,θ,θk,model)

r1,r2 = divorce_standard_counterfactual(dat,model,θk,stats_baseline)
# do we have to reset the divorce standard after this? I think so.

r1,r2 = custody_counterfactual(dat,model,θk,stats_baseline)


r = child_support_counterfactual(dat,model,θk,stats_baseline)

break

dsim = gen_data_frame(model,dat)

using StatsPlots

# shows divorce policy, unlikely for couples without kids
@chain dsim begin
    @subset :AK.<18
    @transform :kids = :AK.>=0
    groupby([:kids,:omega])
    @combine :pD = mean(:pD)
    #@transform :eW = Int64.(:eW)
    @df plot(:omega,:pD,group=:kids)#,color = :kids)
end

# divorce probability with experience
@chain dsim begin
    @subset :AK.<18
    @transform :kids = :AK.>=0
    groupby([:kids,:exp])
    @combine :pD = mean(:pD)
    #@transform :eW = Int64.(:eW)
    @df plot(:exp,:pD,group=:kids)#,color = :kids)
end

# divorce probability with age (shows that all action without kids is to do with age)
@chain dsim begin
    @subset :AK.<18 :AW.<50
    @transform :kids = :AK.>=0
    groupby([:kids,:AW])
    @combine :pD = mean(:pD)
    #@transform :eW = Int64.(:eW)
    @df plot(:AW,:pD,group=:kids)#,color = :kids)
end


@chain dsim begin
    @subset :AK.<0 :AW.<35 #:pD.>0.1
    #@combine :m = maximum(:pD)
end


@chain dsim begin
    groupby([:eW,:omega])
    @combine :FT = mean(:FT)
    @transform :eW = Int64.(:eW)
    @df plot(:omega,:FT,group=:eW,color = :eW)
end

@chain dsim begin
    @subset :AK.<18
    @transform :kids = :AK.>=0
    groupby([:kids,:omega])
    @combine :FT = mean(:FT)
    #@transform :eW = Int64.(:eW)
    @df plot(:omega,:FT,group=:kids)#,color = :kids)
end



# shows an increase in labor supply approach divorce
@chain dsim begin
    @subset :tsd.>-5 :tsd.<5
    groupby([:eW,:tsd])
    @combine :FT = mean(:FT)
    @transform :eW = Int64.(:eW)
    @df plot(:tsd,:FT,group=:eW,color = :eW)
end

# same figure by fertility status instead of education
@chain dsim begin
    @transform :kids = :AK.>=0
    @subset :tsd.>-5 :tsd.<5
    groupby([:kids,:tsd])
    @combine :FT = mean(:FT)
    #@transform :eW = Int64.(:eW)
    @df plot(:tsd,:FT,group=:kids)#,color = :kids)
end


# shows a child penalty in full-time work
@chain dsim begin
    @subset :tsf.>-8 :tsf.<10
    groupby([:eW,:tsf])
    @combine :FT = mean(:FT)
    @transform :eW = Int64.(:eW)
    @df plot(:tsf,:FT,group=:eW,color = :eW)
end

@chain dsim begin
    @subset :tsf.>-8 :tsf.<10
    groupby([:eW,:eH,:tsf])
    @combine :FT = mean(:FT)
    @transform :eW = Int64.((:eW .- 1)*2 .+ :eH)
    @df plot(:tsf,:FT,group=:eW)# ,color = :eW)
end


function model_stats(mod,dat;seed=1234)
    (;values,θ,F) = mod
    (;Π_ϵ,Π_ω) = θ
    (;N_ϵ,N_t,N_ω,A_bar,N_d,T_f) = F
    (;VW1,VH1) = values
    (;γ_YW) = θ

    # = MersenneTwister(seed)
    Random.seed!(seed)

    # interpolate choice probabilities
    pL1,pL2,pL3,pL4,pL5,pD1,pD2,pD3,pF = interpolate_probs(values,F)
    
    # interpolate stage 1 values
    etpW = [interpolate3(VW1,F,t) for t in 1:24]
    etpH = [interpolate3(VH1,F,t) for t in 1:24]

    # create an indexing rule
    idx = LinearIndices((N_ϵ,N_ω,N_d,2,2,2))

    # create distribution objects (including initial)
    πϵ = [Categorical(Π_ϵ[ϵi,:]) for ϵi=1:N_ϵ]
    πϵ0 = eigvecs(I(N_ϵ) .- Π_ϵ')[:,1]
    πϵ0 ./= sum(πϵ0)
    πϵ_init = Categorical(πϵ0)
    πω = [Categorical(Π_ω[ωi,:]) for ωi=1:N_ω]
    πω0 = Categorical(fill(1/N_ω,N_ω))

    # 
    N = size(dat.edW,1)

    # welfare
    X = zeros(N,7)

    nt = 1
    for n in 1:N
        eW = dat.edW[n]
        eH = dat.edW[n] #!!!!!!!
        ad = dat.AD[n]
        adi = 1 + (ad>-5) + (ad>5)
        aw0 = dat.AW0[n]
        maxT = dat.tlength[n]
        κ = dat.κ[n]
        #ωt = N_ω #<- everyone starts at the highest draw (an alternative assumption)
        ωt = rand(πω0) 
        ϵt = rand(πϵ_init) 
        tM = 1
        ttF = 9999
        ttD = 9999
        stage = 1
        AK = -1

        i2 = idx[ϵt,ωt,adi,eH,eW,2]
        i1 = idx[ϵt,ωt,adi,eH,eW,1]
        tt = aw0+1-18
        X[n,1] = etpH[tt](κ,i2) - etpH[tt](κ,i1)
        X[n,2] = etpW[tt](κ,i2) - etpW[tt](κ,i1)
        X[n,3] = κ
        X[n,4] = ωt
        X[n,5] = ϵt
        X[n,6] = eW
        X[n,7] = eH
        nt += maxT
    end
    d = DataFrame(dH = X[:,1],dW = X[:,2],kappa = X[:,3],omega = X[:,4],epsilon = X[:,5], edH = X[:,6],edW = X[:,7])
    return d
end

d = model_stats(model,dat)

mod = lm(@formula(dH ~ omega * (edW==2)),d)

mod = lm(@formula(dH ~ omega + epsilon + (edW==2) + (edH==2)),d)

using DataFramesMeta

@chain d begin
    groupby([:edW,:omega])
    @combine :dm = mean(:dH)
    @transform :edW = Int64.(:edW)
    @df plot(:omega,:dm,color = :edW)
end

@chain d begin
    groupby([:epsilon,:omega])
    @combine :dm = mean(:dH)
    @transform :epsilon = Int64.(:epsilon)
    @df plot(:omega,:dm,color = :epsilon)
end

