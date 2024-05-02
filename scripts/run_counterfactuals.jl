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

mc,unil = divorce_standard_counterfactual(dat,model,θk,stats_baseline)
# do we have to reset the divorce standard after this? I think so.

mcust,ptcust = custody_counterfactual(dat,model,θk,stats_baseline)


childsupp1,childsupp2 = child_support_counterfactual(dat,model,θk,stats_baseline)



break

function test_func()
    return (;a = rand(), b = rand(2))
end

a = [test_func(),test_func()]

push!(a,test_func())
function test_boot()
    R = []
    for b in 1:10
        r = test_func()
        push!(R,r)
    end
    R
end


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
