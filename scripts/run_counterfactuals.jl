# load everything
# load the estimates
include("../src/model/model.jl")
include("../src/estimation/estimation.jl")
include("../src/counterfactuals.jl")
using DelimitedFiles

M = CSV.read("data/MarriageFile.csv",DataFrame,missingstring="NA")
P = CSV.read("data/MarriagePanel.csv",DataFrame,missingstring="NA")
K = CSV.read("data/KidPanelv2.csv",DataFrame,missingstring="NA")
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

# ----- Simulate Data  for baseline
dat = prep_sim_data(M,P;R = 10)
model = (;θ,values=V,F)
solve_all!(model)
sim_data = data_gen(model,dat)
kid_data = prep_child_data(sim_data,dat,cprobs);

# get baseline stats
stats_baseline = counterfactual_statistics(kid_data,dat,θ,θk,model)

# get outcomes under mutual consent and unilateral / no-fault divorce
mc,unil = divorce_standard_counterfactual(dat,model,θk,stats_baseline)

# get outcomes for full-time maternal custody and part-time
mcust,ptcust = custody_counterfactual(dat,model,θk,stats_baseline)

# get outoutcomes for π=0.3 and π=0.5
childsupp1,childsupp2 = child_support_counterfactual(dat,model,θk,stats_baseline)
mc_boot,unil_boot,mcust_boot,pcust_boot,cs1_boot,cs2_boot = bootstrap_counterfactuals(X1b,X2b,X3b,X4b,X5b,model,M,P,K,θk)

# write these results to a table
write_counterfactuals_table((mc,unil),(mc_boot,unil_boot),"divorce_standard.tex")
write_counterfactuals_table((mcust,ptcust),(mcust_boot,pcust_boot),"custody_standard.tex")
write_counterfactuals_table((childsupp1,childsupp2),(cs1_boot,cs2_boot),"child_support.tex")

# now get outcomes when everyone is mututal consent and unilateral
stats_mc,stats_ul = divorce_standard_comparison(dat,model,θk)
div_standard = compare_stats(stats_ul,stats_mc,θ.β)

# Now re-run the π=0.3 counterfactual under both regimes
childsupp1,childsupp2 = extended_child_support_counterfactual(dat,model,θk,stats_mc,stats_ul)

ds_boot,cs1_boot,cs2_boot = bootstrap_extra_counterfactuals(X1b,X2b,X3b,X4b,X5b,model,M,P,K,θk)

# write these results to a table as well
write_counterfactuals_table((div_standard,),(ds_boot,),"divorce_standard_comparison.tex")
write_counterfactuals_table((childsupp1,childsupp2),(cs1_boot,cs2_boot),"extended_child_support.tex")

# save simulation results for figures in the paper
div = [b.Δdiv for b in ds_boot]
fert = [b.Δfert for b in ds_boot]
d = [DataFrame(x = div, y = fert, name = "Change in Fertility Rate (%)") ;
    DataFrame(x = div, y = X4b[11,:], name = "Dispersion of Marriage Taste Shocks")]
CSV.write("output/dissolution.csv",d)
