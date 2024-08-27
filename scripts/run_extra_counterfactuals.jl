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

# ----- Simulate Data for model fit
dat = prep_sim_data(M,P;R = 10)
model = (;θ,values=V,F)
solve_all!(model)
sim_data = data_gen(model,dat)
kid_data = prep_child_data(sim_data,dat,cprobs);

# get baseline stats
stats_baseline = counterfactual_statistics(kid_data,dat,θ,θk,model)

stats_mc,stats_ul = divorce_standard_comparison(dat,model,θk)
div_standard = compare_stats(stats_ul,stats_mc,θ.β)

childsupp1,childsupp2 = extended_child_support_counterfactual(dat,model,θk,stats_mc,stats_ul)

ds_boot,cs1_boot,cs2_boot = bootstrap_extra_counterfactuals(X1b,X2b,X3b,X4b,X5b,model,M,P,K,θk)

# a function for writing everything to a table

function write_counterfactuals_table(results,boot,filename)
    stat_names = ["\\% CEV Husbands","\\% CEV Wives","\$\\Delta\$ Wife log-wage (\$\\times\$ 100)","\$\\Delta\$ Fertility (\\%)","\$\\Delta\$ Divorce (\\%)"]
    field_names = [:ΔwH,:ΔwW,:Δlogwage,:Δfert,:Δdiv,:Δskill]
    file = open("output/tables/"*filename,"w")

    for r in eachindex(stat_names)
        write(file,stat_names[r]) #form(d[r])," & ",formci(dse[r])," \\\\ \n")
        for cr in results
            num = getfield(cr,field_names[r])
            write(file," & ",form(num))
        end
        write(file," \\\\ \n")
        for cb in boot
            lb = quantile([getfield(b,field_names[r]) for b in cb],0.05)
            ub = quantile([getfield(b,field_names[r]) for b in cb],0.95)
            write(file," & ",formci((lb,ub)))
        end
        write(file," \\\\ \n")
    end

    stat_names = ["\$\\Delta\$ Skill (\\% sd)","\\hspace{10pt}\$\\Delta\$ TFP","\\hspace{10pt}\$\\Delta\$ Mother's Time","\\hspace{10pt}\$\\Delta\$ Father's Time"]
    for r in eachindex(stat_names)
        write(file,stat_names[r]) #form(d[r])," & ",formci(dse[r])," \\\\ \n")
        for cr in results
            num = getfield(cr,:Δskill)[r]
            write(file," & ",form(num))
        end
        write(file," \\\\ \n")
        for cb in boot
            lb = quantile([getfield(b,:Δskill)[r] for b in cb],0.05)
            ub = quantile([getfield(b,:Δskill)[r] for b in cb],0.95)
            write(file," & ",formci((lb,ub)))
        end
        write(file," \\\\ \n")
    end

    close(file)
end

write_counterfactuals_table((div_standard,),(ds_boot,),"divorce_standard_comparison.tex")
write_counterfactuals_table((childsupp1,childsupp2),(cs1_boot,cs2_boot),"extended_child_support.tex")

div = [b.Δdiv for b in ds_boot]
fert = [b.Δfert for b in ds_boot]
d = [DataFrame(x = div, y = fert, name = "Change in Fertility Rate (%)") ;
    DataFrame(x = div, y = X4b[11,:], name = "Dispersion of Marriage Taste Shocks")]
CSV.write("output/dissolution.csv",d)

