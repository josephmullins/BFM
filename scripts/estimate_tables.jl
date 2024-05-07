# load everything
# load the estimates
include("../src/model/model.jl")
include("../src/estimation/estimation.jl")
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

se1 = std(X1b,dims=2)[:]
se2 = std(X2b,dims=2)[:]
se3 = std(X3b,dims=2)[:]
se4 = std(X4b,dims=2)[:]
se5 = std(X5b,dims=2)[:]

# calculate standard errors and make a new θ
θse,θkse = update_all((se1,se2,se3,se4,se5),θse,θkse,F);

θk = (;θk...,ρ = θ.ρ)
θkse = (;θkse...,ρ = θse.ρ)

# ---- tables with estimates
write_inc_table(θ,θse)
write_prefs(θ,θse)
write_production(θk,θkse)

# ---- save estimates to file for making figures
d = [DataFrame(Age = 0:17,Input = "Mother's Time", value = θk.δW, se = θkse.δW) ; 
    DataFrame(Age = 0:17,Input = "Father's Time", value = θk.δH, se = θkse.δH) ]
CSV.write("output/factor_shares.csv",d)

tfp_boot = X5b[3,:] .+ X5b[4,:].*F.ω_grid'
se_tfp = std(tfp_boot,dims=1)[:]
d = [DataFrame(MarriageQuality = F.ω_grid,case = "Married",TFP = θk.γ_ψ[2] .+ θk.γ_ψ[3]*F.ω_grid, se = se_tfp) ;
DataFrame(MarriageQuality = F.ω_grid,case = "Divorced",TFP = θk.γ_ψ[1], se = θkse.γ_ψ[1])]
CSV.write("output/tfp.csv",d)

tfp_boot = X5b[3,:] .+ X5b[4,:].*F.ω_grid' .- X5b[2,:]
se_tfp = std(tfp_boot,dims=1)[:]
d = DataFrame(MarriageQuality = F.ω_grid,TFP = θk.γ_ψ[2] .+ θk.γ_ψ[3]*F.ω_grid .- θk.γ_ψ[1], se = se_tfp)
CSV.write("output/rel_tfp.csv",d)


# ----- Simulate Data for model fit
dat = prep_sim_data(M,P;R = 10)
mod = (;θ,values=V,F)
solve_all!(mod)
sim_data = data_gen(mod,dat)
kid_data = prep_child_data(sim_data,dat,cprobs);
S = zeros(length(kid_data.ΩK))

# ------ Calculate moments
moms0 = data_moms(M,P)
moms1 = get_moments(θ,V,F,dat)
kmoms0 = kidmoms_data(K)
kmoms1 = kid_moments(S,θk,θ,F,kid_data)

# ------ Write tables and data frames
# write tables with model fit from stage 4 and stage 5
write_modelfit_s4(moms0,moms1)
write_modelfit_s5(kmoms0[3],kmoms1[3])

# save a dataframe with the test score fit (and sd fit)
d = [DataFrame(Age = repeat(4:15,3),Dgroup=repeat(1:3,inner=12),value=kmoms0[1],case="Data") ;
DataFrame(Age = repeat(4:15,3),Dgroup=repeat(1:3,inner=12),value=kmoms1[1],case="Model")]
CSV.write("output/modelfit_testcores.csv",d)

rel1 = [kmoms0[1][13:24] .- kmoms0[1][1:12] ; kmoms0[1][25:36] .- kmoms0[1][1:12]]
rel2 = [kmoms1[1][13:24] .- kmoms1[1][1:12] ; kmoms1[1][25:36] .- kmoms1[1][1:12]]
d = [DataFrame(Age = repeat(4:15,2),Dgroup=repeat(2:3,inner=12),value= rel1,case="Data") ;
DataFrame(Age = repeat(4:15,2),Dgroup=repeat(2:3,inner=12),value=rel2,case="Model")]
CSV.write("output/modelfit_testcores_relative.csv",d)


CSV.write("output/modelfit_sd.csv",[DataFrame(Age = 4:15,sd = kmoms0[2],case="Data") ; DataFrame(Age = 4:15,sd = kmoms1[2],case="Model")])

# calculate the input decomposition
function divorce_comparison(θk,θ,F,kid_data)
    inputs,mstat = input_decomposition(θk,θ,F,kid_data)
    d = [mean(inputs[:,mstat.==0],dims=2) mean(inputs[:,mstat.==1],dims=2)]
    d_sum = sum(d,dims=1)
    sd_skills = std(sum(inputs,dims=1))
    d_gap = (d_sum[1] - d_sum[2]) / sd_skills
    d ./= sd_skills
    return [d_gap;d[:,1] .- d[:,2]]
end

d = divorce_comparison(θk,θ,F,kid_data)

function bootstrap_skilldecomp(Xb,mod,θk)
    (;θ,F) = mod
    X1b,X2b,X3b,X4b,X5b = Xb
    B = size(X1b,2)
    Db = zeros(4,B)
    θ,θk = update_all((X1b[:,1],X2b[:,1],X3b[:,1],X4b[:,1],X5b[:,1]),θ,θk,F);
    for b in axes(X1b,2)
        println(b)
        θ,θk = update_all((X1b[:,b],X2b[:,b],X3b[:,b],X4b[:,b],X5b[:,b]),θ,θk,F);
        mod = (;mod...,θ)
        sim_data,kid_data = full_simulation(dat,mod,θ.cprobs)
        Db[:,b] = divorce_comparison(θk,θ,F,kid_data)     
    end
    #return std(Db,dims=2)
    #return [(quantile(Db[s,:],0.025),quantile(Db[s,:],0.975)) for s in axes(Db,1)]
    return Db
end

dse = bootstrap_skilldecomp((X1b,X2b,X3b,X4b,X5b),mod,θk)

ci = [(quantile(dse[s,:],0.05),quantile(dse[s,:],0.95)) for s in axes(dse,1)]
write_decomposition(d,ci)

# calculate outcomes for bilateral vs unilateral (easy, just through L)

# calculate outcomes for child support