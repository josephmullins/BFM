# this script re-estimates stage 5 (production parameters fixing everything else)
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

function rerun_stage5(x,x5_0,mod,θk,data,panel,kid_data; show_trace=true, seed = 20240220)
    (x1,x2,x3,x4,x5) = x
    (;θ,F) = mod
    θ,θk = update_all((x1,x2,x3,x4,x5),θ,θk,F)
    mod = (;mod...,θ)
    solve_all!(mod)
    dat = prep_sim_data(data,panel;R = 10)
    θk = stage5(x5_0, θk, mod, dat, kid_data; num_iter=500, show_trace, seed)
    return θk
end

# get new estimates
model = (;θ,values=V,F)

x5_0 = [-15.,0.75,0.75,0.05,0.00,2.5,log(0.2),log(0.7),log(12.)]

θk = rerun_stage5((x1,x2,x3,x4,x5),x5_0,model,θk,M,P,K)
x1,x2,x3,x4,x5 = stack_ests(θ,θk)
x5_1 = get_xk(θk) #<- for the remaining bootstrap trials

writedlm("output/est_stage5",x5)

for b in axes(X5b,2)
    @show "doing trial $b"
    Random.seed!(1010+b)
    Mb,Pb,Kb = draw_boot_sample(M,P,K)
    θk = rerun_stage5((X1b[:,b],X2b[:,b],X3b[:,b],X4b[:,b],X5b[:,b]),x5_1,model,θk,Mb,Pb,Kb ; show_trace = false, seed = 3030+b)
    _,_,_,_,X5b[:,b] = stack_ests(θ,θk)
end

writedlm("output/boot_stage5",X5b)
