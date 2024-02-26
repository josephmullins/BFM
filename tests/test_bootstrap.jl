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

x1 = readdlm("output/est_stage1")[:]
x2 = readdlm("output/est_stage2")[:]
x3 = readdlm("output/est_stage3")[:]
x4 = readdlm("output/est_stage4")[:]
x5 = readdlm("output/est_stage5")[:]

θ,θk = update_all((x1,x2,x3,x4,x5),θ,θk,F);

# practice getting a bootstrap sample
# require the K and P be ordered properly
function draw_boot_sample(M,P,K)
    N = nrow(M)
    Mb = M[rand(1:N,N),:]
    Mb[!,:MIDb] = 1:N

    Pb = @chain Mb begin
        @select :MID :MIDb
        innerjoin(P,on=:MID)
        @orderby :MIDb :YEAR
    end
    Kb = @chain Mb begin
        @select :MID :MIDb
        innerjoin(K,on=:MID)
        @orderby :MIDb :KID :YEAR
    end
    return Mb,Pb,Kb
end

#OK: problem (?) do we want to factor simulation noise into this? yes, we do?
Mb,Pb,Kb = draw_boot_sample(M,P,K)

x4_0 = get_x(θ)
x5_0 = get_xk(θk)

θ,θk = estimate_model(x4_0, x5_0, θ, θk, V, F, Mb, Pb, Kb ; R = 10, num_iter = 5, show_trace = true)

x4_0 = readdlm("output/Xsave_round2")[:,1]
x4_1 = readdlm("output/Xsave_round2")[:,2]

dat = prep_sim_data(M,P;R = 10)
moms0 = data_moms(M,P)
ssq(update(x4_0,θ,F),V,F,moms0,dat)
ssq(update(x4_1,θ,F),V,F,moms0,dat)

θ = update(x4_0,θ,F)
moms1 = get_moments(θ,V,F,dat)
