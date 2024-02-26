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

