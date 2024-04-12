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

Mb = zeros()

function time_moms(d)
    phi_m = @chain d begin
        @subset :DIV.==0 :AGE.<=17 :AGE.>=7
        groupby(:AGE)
        @combine :PHIm = mean(skipmissing(:tau_f)) :wght = sum(.!ismissing.(:tau_f))
    end
    phi_d = @chain d begin
        @subset :DIV.==1 :AGE.<=17 :AGE.>=7
        groupby(:AGE)
        @combine :PHId = mean(skipmissing(:tau_f)) #:wght = sum(.!ismissing.(:tau_f))
    end
    m = @chain phi_d begin
        innerjoin(phi_m,on=:AGE)
        @transform :r = :PHId ./ :PHIm
        @combine :m=sum(:wght .* :r) / sum(:wght)
        _.m[1]
    end
    return m
end

function time_moms_alt(d)
    phi_m = @chain d begin
        @subset :DIV.==0 :AGE.<=17 #:AGE.>=7
        @combine :m = mean(skipmissing(:tau_f))
    end
    phi_d = @chain d begin
        @subset :DIV.==1 :AGE.<=17 #:AGE.>=7
        @combine :m = mean(skipmissing(:tau_f))
    end
    return phi_d.m[1] / phi_m.m[1]
end


m0 = time_moms(K)

break

num_boot = 100
mb = zeros(num_boot)
#mb = zeros(11,num_boot)
ρb = zeros(num_boot)
ρb2 = zeros(num_boot)

modratio(ρ,τgrid,p) = sum((1 ./ (1 .+ ρ.*(1 .- τgrid))).*p)
modratio2(ρ,τgrid,p) = sum((τgrid ./ (1 .+ ρ.*(1 .- τgrid))).*p)

modratio3(ρ,τgrid,p) = sum(((τgrid .+ τgrid.^2) ./ (1 .+ ρ.*(1 .- τgrid))).*p)



for b in 1:num_boot
    println("Performing bootstrap trial $b")
    Random.seed!(1010 + b)
    Mb,Pb,Kb = draw_boot_sample(M,P,K)
    mb[b] = time_moms_alt(Kb)
    mb2 = time_moms(Kb)
    #mb[:,b] = time_moms(Kb)
    #θ = stage3(θ, Kb, F.τgrid, θ.cprobs)
    #ρb[b] = θ.ρ
    res = optimize(ρ->(mb2 - modratio(ρ,F.τgrid,θ.cprobs))^2,0,20)
    ρb[b] = res.minimizer
    res = optimize(ρ->(mb2 - modratio3(ρ,F.τgrid,θ.cprobs))^2,0,20)
    ρb2[b] = res.minimizer
    
end

#histogram(mb)

ρgrid = LinRange(0.,10.,100)
plot(ρgrid, [modratio(x,F.τgrid,cprobs) for x in ρgrid])