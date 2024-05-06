# load the estimates
include("../src/model/model.jl")
include("../src/estimation/estimation.jl")
using DelimitedFiles

# load the data
M = CSV.read("data/MarriageFile.csv",DataFrame,missingstring="NA")
P = CSV.read("data/MarriagePanel.csv",DataFrame,missingstring="NA")
K = CSV.read("data/KidPanel.csv",DataFrame,missingstring="NA")
cprobs = CSV.read("data/CustodyMomentsSimple.csv",DataFrame).p[:]


m = testscore_averages(K)

function testcore_diffs(m)
    m0 = @chain m begin
        @subset :dgroup.==1
        @select :AGE :S
        @rename :S0 = :S
    end
    m1 = @chain m begin
        @subset :dgroup.>1
        innerjoin(m0,on=:AGE)
        @transform :S = :S .- :S0
    end
    return m1
end

using StatsPlots
dm = testcore_diffs(m)
@df dm plot(:AGE,:S,group=:dgroup)

function boot_diffs(M,P,K ; B = 50)
    db = []
    for b in 1:B
        Mb,Pb,Kb = draw_boot_sample(M, P, K)
        m = testscore_averages(Kb)
        dm = testcore_diffs(m)
        dm[!,:b] .= b
        db = [db; dm]
    end
    return vcat(db...)
end

db = boot_diffs(M,P,K)

@chain db begin
    groupby([:dgroup,:AGE])
    @combine lb = quantile(:S,0.05) ub = quantile(:S,0.95)
    innerjoin(dm,on=[:dgroup,:AGE])
    @df begin
        plot(:AGE,:S,group=:dgroup)
        plot!(:AGE,:lb,group=:dgroup)
        plot!(:AGE,:ub,group=:dgroup)
    end
end

@chain db begin
    groupby([:dgroup,:AGE])
    @combine se = std(:S)
    innerjoin(dm,on=[:dgroup,:AGE])
    @transform :lb = :S .- :se :ub = :S .+ :se
    @df begin
        plot(:AGE,:S,group=:dgroup)
        plot!(:AGE,:lb,group=:dgroup)
        plot!(:AGE,:ub,group=:dgroup)
    end
end


# just look at conventional measures of differences here. 
# they are there but marginal

@chain K begin
    @subset :AGE.>14
    @subset .!ismissing.(:AP_std)
    groupby(:dgroup)
    @combine :m = mean(:AP_std) :se = std(:AP_std) :N = sum(.!ismissing.(:AP_raw))
    @transform :sd = :se ./ sqrt.(:N)
end

# NEXT: review how we constructed the dataset. there seem to be way to few observations here.

dd = @chain db begin
    groupby([:dgroup,:AGE])
    @combine lb = quantile(:S,0.05) ub = quantile(:S,0.95)
    #@combine se = std(:S)
    innerjoin(dm,on=[:dgroup,:AGE])
    #@transform :lb = :S .- :se :ub = :S .+ :se
end

plot()
for g in 2:3
    @chain dd begin
        @subset :dgroup.==g
        @df begin 
            plot!(:AGE,:S)
            errorline!(:AGE,:lb,:ub)
        end
    end
end

d2 = @subset(dd,:dgroup.==2)
d3 = @subset(dd,:dgroup.==3)

plot(d2.AGE,d2.S)
errorline!(d2.AGE,[d2.lb d2.ub])
plot!(d3.AGE,d3.S)
errorline!(d3.AGE,[d3.lb d3.ub])