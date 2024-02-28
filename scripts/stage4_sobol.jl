include("../src/model/model.jl")
include("../src/estimation/estimation.jl")

M = CSV.read("data/MarriageFile.csv",DataFrame,missingstring="NA")
P = CSV.read("data/MarriagePanel.csv",DataFrame,missingstring="NA")
K = CSV.read("data/KidPanel.csv",DataFrame,missingstring="NA")
cprobs = CSV.read("data/CustodyMomentsSimple.csv",DataFrame).p[:]

F = FixedParams()
θ = Params(F)
V = values(F);


# this works, but we need to set up the sequence.
θ = stage1(θ,F,P)
θ = stage2(θ,K)
θ = stage3(θ,K,F.τgrid,cprobs)

mod = (;θ,values=V,F)

dat = prep_sim_data(M,P;R = 10)
moms0 = data_moms(M,P)

θ_lb = (σ_L = 0.5,α_l = 0.1,α_νd = -2.,σ_F = 0.5,α_ω = [-1.,0.1],σ_ω = 0.5,π_ω = 0.5,α_νH0 = [-5.,-1.],α_νH1 = [0.01,-1.])
θ_ub = (σ_L = 5,α_l = 5.,α_νd = 2.,σ_F = 20.,α_ω = [5.,10.],σ_ω = 20.,π_ω = 0.999,α_νH0 = [5.,1.],α_νH1 = [5.,1.])

x_lb = get_x(θ_lb)
x_ub = get_x(θ_ub)

using Sobol, DelimitedFiles
s = SobolSeq(12)

function sobol_search(x_lb,x_ub,mod,moms0,dat;N=10_000)
    (;θ,values,F) = mod
    Qsave = zeros(N)
    Xsave = zeros(12,N)
    for n in 1:N
        u = next!(s)
        x = x_lb .+ u .* (x_ub .- x_lb)
        Xsave[:,n] .= x
        Qsave[n] = ssq(update(x,θ,F),values,F,moms0,dat)
        #println(n," ",Qsave[n])
    end
    ii = sortperm(Qsave)[1:10]
    X = Xsave[:,ii]
    return X,Qsave[ii]
end

X,Q = sobol_search(x_lb,x_ub,mod,moms0,dat ; N=20_000)

writedlm("output/Xsave_round1",X)
writedlm("output/Qsave_round1",Q)

# x0 = Xsave[:,argmin(Qsave)]
# res1 = optimize(x->ssq(update(x,θ,F),V,F,moms0,dat,legal),x0,Optim.Options(iterations=200,show_trace=true))
# θ = update(res1.minimizer,θ,F)
# m1 = get_moments(θ,V,F,dat,legal)
# display([moms0 m1])

x_lb = [minimum(X[i,:]) for i in axes(X,1)]
x_ub = [maximum(X[i,:]) for i in axes(X,1)]

X,Q = sobol_search(x_lb,x_ub,mod,moms0,dat ; N=20_000)

writedlm("output/Xsave_round2",X)
writedlm("output/Qsave_round2",Q)
