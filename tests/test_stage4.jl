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

L = convert(Vector{Int64},P.L)
moms0 = data_moms(M,P)

dat = repeat(M,10)
legal = repeat(L,10)


θ = (;θ...,
    σ_L = 2., α_l = 1.,
    α_νd = 0.,σ_F = 10.,
    α_ω = [0., 5.],σ_ω = 3., # σ_ω = 2.
    π_ω = 0.7, Π_ω = transmat_ω(0.7,F.N_ω),
    α_νH0 = [3.,-0.1],α_νH1 = [7.,0.2])


x0 = get_x(θ)
mod = (;θ,values=V,F)
data_gen(mod,dat,legal);
@time data_gen(mod,dat,legal)

break


wght = [[1.,0.,1.,0.];ones(10)]

res1 = optimize(x->ssq(update(x,θ,F),wght,V,F,moms0,dat,legal),x0,Optim.Options(iterations=200,show_trace=true))
θ = update(res1.minimizer,θ,F)
m1 = get_moments(θ,V,F,dat,legal)
display([moms0 m1])


break

# blocks = [
#     [5,6,7,9,11,12],
#     [3,4,5,7],
#     [3,5,7,8,9,10,11],
#     [1,2]
# ]
# θ = estimate_blocks(θ,blocks,wght,V,F,moms0,M,L)
# θ = estimate_blocks(θ,blocks,wght,V,F,moms0,M,L)


x0 = get_x(θ)

wght = ones(14); wght[11] = 50.
res1 = optimize(x->ssq(update(x,θ,F),wght,V,F,moms0,M,L),x0,Optim.Options(iterations=200,show_trace=true))
res1 = optimize(x->ssq(update(x,θ,F),wght,V,F,moms0,M,L),res1.minimizer,Optim.Options(iterations=1000,show_trace=true))
θ = update(res1.minimizer,θ,F)
m1 = get_moments(θ,V,F,M,L)
#display([moms0 m1])

res2 = optimize(x->ssq(update(x,θ,F),V,F,moms0,M,L),x0,Optim.Options(iterations=200,show_trace=true))
res2 = optimize(x->ssq(update(x,θ,F),V,F,moms0,M,L),res2.minimizer,Optim.Options(iterations=1000,show_trace=true))
θ = update(res2.minimizer,θ,F)
m2 = get_moments(θ,V,F,M,L)

mstring = ["FT - married","PT - married","FT - divorced","PT - divorced"]
mstring = [mstring;["divorce < $x years" for x in (5,10,15)]]
mstring = [mstring;["birth < $x years" for x in (2,4,6)]]
mstring = [mstring;["time to divorce - time to birth < $x years" for x in (0,5,10,15)]]
display([mstring moms0 m1 m2])

