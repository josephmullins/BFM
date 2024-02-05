include("../src/model/model.jl")
using Random
F = FixedParams()
θ = Params(F)
V = values(F);

mod = (;F,θ,values=V);

solve_all!(mod)


@views interpolate(V.pL5[:,:,:,1],(NoInterp(),NoInterp(),BSpline(Cubic(Line(OnGrid())))))

itp = interpolate_probs(V,F);

@time itp = interpolate_probs(V,F);

# eW = dat.eW[n]
# eH = dat.eW[n]
# ad = dat.AD[n]
# AW0 = max(dat.AW0[n],18)
# maxT = min(dat.maxT[n],60-AW0)
# κ = dat.κ0[n]
# l = dat.L[n]
N = 2000
dat = (eW = rand(1:2,N),eH = rand(1:2,N),AD = rand(1:3,N),AW0 = rand(18:30,N),maxT = fill(15,N),κ0 = fill(0,N),L = rand(1:2,N))

# very nice
@time d = data_gen(mod,dat)
