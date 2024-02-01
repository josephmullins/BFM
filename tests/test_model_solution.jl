include("../src/model/model.jl")

F = FixedParams()
θ = Params(F)
V = values(F)

mod = (;F,θ,values=V);


@show "Stage five testing"
iterateW5!(mod,10)
iterateH5!(mod,10)

@time iterateW5!(mod,10)
@time iterateH5!(mod,10)

solve5!(mod)
@show solve5!(mod)

@show "Stage four testing"
iterateW4!(mod,10)
iterateH4!(mod,10)

@time iterateW4!(mod,10)
@time iterateH4!(mod,10)

solve4!(mod)
@show solve4!(mod)


@show "Stage three testing"
# these don't work yet
# iterateW4!(mod,10)
# iterateH4!(mod,10)

# @time iterateW4!(mod,10)
# @time iterateH4!(mod,10)

# solve4!(mod)
# @show solve4!(mod)


