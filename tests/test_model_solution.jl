include("../src/model/model.jl")

F = FixedParams()
θ = Params(F)
V = values(F);

mod = (;F,θ,values=V);


@show "Stage five testing"
iterateW5!(mod,10)
iterateH5!(mod,10)

@time iterateW5!(mod,10)
@time iterateH5!(mod,10)

solve5!(mod)
@time solve5!(mod)

@show "Stage four testing"
iterateW4!(mod,10)
iterateH4!(mod,10)

@time iterateW4!(mod,10)
@time iterateH4!(mod,10)

solve4!(mod)
@time solve4!(mod)


@show "Stage three testing"
# these don't work yet
iterate3!(mod,10)
@time iterate3!(mod,10)

solve3!(mod)
@time solve3!(mod)

@show "Stage two testing"
# these don't work yet
iterate2!(mod,10);
@time iterate2!(mod,10);

solve2!(mod)
@time solve2!(mod)

