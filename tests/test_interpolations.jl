# it doesn't seem like the order of interpolation dimensions matters too much

using Interpolations

N_grid = 30

function interp2(V)
    itp = interpolate(V, (BSpline(Cubic(Line(OnGrid()))),NoInterp()))
end

function interp3(V)
    itp = interpolate(V, (NoInterp(),BSpline(Cubic(Line(OnGrid()))),NoInterp(),NoInterp()))
end


V2 = rand(N_grid,8)
V3 = rand(2,N_grid,2,2)
interp3(V3)
interp2(V2)
@time interp3(V3);
@time interp2(V2);

