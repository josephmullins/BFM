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


# now let's compare interpolating in a loop vs a big thing

V = zeros(8,10_000)
V2 = zeros(10,10,10,8,10)
grid = LinRange(0,3,8)
function looped_interp(V,grid)
    Interpolations.scale(interpolate(V,(NoInterp(),NoInterp(),NoInterp(),BSpline(Cubic(Line(OnGrid()))),NoInterp())),1:10,1:10,1:10,grid,1:10)
    # for i in axes(V,2)
    #     @views Interpolations.scale(interpolate(V[:,i],BSpline(Cubic(Line(OnGrid())))),grid)
    # end
end

# itp = interpolate(V[:,1],BSpline(Cubic(Line(OnGrid()))))
# etp = scale(itp,LinRange(0,3,8))

function big_interp(V,grid)
    N = size(V,2)
    Interpolations.scale(interpolate(V,(BSpline(Cubic(Line(OnGrid()))),NoInterp())),grid,1:N)
end

looped_interp(V2,grid)
big_interp(V,grid)

@time looped_interp(V2,grid);
@time big_interp(V,grid);