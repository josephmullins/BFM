
function interpolateW5(V,F,t)
    (;κ_W_grid) = F
    @views etp = extrapolate(
        Interpolations.scale(
            interpolate(V[:,:,t+1], (NoInterp(),BSpline(Cubic(Line(OnGrid()))))),
        1:2,κ_W_grid[t+1]),
    Interpolations.Flat())
    return etp
end

function interpolateW4(V,F,t)
    (;N_d,N_ϵ,κ_W_grid,A_grid) = F

    @views etp = extrapolate(
        Interpolations.scale(
            interpolate(
                V[:,:,:,:,:,:,t+1], (NoInterp(),NoInterp(),NoInterp(), NoInterp(), BSpline(Cubic(Line(OnGrid()))), BSpline(Cubic(Line(OnGrid()))))), 
                1:2,1:2,1:N_d, 1:N_ϵ, κ_W_grid[t+1],A_grid),
    Interpolations.Flat())
    return etp
end

function interpolate3(V,F,t)
    (;N_d,N_ϵ,κ_W_grid,N_ω) = F
    @views etp = extrapolate(
        Interpolations.scale(
            interpolate(
                V[:,:,:,:,:,:,:,t+1], (NoInterp(),NoInterp(),NoInterp(),NoInterp(),NoInterp(),BSpline(Cubic(Line(OnGrid()))),NoInterp())),
            1:2,1:2,1:2, 1:N_d, 1:N_ϵ, κ_W_grid[t+1], 1:N_ω),
        Interpolations.Flat())                           # N_d x N_ϵ x N_κ x N_ω
    return etp
end


function interpolate2(V,F,t)
    (;N_d,N_ϵ,κ_W_grid,N_ω,A_grid) = F
    @views etp = extrapolate(
        Interpolations.scale(
            interpolate(
                V[:,:,:,:,:,:,:,:,t+1], (NoInterp(),NoInterp(),NoInterp(),NoInterp(),NoInterp(),BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid()))),NoInterp())),
            1:2, 1:2, 1:2, 1:N_d, 1:N_ϵ, κ_W_grid[t+1], A_grid, 1:N_ω),
        Interpolations.Flat())                           # N_d x N_ϵ x N_κ x N_ω
    return etp
end
