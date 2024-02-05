
function interpolateW5(V,F,t)
    (;κ_W_grid) = F
    @views etp = extrapolate(
        Interpolations.scale(
            interpolate(V[:,:,t+1], (NoInterp(),BSpline(Cubic(Line(OnGrid()))))),
        1:2,κ_W_grid[t]),
    Interpolations.Flat())
    return etp
end

function interpolate4(V,F,t)
    (;N_d,N_ϵ,κ_W_grid,A_grid) = F

    @views etp = extrapolate(
        Interpolations.scale(
            interpolate(
                V[:,:,:,:,:,:,t+1], (NoInterp(),NoInterp(),NoInterp(), NoInterp(), BSpline(Cubic(Line(OnGrid()))), BSpline(Cubic(Line(OnGrid()))))), 
                1:2,1:2,1:N_d, 1:N_ϵ, κ_W_grid[t],A_grid),
    Interpolations.Flat())
    return etp
end


function interpolate3(V,F,t)
    (;N_d,N_ϵ,κ_W_grid,N_ω) = F
    NN = 2 * 2 * 2 * N_d * N_ϵ * N_ω
    @views etp = extrapolate(
        Interpolations.scale(
            interpolate(
                V[:,:,t+1], (BSpline(Cubic(Line(OnGrid()))),NoInterp())),
                κ_W_grid[t], 1:NN),
        Interpolations.Flat())                          
    return etp
end


function interpolate2(V,F,t)
    (;N_d,N_ϵ,κ_W_grid,N_ω,A_grid) = F
    NN = 2 * 2 * 2 * N_d * N_ϵ * N_ω
    @views etp = extrapolate(
        Interpolations.scale(
            interpolate(
                V[:,:,:,t+1], (BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid()))),NoInterp())),
            κ_W_grid[t], A_grid, 1:NN),
        Interpolations.Flat())    # N_d x N_ϵ x N_κ x N_ω
    return etp
end

# why not interpolate probabilities instead of value differences?
# - I should try and write it this way and see if it works.
function interpolate_probs(values,F)
    (;pL5,pL4,pL3,pL2,pL1,pD3,pD2,pD1,pF) = values
    (;A_grid,N_d,N_ω,N_ϵ,κ_W_grid,N_t,T_f) = F

    # interpolate work probabilities
    @views pL5_itp = [extrapolate(
            Interpolations.scale(
                interpolate(pL5[:,:,:,t],(NoInterp(),NoInterp(),BSpline(Cubic(Line(OnGrid()))))),
                1:2,1:2,κ_W_grid[t]),
                Interpolations.Flat()) for t in 1:N_t]
    @views pL4_itp = [extrapolate(
        Interpolations.scale(
            interpolate(
                pL4[:,:,:,:,:,:,:,t],
                (NoInterp(),NoInterp(),NoInterp(),NoInterp(), NoInterp(), BSpline(Cubic(Line(OnGrid()))), BSpline(Cubic(Line(OnGrid()))))), 
    1:2,1:2,1:2,1:N_d, 1:N_ϵ, κ_W_grid[t],A_grid
    ),
        Interpolations.Flat()) for t in 1:N_t]
    NN = 2 * 2 * 2 * N_d * N_ϵ * N_ω
    @views pL3_itp = [extrapolate(
        Interpolations.scale(
            interpolate(
                pL3[:,:,:,t], (NoInterp(),BSpline(Cubic(Line(OnGrid()))),NoInterp())),
                1:2, κ_W_grid[t], 1:NN),
        Interpolations.Flat()) for t in 1:N_t]
    @views pL2_itp = [extrapolate(
        Interpolations.scale(
            interpolate(
                pL2[:,:,:,:,t], (NoInterp(),BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid()))),NoInterp())),
                1:2, κ_W_grid[t], A_grid, 1:NN),
        Interpolations.Flat()) for t in 1:N_t]
    @views pL1_itp = [extrapolate(
        Interpolations.scale(
            interpolate(
                pL1[:,:,:,t], (NoInterp(),BSpline(Cubic(Line(OnGrid()))),NoInterp())),
                1:2, κ_W_grid[t], 1:NN),
        Interpolations.Flat()) for t in 1:T_f-1]

    # interpolate divorce probabilities
    @views pD1_itp = [extrapolate(
        Interpolations.scale(
            interpolate(
                pD1[:,:,t], (BSpline(Cubic(Line(OnGrid()))),NoInterp())),
                κ_W_grid[t], 1:NN),
        Interpolations.Flat()) for t in 1:T_f-1]
    @views pD2_itp = [extrapolate(
        Interpolations.scale(
            interpolate(
                pD2[:,:,:,t], (BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid()))),NoInterp())),
                κ_W_grid[t], A_grid, 1:NN),
        Interpolations.Flat()) for t in 1:N_t]
    @views pD3_itp = [extrapolate(
        Interpolations.scale(
            interpolate(
                pD3[:,:,t], (BSpline(Cubic(Line(OnGrid()))),NoInterp())),
                κ_W_grid[t], 1:NN),
        Interpolations.Flat()) for t in 1:N_t]
    # fertility choice probability
    @views pF_itp = [extrapolate(
        Interpolations.scale(
            interpolate(
                pF[:,:,t], (BSpline(Cubic(Line(OnGrid()))),NoInterp())),
                κ_W_grid[t], 1:NN),
        Interpolations.Flat()) for t in 1:T_f-1]
    return pL1_itp,pL2_itp,pL3_itp,pL4_itp,pL5_itp,pD1_itp,pD2_itp,pD3_itp,pF_itp
end