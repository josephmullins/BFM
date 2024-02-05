
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
    NN = 2 * 2 * 2 * N_d * N_ϵ * N_ω
    @views etp = extrapolate(
        Interpolations.scale(
            interpolate(
                V[:,:,t+1], (BSpline(Cubic(Line(OnGrid()))),NoInterp())),
                κ_W_grid[t+1], 1:NN),
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
            κ_W_grid[t+1], A_grid, 1:NN),
        Interpolations.Flat())    # N_d x N_ϵ x N_κ x N_ω
    return etp
end

# why not interpolate probabilities instead of value differences?
# - I should try and write it this way and see if it works.
function interpolate_vd(values,F)
    (;vdL5,vdL4,vdL3,vdL2,vdL1,vdWD3,vdWD2,vdWD1,vdHD3,vdHD2,vdHD1,vdWF1,vdHF1) = values
    (;N_c,A_d,A_grid,ω_grid,κ_W_grid,N_t,T_f = F)

    @views vdL5 = [extrapolate(Interpolations.scale(interpolate(mod.vd5[:,t,:],(NoInterp(),BSpline(Cubic(Line(OnGrid()))))),1:2,κ_W_grid[t]),Interpolations.Flat()) for t=1:N_t]

    @views vdL4 = [extrapolate(Interpolations.scale(interpolate(mod.vd4[:,:,t,:,:,:,:],(NoInterp(),NoInterp(),BSpline(Linear()),NoInterp(),BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid()))))),1:2,1:2,A_d, 1:5, κ_W_grid[t],A_grid),Interpolations.Flat()) for t=1:N_t]

    @views vdL3 = [extrapolate(Interpolations.scale(interpolate(mod.vdL3[:,:,:,t,:,:,:,:],(NoInterp(),NoInterp(),NoInterp(),BSpline(Linear()),NoInterp(),BSpline(Cubic(Line(OnGrid()))),NoInterp())),1:2,1:2,1:2,A_d, 1:5, κ_W_grid[t],1:5),Interpolations.Flat()) for t=1:N_t]

    @views vdL2 = [extrapolate(Interpolations.scale(interpolate(mod.vdL2[:,:,:,:,t,:,:,:,:,:],(NoInterp(),NoInterp(),NoInterp(),NoInterp(),BSpline(Linear()),NoInterp(),BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid()))),NoInterp())),1:N_c,1:2,1:2,1:2,A_d, 1:5, κ_W_grid[t],A_grid,1:5),Interpolations.Flat()) for t=1:N_t]

    @views vdL1 = [extrapolate(Interpolations.scale(interpolate(mod.vdL1[:,:,:,:,t,:,:,:,:],(NoInterp(),NoInterp(),NoInterp(),NoInterp(),BSpline(Linear()),NoInterp(),BSpline(Cubic(Line(OnGrid()))),NoInterp())),1:N_c,1:2,1:2,1:2,A_d, 1:5, κ_W_grid[t],1:5),Interpolations.Flat()) for t=1:T_f-1]

    @views vdWD3 = [extrapolate(Interpolations.scale(interpolate(mod.vdWD3[:,:,:,t,:,:,:,:],(NoInterp(),NoInterp(),NoInterp(),BSpline(Linear()),NoInterp(),BSpline(Cubic(Line(OnGrid()))),NoInterp())),1:2,1:2,1:2,A_d, 1:5, κ_W_grid[t],1:5),Interpolations.Flat()) for t=1:N_t]

    @views vdHD3 = [extrapolate(Interpolations.scale(interpolate(mod.vdHD3[:,:,:,t,:,:,:,:],(NoInterp(),NoInterp(),NoInterp(),BSpline(Linear()),NoInterp(),BSpline(Cubic(Line(OnGrid()))),NoInterp())),1:2,1:2,1:2,A_d, 1:5, κ_W_grid[t],1:5),Interpolations.Flat()) for t=1:N_t]

    @views vdWD2 = [extrapolate(Interpolations.scale(interpolate(mod.vdWD2[:,:,:,:,t,:,:,:,:,:],(NoInterp(),NoInterp(),NoInterp(),NoInterp(),BSpline(Linear()),NoInterp(),BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid()))),NoInterp())),1:N_c,1:2,1:2,1:2,A_d, 1:5, κ_W_grid[t], A_grid,1:5),Interpolations.Flat()) for t=1:N_t]

    @views vdHD2 = [extrapolate(Interpolations.scale(interpolate(mod.vdWD2[:,:,:,:,t,:,:,:,:,:],(NoInterp(),NoInterp(),NoInterp(),NoInterp(),BSpline(Linear()),NoInterp(),BSpline(Cubic(Line(OnGrid()))),BSpline(Cubic(Line(OnGrid()))),NoInterp())),1:N_c,1:2,1:2,1:2,A_d, 1:5, κ_W_grid[t], A_grid,1:5),Interpolations.Flat()) for t=1:N_t]

    @views vdWD1 = [extrapolate(Interpolations.scale(interpolate(mod.vdWD1[:,:,:,:,t,:,:,:,:],(NoInterp(),NoInterp(),NoInterp(),NoInterp(),BSpline(Linear()),NoInterp(),BSpline(Cubic(Line(OnGrid()))),NoInterp())),1:N_c,1:2,1:2,1:2,A_d, 1:5, κ_W_grid[t],1:5),Interpolations.Flat()) for t=1:T_f-1]

    @views vdHD1 = [extrapolate(Interpolations.scale(interpolate(mod.vdHD1[:,:,:,:,t,:,:,:,:],(NoInterp(),NoInterp(),NoInterp(),NoInterp(),BSpline(Linear()),NoInterp(),BSpline(Cubic(Line(OnGrid()))),NoInterp())),1:N_c,1:2,1:2,1:2,A_d, 1:5, κ_W_grid[t],1:5),Interpolations.Flat()) for t=1:T_f-1]

    @views vdWF1 = [extrapolate(Interpolations.scale(interpolate(mod.vdWF1[:,:,:,:,t,:,:,:,:],(NoInterp(),NoInterp(),NoInterp(),NoInterp(),BSpline(Linear()),NoInterp(),BSpline(Cubic(Line(OnGrid()))),NoInterp())),1:N_c,1:2,1:2,1:2,A_d, 1:5, κ_W_grid[t],1:5),Interpolations.Flat()) for t=1:T_f-1]

    @views vdHF1 = [extrapolate(Interpolations.scale(interpolate(mod.vdHF1[:,:,:,:,t,:,:,:,:],(NoInterp(),NoInterp(),NoInterp(),NoInterp(),BSpline(Linear()),NoInterp(),BSpline(Cubic(Line(OnGrid()))),NoInterp())),1:N_c,1:2,1:2,1:2,A_d, 1:5, κ_W_grid[t],1:5),Interpolations.Flat()) for t=1:T_f-1]

    return vdL5,vdL4,vdL3,vdL2,vdL1,vdWD3,vdHD3,vdWD2,vdHD2,vdWD1,vdHD1,vdWF1,vdHF1
end