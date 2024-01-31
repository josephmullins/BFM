
# ----- Work choices ----- #
function work_probs(v) #<- value of part-time and full-time work relative to not-working
    denom = 1 + exp(v[1]) + exp(v[2])
    p1 = exp(v[1]) / denom
    p2 = exp(v[2]) / denom
    return p1, p2
end

# ----- Stage 5 ----- # 
function solve5!(mod)
    (;F) = mod
    (;N_t) = F
    for t in reverse(1:N_t)
        iterateW5!(mod,t)
        iterateH5!(mod,t)
    end
end
function iterateW5!(mod,t)
    (; θ, F, values) = mod
    (;VW5,vd5) = values
    (; σ_L) = θ
    (;A_W,κ_W_grid,β) = F
    @views itp = interpolate(mod.VW5[:,:,t+1], (NoInterp(),BSpline(Cubic(Line(OnGrid())))))
    itp = Interpolations.scale(itp,1:2,κ_W_grid[t+1])
    etp = extrapolate(itp,Interpolations.Flat())
    AW = A_W[t]
    for κi in axes(VW5,2), eW in axes(VW5,1)
        κ = κ_W_grid[t][κi]
        U0,U1,U2 = IUW5(θ,F,eW,AW,κ)
        v0 = U0 + β*etp(eW,κ)
        v1 = U1 + β*etp(eW,κ+0.5)
        v2 = U2 + β*etp(eW,κ+1)
        vmax = max(v0,v1,v2)
        
        VW5[eW,κi,t] = vmax + σ_L*(log(exp((v0-vmax)/σ_L) + exp((v1-vmax)/σ_L) + exp((v2-vmax)/σ_L))) 
        vd5[1,eW,κi,t] = (v1 - v0)/σ_L # difference in choice utility of not working and working
        vd5[2,eW,κi,t] = (v2 - v0)/σ_L # difference in choice utility of not working and working
    end
end

function iterateH5!(mod,t)
    (; θ, F, values) = mod
    (;VH5) = values
    (; Λ_ϵ,Π_ϵ) = θ
    (; A_W,A_d,β) = F
    AW = A_W[t]
    for ϵi in axes(VH5,3), di in axes(VH5,2), eH in axes(VH5,1)
        ϵH = Λ_ϵ[ϵi]
        AH = AW + A_d[di]
        U = IUH5(θ,eH,AH,ϵH)
        @views VH5[eH,t,di,ϵi] = U + β*dot(VH5[eH,t+1,di,:],Π_ϵ[ϵi,:])
    end
end

# ------ Stage 4 ----- #
function solve4!(mod)
    for s = 0:N_t-1 # current time index: N_t - s
        iterateW4!(mod,N_t-s)
        iterateH4!(mod,N_t-s)
    end
end
function iterateW4!(mod,t)
    (;θ, F, values) = mod
    (;σ_L,Λ_ϵ,Π_ϵ) = θ
    (;vd4,VW4,VW5) = values
    (;N_d,N_ϵ,A_W,A_d,κ_W_grid,A_grid,A_bar,β) = F
    
    # interpolate VW4 next period
    @views itp = interpolate(VW4[:,:,:,:,:,:,t+1], (NoInterp(),NoInterp(),NoInterp(), NoInterp(), BSpline(Cubic(Line(OnGrid()))), BSpline(Cubic(Line(OnGrid())))))
    itp = Interpolations.scale(itp, 1:2,1:2,1:N_d, 1:N_ϵ, κ_W_grid[t+1],A_grid)
    etp = extrapolate(itp,Interpolations.Flat())

    # interpolate VW5
    @views itp = interpolate(VW5[:,:,t+1], (NoInterp(),BSpline(Cubic(Line(OnGrid())))))
    itp = Interpolations.scale(itp,1:2,κ_W_grid[t+1])
    etp5 = extrapolate(itp,Interpolations.Flat())

    # iterate over states
    AW = A_W[t]
    for ai in axes(VW4,6), κi in axes(VW4,5), ϵi in axes(VW4,4), di in axes(VW4,3),eH in axes(VW4,2), eW in axes(VW4,1)
        AH = AW + A_d[di]
        ϵH = Λ_ϵ[ϵi]
        AK = Int(A_grid[ai])
        κ = κ_W_grid[t][κi]
        U0,U1,U2 = IUW4(θ,F,eW,eH,AW,κ,AH,ϵH,AK)
        if AK+1<A_bar
            v0 = U0
            v1 = U1
            v2 = U2
            for ϵii in axes(Π,2)
                v0 += β*Π_ϵ[ϵi,ϵii]*etp(eW,eH,di,ϵii,κ,AK+1)
                v1 += β*Π_ϵ[ϵi,ϵii]*etp(eW,eH,di,ϵii,κ+0.5,AK+1)
                v2 += β*Π_ϵ[ϵi,ϵii]*etp(eW,eH,di,ϵii,κ+1,AK+1)
            end
        else
            v0 = U0 + β*etp5(eW,κ)
            v1 = U1 + β*etp5(eW,κ+0.5)
            v2 = U2 + β*etp5(eW,κ+1)
        end

        vd4[1,eW,eH,di,ϵi,κi,ai,t] = (v1 - v0)/σ_L # diffrence in choice utility of not working and working
        vd4[2,eW,eH,di,ϵi,κi,ai,t] = (v2 - v0)/σ_L # diffrence in choice utility of not working and working
        vmax = max(v0,v1,v2)
        VW4[eW,eH,di,ϵi,κi,ai,t] = vmax + σ_L*(log(exp((v0-vmax)/σ_L) + exp((v1-vmax)/σ_L) + exp((v2-vmax)/σ_L))) 
    end
end

function iterateH4!(mod,t)
    (;θ, F, values) = mod
    (;Λ_ϵ,Π_ϵ) = θ
    (;VH4,VH5,vd4) = values
    (;N_d,N_ϵ,A_W,A_d,κ_W_grid,A_grid,A_bar,β) = F
    @views itp = interpolate(VH4[:,:,:,:,:,:,t+1], (NoInterp(),NoInterp(),NoInterp(), NoInterp(), BSpline(Cubic(Line(OnGrid()))), BSpline(Cubic(Line(OnGrid())))))
    itp = Interpolations.scale(itp, 1:2,1:2,1:N_d, 1:N_ϵ, κ_W_grid[t+1], A_grid)
    etp = extrapolate(itp,Interpolations.Flat())
    AW = A_W[t]
    for ai in axes(VW4,6), κi in axes(VW4,5), ϵi in axes(VW4,4), di in axes(VW4,3),eH in axes(VW4,2), eW in axes(VW4,1)
        AH = AW + A_d[di]
        ϵH = Λ_ϵ[ϵi]
        AK = Int(A_grid[ai])
        κ = κ_W_grid[t][κi]
        U0,U1,U2 = IUH4(θ,F,eH,AH,ϵH,AK)
        if AK+1<A_bar
            v0 = U0
            v1 = U1
            v2 = U2
            for ϵii in axes(Π_ϵ,2)
                v0 += β*Π_ϵ[ϵi,ϵii]*etp(eW,eH,di,ϵii,κ,AK+1)
                v1 += β*Π_ϵ[ϵi,ϵii]*etp(eW,eH,di,ϵii,κ+0.5,AK+1)
                v2 += β*Π_ϵ[ϵi,ϵii]*etp(eW,eH,di,ϵii,κ+1,AK+1)
            end
        else
            @views cv = dot(VH5[eH,di,:,t+1],Π_ϵ[ϵi,:])
            v0 = U0 + cv
            v1 = U1 + cv
            v2 = U2 + cv
        end
        @views p1, p2 = work_probs(vd4[:,eW,eH,di,ϵi,κi,ai,t])
        mod.VH4[eW,eH,di,ϵi,κi,ai,t] = (1-p1-p2)*v0 + p1*v1 + p2*v2
    end
end

# ------ Stage 3 ----- #

# ------ Stage 2 ----- #

# ------ Stage 1 ----- #
