
# ----- Work choices ----- #
function work_probs(v1,v2,σ_L) #<- value of part-time and full-time work relative to not-working
    denom = 1 + exp(v1/σ_L) + exp(v2/σ_L)
    p1 = exp(v1/σ_L) / denom
    p2 = exp(v2/σ_L) / denom
    return p1, p2
end

# ----- Stage 5 ----- # 
function solve5!(mod)
    (;F) = mod
    (;N_t) = F
    for t in reverse(1:N_t) # current time index: N_t - s
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
        U0,U1 = IUW5(θ,F,eW,AW,κ)
        v0 = U0 + β*etp(eW,κ)
        v1 = U1 + β*etp(eW,κ+1)
        VW5[eW,κi,t] = σ_L*(log(exp((v0-v1)/σ_L) + 1) + v1/σ_L) #+ γ)
        vd5[1,eW,κi,t] = (v0 - v1)/σ_L # difference in choice utility of not working and working
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
    @unpack θ, F = mod
    @unpack σ_L,Λ_ϵ,Π_ϵ = θ
    @unpack N_t,N_d,N_ϵ,N_Agrid,N_κ,A_W,A_d,κ_W_grid,A_grid,A_bar,γ,β = F
    @views itp = interpolate(mod.VW4[:,:,t+1,:,:,:,:], (NoInterp(),NoInterp(),NoInterp(), NoInterp(), BSpline(Cubic(Line(OnGrid()))), BSpline(Cubic(Line(OnGrid())))))
    itp = Interpolations.scale(itp, 1:2,1:2,1:3, 1:5, κ_W_grid[t+1],A_grid)
    etp = extrapolate(itp,Interpolations.Flat())
    @views itp = interpolate(mod.VW5[:,t+1,:], (NoInterp(),BSpline(Cubic(Line(OnGrid())))))
    itp = Interpolations.scale(itp,1:2,κ_W_grid[t+1])
    etp5 = extrapolate(itp,Interpolations.Flat())
    AW = A_W[t]
    for eW=1:2,eH=1:2,di=1:N_d,ϵi=1:N_ϵ,κi=1:N_κ,ai=1:N_Agrid
        AH = AW + A_d[di]
        ϵH = Λ_ϵ[ϵi]
        AK = Int(A_grid[ai])
        κ = κ_W_grid[t][κi]
        U0,U1 = IUW4(θ,F,eW,eH,AW,κ,AH,ϵH,AK)
        if AK+1<A_bar
            v0 = U0
            v1 = U1
            for ϵii=1:N_ϵ
                v0 += β*Π_ϵ[ϵi,ϵii]*etp(eW,eH,di,ϵii,κ,AK+1)
                v1 += β*Π_ϵ[ϵi,ϵii]*etp(eW,eH,di,ϵii,κ+1,AK+1)
            end
        else
            v0 = U0 + β*etp5(eW,κ)
            v1 = U1 + β*etp5(eW,κ+1)
        end
        mod.vd4[eW,eH,t,di,ϵi,κi,ai] = (v0 - v1)/σ_L # diffrence in choice utility of not working and working
        mod.VW4[eW,eH,t,di,ϵi,κi,ai] = σ_L*(log(exp((v0-v1)/σ_L) + 1) + v1/σ_L) #+ γ)
    end
end

function iterateH4!(mod,t)
    @unpack θ, F = mod
    @unpack σ_L,Λ_ϵ,Π_ϵ = θ
    @unpack N_t,N_d,N_ϵ,N_Agrid,N_κ,A_W,A_d,κ_W_grid,A_grid,A_bar,γ,β = F
    @views itp = interpolate(mod.VH4[:,:,t+1,:,:,:,:], (NoInterp(),NoInterp(),NoInterp(), NoInterp(), BSpline(Cubic(Line(OnGrid()))), BSpline(Cubic(Line(OnGrid())))))
    itp = Interpolations.scale(itp, 1:2,1:2,1:3, 1:5, κ_W_grid[t+1], A_grid)
    etp = extrapolate(itp,Interpolations.Flat())
    AW = A_W[t]
    for eW=1:2,eH=1:2,di=1:N_d,ϵi=1:N_ϵ,κi=1:N_κ,ai=1:N_Agrid
        AH = AW + A_d[di]
        ϵH = Λ_ϵ[ϵi]
        AK = Int(A_grid[ai])
        κ = κ_W_grid[t][κi]
        U0,U1 = IUH4(θ,F,eH,AH,ϵH,AK)
        if AK+1<A_bar
            v0 = U0
            v1 = U1
            for ϵii=1:N_ϵ
                v0 += β*Π_ϵ[ϵi,ϵii]*etp(eW,eH,di,ϵii,κ,AK+1)
                v1 += β*Π_ϵ[ϵi,ϵii]*etp(eW,eH,di,ϵii,κ+1,AK+1)
            end
        else
            cv = dot(mod.VH5[eH,t+1,di,:],Π_ϵ[ϵi,:])
            v0 = U0 + cv
            v1 = U1 + cv
        end
        pW = 1/(1+exp(mod.vd4[eW,eH,t,di,ϵi,κi,ai]))
        mod.VH4[eW,eH,t,di,ϵi,κi,ai] = pW*v1 + (1-pW)*v0
    end
end

# ------ Stage 3 ----- #

# ------ Stage 2 ----- #

# ------ Stage 1 ----- #
