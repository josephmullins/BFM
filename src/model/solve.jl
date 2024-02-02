
# ----- Work choices ----- #
function work_probs(v) #<- value of part-time and full-time work relative to not-working
    denom = 1. + exp(v[1]) + exp(v[2])
    p1 = exp(v[1]) / denom
    p2 = exp(v[2]) / denom
    return p1, p2
end
function inclusive_value(σ,values)
    vmax = max(values...)
    return vmax + σ * log(sum((exp((v - vmax)/σ) for v in values)))
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
    AW = A_W[t]
    etp = interpolateW5(VW5,F,t)
    for κi in axes(VW5,2), eW in axes(VW5,1)
        κ = κ_W_grid[t][κi]
        U0,U1,U2 = IUW5(θ,F,eW,AW,κ)
        v0 = U0 + β*etp(eW,κ)
        v1 = U1 + β*etp(eW,κ+0.5)
        v2 = U2 + β*etp(eW,κ+1)
        
        VW5[eW,κi,t] = inclusive_value(σ_L,(v0,v1,v2)) 
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
        @views VH5[eH,di,ϵi,t] = U + β*dot(VH5[eH,di,:,t],Π_ϵ[ϵi,:])
    end
end

# ------ Stage 4 ----- #
function solve4!(mod)
    (;F) = mod
    (;N_t) = F
    for t in reverse(1:N_t) # current time index: N_t - s
        iterateW4!(mod,t)
        iterateH4!(mod,t)
    end
end
function iterateW4!(mod,t)
    (;θ, F, values) = mod
    (;σ_L,Λ_ϵ,Π_ϵ) = θ
    (;vd4,VW4,VW5) = values
    (;N_d,N_ϵ,A_W,A_d,κ_W_grid,A_grid,A_bar,β) = F
    
    # interpolate VW4 next period
    etp = interpolateW4(VW4,F,t)

    # interpolate VW5
    etp5 = interpolateW5(VW5,F,t)

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
            for ϵii in axes(Π_ϵ,2)
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
        VW4[eW,eH,di,ϵi,κi,ai,t] = inclusive_value(σ_L,(v0,v1,v2))
    end
end

function iterateH4!(mod,t)
    (;θ, F, values) = mod
    (;Λ_ϵ,Π_ϵ) = θ
    (;VH4,VH5,vd4) = values
    (;N_d,N_ϵ,A_W,A_d,κ_W_grid,A_grid,A_bar,β) = F
    @views itp = interpolate(VH4[:,:,:,:,:,:,t+1], (NoInterp(),NoInterp(),NoInterp(), NoInterp(), BSpline(Cubic(Line(OnGrid()))), BSpline(Cubic(Line(OnGrid())))))
    itp_scale = Interpolations.scale(itp, 1:2,1:2,1:N_d, 1:N_ϵ, κ_W_grid[t+1], A_grid)
    etp = extrapolate(itp_scale,Interpolations.Flat())
    AW = A_W[t]
    for ai in axes(VH4,6), κi in axes(VH4,5), ϵi in axes(VH4,4), di in axes(VH4,3),eH in axes(VH4,2), eW in axes(VH4,1)
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
        VH4[eW,eH,di,ϵi,κi,ai,t] = (1-p1-p2)*v0 + p1*v1 + p2*v2
    end
end

# ------ Stage 3 ----- #
function solve3!(mod)
    (;F) = mod
    (;N_t) = F
    for t in reverse(1:N_t) # current time index: N_t - s
        iterate3!(mod,t)
    end
end

function iterate3!(mod,t)
    (;θ,F,values) = mod
    (;α_ω, σ_L, σ_ω, Λ_ϵ, α_ω, Π_ϵ, Π_ω) = θ
    (; N_d, N_ϵ,  N_ω, κ_W_grid, A_W, A_d, ω_grid, β, N_κ) = F
    (;VW5,VH5,VW3,VH3,vdL3,vdWD3,vdHD3) = values

    etpW = interpolate3(VW3,F,t)
    etpH = interpolate3(VH3,F,t)

    AW = A_W[t]
    @Threads.threads for idx in CartesianIndices((2,2,2,N_d,N_ϵ,N_κ,N_ω))
    l,eW,eH,di,ϵi,κi,ωi = Tuple(idx)
        κ = κ_W_grid[t][κi]
        ω = ω_grid[ωi]
        ϵH = Λ_ϵ[ϵi]
        AH = AW + A_d[di]
        vW0,vW1,vW2,vH0,vH1,vH2 = IU3(θ,eW,eH,AW,AH,ϵH,κ,ω)
        # --- work decision
        for ωii in 1:N_ω, ϵii in 1:N_ϵ
            vW0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW(l,eW,eH,di,ϵi,κ,ωii)
            vW1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW(l,eW,eH,di,ϵi,κ+0.5,ωii)
            vW2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW(l,eW,eH,di,ϵi,κ+1,ωii)
            vH0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH(l,eW,eH,di,ϵi,κ,ωii)
            vH1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH(l,eW,eH,di,ϵi,κ+0.5,ωii)
            vH2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW(l,eW,eH,di,ϵi,κ+1,ωii)
        end
        vdL3[1,l,eW,eH,di,ϵi,κi,ωi,t]  = (vW1 - vW0)/σ_L
        vdL3[2,l,eW,eH,di,ϵi,κi,ωi,t]  = (vW2 - vW0)/σ_L
        @views p1,p2 = work_probs(vdL3[:,l,eW,eH,di,ϵi,κi,ωi,t])

        # continuation values for staying married
        VW_tilde = inclusive_value(σ_L,(vW0,vW1,vW2))
        VH_tilde = vH0 + p1*(vH1 - vH0) + p2*(vH2 - vH0)

        # continuation values for getting divorced
        vw5 = VW5[eW,κi,t]
        vh5 = VH5[eH,di,ϵi,t]

        # Emax values if either agent can choose their preferred option:
        EVW3 = inclusive_value(σ_ω,(vw5,VW_tilde))
        EVH3 = inclusive_value(σ_ω,(vh5,VH_tilde))

        # probability that each agent prefers divorce:

        vd_wd = (VW_tilde - vw5)  / σ_ω
        vd_hd = (VH_tilde - vh5) / σ_ω
        pwd = 1 / (1 + exp(vd_wd)) # wife's prob of divorce
        phd = 1 / (1 + exp(vd_hd)) # husband's prob of divorce

        # Finally, Emax values of arriving in this state
        if l==1 # mutual consent
            VW3[l,eW,eH,di,ϵi,κi,ωi,t] = (1 - phd)*VW_tilde + phd*EVW3
            VH3[l,eW,eH,di,ϵi,κi,ωi,t] = (1 - pwd)*VH_tilde + pwd*EVH3
        else            # unilateral
            VW3[l,eW,eH,di,ϵi,κi,ωi,t] = phd*vw5   + (1 - phd)*EVW3
            VH3[l,eW,eH,di,ϵi,κi,ωi,t] = pwd*vh5 + (1 - pwd)*EVH3
        end
        # store the value differences for simulation
        vdWD3[l,eW,eH,di,ϵi,κi,ωi,t] = vd_wd
        vdHD3[l,eW,eH,di,ϵi,κi,ωi,t] = vd_hd
    end
end


# ------ Stage 2 ----- #
function solve2!(mod)
    (;F) = mod
    (;N_t) = F
    for t = reverse(1:N_t) # current time index: N_t - s
        iterate2!(mod,t)
    end
end


function iterate2!(mod,t)
    (;θ,F,values) = mod
    (;σ_L, σ_ω, Π_ϵ, Π_ω, Λ_ϵ, Cτ) = θ
    (;N_ϵ, N_ω, A_bar, A_grid, ω_grid, κ_W_grid, β,A_W,A_d,N_d,N_a,N_κ) = F
    (;VW2,VH2,VW3,VH3,VW4,VH4,vdL2,vdWD2,vdHD2) = values

    etpW2 = interpolate2(VW2,F,t)
    etpH2 = interpolate2(VH2,F,t)

    etpW3 = interpolate3(VW3,F,t)
    etpH3 = interpolate3(VH3,F,t)

    AW = A_W[t]
    @Threads.threads for idx in CartesianIndices((2,2,2,N_d,N_ϵ,N_κ,N_a,N_ω))
        l,eW,eH,di,ϵi,κi,ai,ωi = Tuple(idx)
        #ωi,ai,κi,ϵi,di,eH,eW,l = Tuple(idx)
        AK = Int(A_grid[ai])
        κ = κ_W_grid[t][κi]
        ω = ω_grid[ωi]
        ϵH = Λ_ϵ[ϵi]
        AH = AW + A_d[di]
        vW0,vW1,vW2,vH0,vH1,vH2 = IU2(θ,eW,eH,AW,AH,ϵH,κ,ω,AK)
        if AK+1<A_bar
            for ωii=1:N_ω,ϵii=1:N_ϵ
                vW0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW2(l,eW,eH,di,ϵi,κ,AK+1,ωii)
                vW1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW2(l,eW,eH,di,ϵi,κ+0.5,AK+1,ωii)
                vW2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW2(l,eW,eH,di,ϵi,κ+1,AK+1,ωii)
                vH0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH2(l,eW,eH,di,ϵi,κ,AK+1,ωii)
                vH1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH2(l,eW,eH,di,ϵi,κ+0.5,AK+1,ωii)
                vH2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH2(l,eW,eH,di,ϵi,κ+1,AK+1,ωii)
            end
        else
            for ωii=1:N_ω,ϵii=1:N_ϵ
                vW0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW3(l,eW,eH,di,ϵi,κ,ωii)
                vW1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW3(l,eW,eH,di,ϵi,κ+0.5,ωii)
                vW2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW3(l,eW,eH,di,ϵi,κ+1,ωii)
                vH0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH3(l,eW,eH,di,ϵi,κ,ωii)
                vH1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH3(l,eW,eH,di,ϵi,κ+0.5,ωii)
                vH2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH3(l,eW,eH,di,ϵi,κ+1,ωii)
            end
        end
        vdL2[1,l,eW,eH,di,ϵi,κi,ai,ωi,t]  = (vW1 - vW0)/σ_L
        vdL2[2,l,eW,eH,di,ϵi,κi,ai,ωi,t]  = (vW2 - vW0)/σ_L
        @views p1,p2 = work_probs(vdL2[:,l,eW,eH,di,ϵi,κi,ai,ωi,t])

        # continuation values for staying married
        VW_tilde = inclusive_value(σ_L,(vW0,vW1,vW2))
        VH_tilde = vH0 + p1*(vH1 - vH0) + p2*(vH2 - vH0)

        # continuation values for getting divorced
        νi = Int(AK+1)  # Index for expected value from the custody decision
        vw4 = VW4[eW,eH,di,ϵi,κi,ai,t] - Cτ[νi]
        vh4 = VH4[eW,eH,di,ϵi,κi,ai,t] - Cτ[νi]

        # Emax values if either agent can choose:
        #EVW2 = inclusive_value(σ_ω,(1.,1.))
        EVW2 = inclusive_value(σ_ω,(VW_tilde,vw4))
        EVH2 = inclusive_value(σ_ω,(VH_tilde,vh4))

        # probability that each agent prefers divorce:

        vd_wd = (VW_tilde - vw4)  / σ_ω
        vd_hd = (VH_tilde - vh4) / σ_ω
        pwd = 1 / (1 + exp(vd_wd)) # wife's prob of divorce
        phd = 1 / (1 + exp(vd_hd)) # husband's prob of divorce

        if l==1 # mutual consent
            VW2[l,eW,eH,di,ϵi,κi,ai,ωi,t] = (1 - phd)*VW_tilde + phd*EVW2
            VH2[l,eW,eH,di,ϵi,κi,ai,ωi,t] = (1 - pwd)*VH_tilde + pwd*EVH2
        else            # unilateral
            VW2[l,eW,eH,di,ϵi,κi,ai,ωi,t] = phd*vw4 + (1 - phd)*EVW2
            VH2[l,eW,eH,di,ϵi,κi,ai,ωi,t] = pwd*vh4 + (1 - pwd)*EVH2
        end
        vdWD2[l,eW,eH,di,ϵi,κi,ai,ωi,t] = vd_wd
        vdHD2[l,eW,eH,di,ϵi,κi,ai,ωi,t] = vd_hd
    end
end

# ------ Stage 1 ----- #
