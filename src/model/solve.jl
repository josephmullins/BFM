
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

# ---- Expected values from work decision
function work_decision(vw0,vw1,vw2,vh0,vh1,vh2,σ)
    p1, p2 = work_probs(((vw1-vw0) / σ , (vw2-vw0) / σ))
    evw = inclusive_value(σ,(vw0,vw1,vw2))
    evh = vh0 + p1 * (vh1 - vh0) + p2 * (vh2 - vh0)
    return evw, evh
end

# ---- Function to return the value and probability of a unilateral decision
#   - choice 1 is the default unless both agents can agree on 0
function unilateral_decision(vw0,vw1,vh0,vh1,σ)
    # inclusive values (when get to choose)
    evw = inclusive_value(σ,(vw0,vw1))
    evh = inclusive_value(σ,(vh0,vh1))
    
    # probabilities
    pw = 1 / (1 + exp(vw0-vw1))
    ph = 1 / (1 + exp(vh0-vh1))
    pr = pw + (1-pw)*ph

    # expected values
    vw = ph * vw1 + (1-ph) * evw
    vh = pw * vh1 + (1-pw) * evh
    return vw,vh,pr
end

# ---- Function to return the value and probability of a mutual decision
#   - choice 0 is the default unless agents can agree on 1
function bilateral_decision(vw0,vw1,vh0,vh1,σ)
    # inclusive values (when get to choose)
    evw = inclusive_value(σ,(vw0,vw1))
    evh = inclusive_value(σ,(vh0,vh1))
    
    # probabilities
    pw = 1 / (1 + exp(vw0-vw1))
    ph = 1 / (1 + exp(vh0-vh1))
    pr = pw * ph 
    
    # expected values
    vw = ph * evw + (1-ph) * vw0
    vh = pw * evh + (1-pw) * vh0
    return vw,vh,pr
end

# ---- function to solve each stage
function solve_all!(mod)
    solve5!(mod)
    solve4!(mod)
    solve3!(mod)
    solve2!(mod)
    solve1!(mod)
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
    (;VW5,pL5) = values
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
        @views pL5[:,eW,κi,t] .= work_probs(((v1 - v0)/σ_L,(v2 - v0)/σ_L))
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
    (;pL4,VW4,VW5) = values
    (;N_d,N_ϵ,A_W,A_d,κ_W_grid,A_grid,A_bar,β) = F
    
    # interpolate VW4 next period
    etp = interpolate4(VW4,F,t)

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

        @views pL4[:,eW,eH,di,ϵi,κi,ai,t] .= work_probs(((v1 - v0)/σ_L,(v2 - v0)/σ_L)) # diffrence in choice utility of not working and working
        VW4[eW,eH,di,ϵi,κi,ai,t] = inclusive_value(σ_L,(v0,v1,v2))
    end
end

function iterateH4!(mod,t)
    (;θ, F, values) = mod
    (;Λ_ϵ,Π_ϵ) = θ
    (;VH4,VH5,pL4) = values
    (;A_W,A_d,κ_W_grid,A_grid,A_bar,β) = F

    etp = interpolate4(VH4,F,t)

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
        @views p1, p2 = pL4[:,eW,eH,di,ϵi,κi,ai,t]
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
    (;VW5,VH5,VW3,VH3,pL3,pD3) = values

    etpW = interpolate3(VW3,F,t)
    etpH = interpolate3(VH3,F,t)

    AW = A_W[t]
    idx = CartesianIndices((N_ϵ,N_ω,N_d,2,2,2))
    lin_idx = LinearIndices(idx)
    @Threads.threads for i in axes(VW3,2)
        ϵi,ωi,di,eH,eW,l = Tuple(idx[i])
        for κi in axes(VW3,1)
            κ = κ_W_grid[t][κi]
            ω = ω_grid[ωi]
            ϵH = Λ_ϵ[ϵi]
            AH = AW + A_d[di]
            vW0,vW1,vW2,vH0,vH1,vH2 = IU3(θ,eW,eH,AW,AH,ϵH,κ,ω)
            # --- work decision
            for ωii in 1:N_ω, ϵii in 1:N_ϵ
                i_next = lin_idx[ϵii,ωii,di,eH,eW,l]
                vW0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW(κ,i_next)
                vW1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW(κ+0.5,i_next)
                vW2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW(κ+1,i_next)
                vH0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH(κ,i_next)
                vH1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH(κ+0.5,i_next)
                vH2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH(κ+1,i_next)
            end
            # ---- Work Decision
            @views pL3[:,κi,i,t]  .= work_probs(((vW1 - vW0)/σ_L,(vW2 - vW0)/σ_L))
            VW,VH = work_decision(vW0,vW1,vW2,vH0,vH1,vH2,σ_L)

            # ---- Divorce Decision
            # continuation values for getting divorced
            vw5 = VW5[eW,κi,t]
            vh5 = VH5[eH,di,ϵi,t]

            # Finally, Emax values of arriving in this state
            if l==1 # mutual consent
                VW3[κi,i,t],VH3[κi,i,t],pD3[κi,i,t] = bilateral_decision(VW,vw5,VH,vh5,σ_ω)
            else            # unilateral
                VW3[κi,i,t],VH3[κi,i,t],pD3[κi,i,t] = unilateral_decision(VW,vw5,VH,vh5,σ_ω)
            end
            # store the value differences for simulation
        end
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
    (;VW2,VH2,VW3,VH3,VW4,VH4,pL2,pD2) = values

    etpW2 = interpolate2(VW2,F,t)
    etpH2 = interpolate2(VH2,F,t)

    etpW3 = interpolate3(VW3,F,t)
    etpH3 = interpolate3(VH3,F,t)

    AW = A_W[t]
    idx = CartesianIndices((N_ϵ,N_ω,N_d,2,2,2))
    lin_idx = LinearIndices(idx)
    @Threads.threads for i in axes(VW2,3)
        ϵi,ωi,di,eH,eW,l = Tuple(idx[i])
        for ai in axes(VW2,2)
            for κi in axes(VW2,1)
                AK = Int(A_grid[ai])
                κ = κ_W_grid[t][κi]
                ω = ω_grid[ωi]
                ϵH = Λ_ϵ[ϵi]
                AH = AW + A_d[di]
                vW0,vW1,vW2,vH0,vH1,vH2 = IU2(θ,eW,eH,AW,AH,ϵH,κ,ω,AK)
                if AK+1<A_bar
                    for ωii in 1:N_ω, ϵii in 1:N_ϵ
                        i_next = lin_idx[ϵii,ωii,di,eH,eW,l]
                        vW0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW2(κ,AK+1,i_next)
                        vW1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW2(κ+0.5,AK+1,i_next)
                        vW2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW2(κ+1,AK+1,i_next)
                        vH0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH2(κ,AK+1,i_next)
                        vH1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH2(κ+0.5,AK+1,i_next)
                        vH2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH2(κ+1,AK+1,i_next)
                    end
                else
                    for ωii in 1:N_ω, ϵii in 1:N_ϵ
                        i_next = lin_idx[ϵii,ωii,di,eH,eW,l]
                        vW0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW3(κ,i_next)
                        vW1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW3(κ+0.5,i_next)
                        vW2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW3(κ+1,i_next)
                        vH0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH3(κ,i_next)
                        vH1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH3(κ+0.5,i_next)
                        vH2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH3(κ+1,i_next)
                    end
                end
                #  ---- Work Decision
                pL2[:,κi,ai,i,t]  .= work_probs(((vW1 - vW0)/σ_L,(vW2 - vW0)/σ_L))
                VW,VH = work_decision(vW0,vW1,vW2,vH0,vH1,vH2,σ_L)

                # ---- Divorce Decision
                # continuation values for getting divorced
                νi = Int(AK+1)  # Index for expected value from the custody decision
                vw4 = VW4[eW,eH,di,ϵi,κi,ai,t] - Cτ[νi]
                vh4 = VH4[eW,eH,di,ϵi,κi,ai,t] - Cτ[νi]

                if l==1 # mutual consent
                    VW,VH,pD2[κi,ai,i,t] = bilateral_decision(VW,vw4,VH,vh4,σ_ω)
                else            # unilateral
                    VW,VH,pD2[κi,ai,i,t] = unilateral_decision(VW,vw4,VH,vh4,σ_ω)
                end
                VW2[κi,ai,i,t] = VW
                VH2[κi,ai,i,t] = VH
            end
        end
    end
end

# ------ Stage 1 ----- #
function iterate1!(mod,t)
    (;θ,F,values) = mod
    (;σ_L, σ_ω, Π_ϵ, Π_ω, Λ_ϵ, σ_F,α_F) = θ
    (;N_ϵ, N_ω, ω_grid, κ_W_grid, β,A_W,A_d,N_d,T_f) = F
    (;VW1,VH1,VW2,VH2,VW3,VH3,pL1,pD1,pF,VW5,VH5) = values

    #VW & VH : T_f-1 x N_d x N_ϵ x N_κ x N_ω
    etpW = interpolate3(VW1,F,t)
    etpH = interpolate3(VH1,F,t)

    etpW3 = interpolate3(VW3,F,t)
    etpH3 = interpolate3(VH3,F,t)

    AW = A_W[t]
    idx = CartesianIndices((N_ϵ,N_ω,N_d,2,2,2))
    lin_idx = LinearIndices(idx)
    @Threads.threads for i in axes(VW1,2)
        ϵi,ωi,di,eH,eW,l = Tuple(idx[i])
        for κi in axes(VW3,1)
            κ = κ_W_grid[t][κi]
            ω = ω_grid[ωi]
            ϵH = Λ_ϵ[ϵi]
            AH = AW + A_d[di]
            vW0,vW1,vW2,vH0,vH1,vH2 = IU3(θ,eW,eH,AW,AH,ϵH,κ,ω)

            if t+1<T_f
                for ωii in 1:N_ω,ϵii in 1:N_ϵ
                    i_next = lin_idx[ϵii,ωii,di,eH,eW,l]
                    vW0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW(κ,i_next)
                    vW1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW(κ+0.5,i_next)
                    vW2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW(κ+1,i_next)
                    vH0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH(κ,i_next)
                    vH1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH(κ+0.5,i_next)
                    vH2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH(κ+1,i_next)
                end
            else
                for ωii=1:N_ω,ϵii=1:N_ϵ
                    i_next = lin_idx[ϵii,ωii,di,eH,eW,l]
                    vW0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW3(κ,i_next)
                    vW1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW3(κ+0.5,i_next)
                    vW2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpW3(κ+1,i_next)
                    vH0 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH3(κ,i_next)
                    vH1 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH3(κ+0.5,i_next)
                    vH2 += β*Π_ϵ[ϵi,ϵii]*Π_ω[ωi,ωii]*etpH3(κ+1,i_next)
                end
            end

            # ----- Work Decision
            @views pL1[:,κi,i,t]  .= work_probs(((vW1 - vW0)/σ_L,(vW2 - vW0)/σ_L))

            VW,VH = work_decision(vW0,vW1,vW2,vH0,vH1,vH2,σ_L)

            # ----- Divorce Decision
            # continuation values for getting divorced
            vw5 = VW5[eW,κi,t]
            vh5 = VH5[eH,di,ϵi,t]

            if l==1 #<- mutual consent case
                VW,VH,pD1[κi,i,t] = bilateral_decision(VW,vw5,VH,vh5,σ_ω)
            else #<- unilateral case
                VW,VH,pD1[κi,i,t] = unilateral_decision(VW,vw5,VH,vh5,σ_ω)
            end

            # ------ Fertility Decision
            # continuation values for having a kid
            vwF = VW2[κi,1,i,t] + α_F
            vhF = VH2[κi,1,i,t] + α_F

            VW,VH,pF[κi,i,t] = bilateral_decision(VW,vwF,VH,vhF,σ_F)

            # save the value of arriving in this time period at stage 1
            VW1[κi,i,t] = VW
            VH1[κi,i,t] = VH
        end
    end
end

function solve1!(mod)
    (;F) = mod
    (;T_f) = F
    for t in reverse(1:T_f-1) # current period T_f-1-s
        iterate1!(mod,t)
    end
end


