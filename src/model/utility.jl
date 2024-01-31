# -- indirect utility functions in each stage of the model --- #

## Stage 5

# eW: education group of wife
function IUW5(θ,F,eW,AW,κ)
     (;α_C, α_ℓ, γ_YW) = θ
     (;b) = F
    # Log income of wife
    logY_W = γ_YW[eW,1] + γ_YW[eW,2]*AW + γ_YW[eW,3]*κ

    # the set of indirect utility (2: L=0,1; N_t: A_W=18:59; N_κ: κ_W=0:A_W)
    U0 = α_C*log(max(b,1)) + α_ℓ*log(112)
    U1 = α_C*log(max(b,exp(logY_W))) + α_ℓ*log(72)
    return (U0, U1)
end
# eH: education group of husband
function IUH5(θ,eH,AH,ϵH)
    (;α_C, α_ℓ, γ_YH) = θ

    logY_H = γ_YH[eH,1] + γ_YH[eH,2]*AH + ϵH
    U = α_C*logY_H + α_ℓ*log(72)
    return U
end

## Stage 4 (add part-time work decision?)

# utility function
function IUW4(θ,F,edW,edH,AW,κ,AH,ϵH,AK)
    (; α_C, α_ℓ, γ_YW, γ_YH, αΓ_τWa, αΓ_τHa) = θ
    (; π_H = F)

    # Income
    logY_W = γ_YW[edW,1] + γ_YW[edW,2]*AW + γ_YW[edW,3]*κ # N_t x N_κ
    logY_H = γ_YH[edH,1] + γ_YH[edH,2]*AH + ϵH

    U0 =  α_C*(log(π_H) + logY_H) + (α_ℓ + αΓ_τWa[AK+1])*log(112) + αΓ_τHa[AK+1]*log(72)
    U1 = α_C*log(π_H*exp(logY_H) + exp(logY_W)) + (α_ℓ + αΓ_τWa[AK+1])*log(72)
    return (U0,U1)
end

function IUH4(θ,F,edH,AH,ϵH,AK)
    (; α_C, α_ℓ, γ_YW, γ_YH, αΓ_τWa, αΓ_τHa) = θ
    π_H = F

    # Income
    logY_H = γ_YH[edH,1] + γ_YH[edH,2]*AH + ϵH

    U = α_C*(log(1-π_H) + logY_H) + (α_ℓ + αΓ_τHa[AK+1])*log(72)
    U0 = U + αΓ_τWa[AK+1]*log(112)
    U1 = U + αΓ_τWa[AK+1]*log(72)
    return (U0,U1)
end

## Stage 3
function IU3(θ,eW,eH,AW,AH,ϵH,κ,ω)
    (; α_C, α_ℓ, γ_YW, γ_YH,α_ω) = θ

    logY_W = γ_YW[eW,1] + γ_YW[eW,2]*AW + γ_YW[eW,3]*κ # N_t x N_κ
    logY_H = γ_YH[eH,1] + γ_YH[eH,2]*AH + ϵH
    vω = α_ω[1] + α_ω[2]*ω
    uc1 = α_C*log(exp(logY_H) + exp(logY_W))
    UW0 = α_C*logY_H + α_ℓ*log(112) + vω
    UW1 = uc1 + α_ℓ*log(72) + vω
    UH0 = α_C*logY_H + α_ℓ*log(72) + vω
    UH1 = uc1 + α_ℓ*log(72) + vω
    return (UW0,UW1,UH0,UH1)
end

## Stage 2

function IU2(θ,eW,eH,AW,AH,ϵH,κ,ω,AK)
    (; αΓ_τHa,αΓ_τWa,α_νH0,α_νH1) = θ

    UW0,UW1,UH0,UH1 = IU3(θ,eW,eH,AW,AH,ϵH,κ,ω)
    u0 = αΓ_τHa[AK+1]*log(72) + αΓ_τWa[AK+1]*log(112)
    u1 = αΓ_τHa[AK+1]*log(72) + αΓ_τWa[AK+1]*log(72)
    vω = α_νH0[1]+α_νH0[2]*AK + α_νH1[1]+α_νH1[2]*AK*ω

    #uωW = α_kWH*(α_νH0[1]+α_νH0[2]*AK) + α_ω[1] + (α_kWH*(α_νH1[1]+α_νH1[2]*AK) + α_ω[2])*ω_grid[ωi]
    #uωH = (α_νH0[1]+α_νH0[2]*AK) + α_ω[1] + ((α_νH1[1]+α_νH1[2]*AK) + α_ω[2])*ω_grid[ωi]

    UW0 += u0 + vω
    UW1 += u1 + vω
    UH0 += u0 + vω
    UH1 += u1 + vω
    
    return UW0,UW1,UH0,UH1
end

# A note on stage 1: utilities are the same as stage 3
