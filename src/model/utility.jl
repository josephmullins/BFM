# -- indirect utility functions in each stage of the model --- #

## Stage 5

# eW: education group of wife
function IUW5(θ,F,eW,AW,κ)
     (;α_C, α_l, γ_YW) = θ
     (;b) = F
    # Log income of wife
    logY_W = γ_YW[eW,1] + γ_YW[eW,2]*κ + γ_YW[eW,3]*κ^2

    # the set of indirect utility (2: L=0,1,2; N_t: A_W=18:59; N_κ: κ_W=0:A_W)
    U0 = α_C*log(max(b,1)) + α_l*log(112)
    U1 = α_C*log(max(b,exp(logY_W)*20)) + α_l*log(92)
    U2 = α_C*log(max(b,exp(logY_W)*40)) + α_l*log(72)
    return (U0, U1, U2)
end
# eH: education group of husband
function IUH5(θ,eH,AH,ϵH)
    (;α_C, α_l, γ_YH) = θ

    logY_H = γ_YH[eH,1] + γ_YH[eH,2]*AH + γ_YH[eH,3]*AH^2 + ϵH
    U = α_C*(log(40) + logY_H) + α_l*log(72)
    return U
end

## Stage 4

# utility function
function IUW4(θ,F,edW,edH,AW,κ,AH,ϵH,AK)
    (; α_C, α_l, γ_YW, γ_YH, αΓ_τWa, αΓ_τHa) = θ
    (; π_H) = F

    # Income
    logY_W = γ_YW[edW,1] + γ_YW[edW,2]*κ + γ_YW[edW,3]*κ^2 # N_t x N_κ
    logY_H = γ_YH[edH,1] + γ_YH[edH,2]*AH + γ_YH[edH,3]*AH^2 + ϵH + log(40)

    U0 =  α_C*(log(π_H) + logY_H) + α_l * (1 + αΓ_τWa[AK+1])*log(112) + αΓ_τHa[AK+1]*log(72)
    U1 = α_C*log(π_H*exp(logY_H) + exp(logY_W)*20) + α_l * (1 + αΓ_τWa[AK+1])*log(92)
    U2 = α_C*log(π_H*exp(logY_H) + exp(logY_W)*40) + α_l * (1 + αΓ_τWa[AK+1])*log(72)
    return (U0,U1,U2)
end

function IUH4(θ,F,edH,AH,ϵH,AK)
    (; α_C, α_l, γ_YH, αΓ_τWa, αΓ_τHa) = θ
    (;π_H) = F

    # Income
    logY_H = γ_YH[edH,1] + γ_YH[edH,2]*AH + γ_YH[edH,3]*AH^2 + ϵH + log(40)

    U = α_C*(log(1-π_H) + logY_H) + α_l * (1 + αΓ_τHa[AK+1])*log(72)
    U0 = U + α_l * αΓ_τWa[AK+1]*log(112)
    U1 = U + α_l * αΓ_τWa[AK+1]*log(92)
    U2 = U + α_l * αΓ_τWa[AK+1]*log(72)
    return (U0, U1, U2)
end

## Stage 3
function IU3(θ,eW,eH,AW,AH,ϵH,κ,ω)
    (; α_C, α_l, γ_YW, γ_YH,α_ω) = θ

    logY_W = γ_YW[eW,1] + γ_YW[eW,2]*κ + γ_YW[eW,3]*κ^2 # N_t x N_κ
    logY_H = γ_YH[eH,1] + γ_YH[eH,2]*AH + γ_YH[eH,3]*AH^2 + ϵH + log(40)
    vω = α_ω[1] + α_ω[2]*ω
    uc1 = α_C*log(exp(logY_H) + exp(logY_W)*20)
    uc2 = α_C*log(exp(logY_H) + exp(logY_W)*40)

    UW0 = α_C*logY_H + α_l*log(112) + vω
    UW1 = uc1 + α_l*log(72) + vω
    UW2 = uc2 + α_l*log(72) + vω

    UH0 = α_C*logY_H + α_l*log(72) + vω
    UH1 = uc1 + α_l*log(72) + vω
    UH2 = uc2 + α_l*log(72) + vω
    return (UW0,UW1,UW2,UH0,UH1,UH2)
end

## Stage 2

function IU2(θ,eW,eH,AW,AH,ϵH,κ,ω,AK)
    (; α_l, αΓ_τHa,αΓ_τWa,α_νH0,α_νH1) = θ

    UW0,UW1,UW2,UH0,UH1,UH2 = IU3(θ,eW,eH,AW,AH,ϵH,κ,ω)
    u0 = α_l * αΓ_τHa[AK+1]*log(72) + α_l * αΓ_τWa[AK+1]*log(112)
    u1 = α_l * αΓ_τHa[AK+1]*log(72) + α_l * αΓ_τWa[AK+1]*log(92)
    u2 = α_l * αΓ_τHa[AK+1]*log(72) + α_l * αΓ_τWa[AK+1]*log(72)
    vω = α_νH0[1]+α_νH0[2]*AK + α_νH1[1]*ω +α_νH1[2]*AK*ω

    UW0 += u0 + vω
    UW1 += u1 + vω
    UW2 += u2 + vω

    UH0 += u0 + vω
    UH1 += u1 + vω
    UH2 += u2 + vω
    return UW0,UW1,UW2,UH0,UH1,UH2
end

# A note on stage 1: utilities are the same as stage 3
