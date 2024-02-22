
function U1(θ,eW,eH,AW,AH,ϵH,κ,AK)
    (; α_C, α_l, γ_YW, γ_YH,αΓ_τWa) = θ

    logY_W = γ_YW[eW,1] + γ_YW[eW,2]*κ + γ_YW[eW,3]*κ^2 # N_t x N_κ
    logY_H = γ_YH[eH,1] + γ_YH[eH,2]*AH + γ_YH[eH,3]*AH^2 + ϵH + log(40)
    uc1 = α_C*log(exp(logY_H) + exp(logY_W)*20)
    uc2 = α_C*log(exp(logY_H) + exp(logY_W)*40)

    u0 = α_l * αΓ_τWa[AK+1]*log(112)
    u1 = α_l * αΓ_τWa[AK+1]*log(92)
    u2 = α_l * αΓ_τWa[AK+1]*log(72)

    UW0 = α_C*logY_H + α_l*log(112) #+ u0
    UW1 = uc1 + α_l*log(92) #+ u1
    UW2 = uc2 + α_l*log(72) #+ u2
    return UW0,UW1,UW2
end

function U2(θ,F,edW,edH,AW,κ,AH,ϵH,AK)
    (; α_νd, α_C, α_l, γ_YW, γ_YH, αΓ_τWa, αΓ_τHa) = θ
    (; π_H) = F

    # Income
    logY_W = γ_YW[edW,1] + γ_YW[edW,2]*κ + γ_YW[edW,3]*κ^2 # N_t x N_κ
    logY_H = γ_YH[edH,1] + γ_YH[edH,2]*AH + γ_YH[edH,3]*AH^2 + ϵH + log(40)

    U0 =  α_C*(log(π_H) + logY_H) + α_l * (1 + αΓ_τWa[AK+1])*log(112) + αΓ_τHa[AK+1]*log(72)
    U1 = α_C*log(π_H*exp(logY_H) + exp(logY_W)*20) + α_l * (1 + αΓ_τWa[AK+1])*log(92)
    U2 = α_C*log(π_H*exp(logY_H) + exp(logY_W)*40) + α_l * (1 + αΓ_τWa[AK+1])*log(72)
    return (U0,U1,U2)
end

ϵH = θ.Λ_ϵ[1]
κ = 4.
u0,u1,u2 = U1(θ,1,1,20,20,ϵH,κ,1)
v0,v1,v2 = U2(θ,F,1,1,20,κ,20,ϵH,1)