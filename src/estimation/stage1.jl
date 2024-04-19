# function to return the estimates of ρ and σ for the Husband's Income Process
function get_ρ_σ(d)
    r = 0.
    denom = 0
    for i in 2:size(d,1)
        if d.FID[i]==d.FID[i-1]
            r += d.eps[i]*d.eps[i-1]
            denom += 1
        end
    end
    v = var(skipmissing(d.eps))
    ρ = (r / denom) / v
    σ = sqrt(v*(1-ρ^2))
    return ρ,σ
end

function stage1(d)
    # income parameters for Husbands
    dh1 = @subset(d,:edH.==1,:F_wage.>0)
    mod = reg(dh1,@formula(log(F_wage) ~ F_AGE + F_AGE^2))
    γ_YCH = coef(mod)
    dh1[!,:eps] = residuals(mod,dh1)

    dh0 = @subset(d,:edH.==0,:F_wage.>0)
    mod = reg(dh0,@formula(log(F_wage) ~ F_AGE + F_AGE^2))
    γ_YNH = coef(mod)
    dh0[!,:eps] = residuals(mod,dh0)
    
    ρ_ϵ,σ_η = get_ρ_σ([dh0;dh1])

    # ----  income parameters for Wives/Mothers
    γ_YCW = @chain d begin 
        @subset :edW.==1 :M_wage.>0
        reg(@formula(log(M_wage) ~ exp2 + exp2^2))
        coef(_)
    end

    γ_YNW = @chain d begin 
        @subset :edW.==0 :M_wage.>0
        reg(@formula(log(M_wage) ~ exp2 + exp2^2))
        coef(_)
    end

    return γ_YNH,γ_YCH,γ_YNW,γ_YCW,ρ_ϵ,σ_η
end

function stage1(θ,F,d)
    (;N_ϵ) = F
    γ_YNH,γ_YCH,γ_YNW,γ_YCW,ρ_ϵ,σ_η = stage1(d)
    θ.γ_YH[1,:] .= γ_YNH
    θ.γ_YH[2,:] .= γ_YCH
    θ.γ_YW[1,:] .= γ_YNW
    θ.γ_YW[2,:] .= γ_YCW
    Λ_ϵ, Π_ϵ = rouwenhorst(N_ϵ,σ_η,ρ_ϵ)
    θ = (;θ..., ρ_ϵ = ρ_ϵ, σ_η = σ_η, Π_ϵ = Π_ϵ, Λ_ϵ = Λ_ϵ)
    return θ
end

# using DataFrames, DataFramesMeta, CSV
# D = CSV.read("data/MarriagePanel.csv",DataFrame,missingstring="NA")
# using GLM

# @chain D begin 
#     @subset :edW.==0 :M_wage.>0
#     reg(@formula(log(M_wage) ~ exp2 + exp2^2),_)
# end

# comparison with fixed effects
# @chain D begin 
#     @subset :edW.==0 :M_wage.>0
#     reg(_,@formula(log(M_wage) ~ exp2 + exp2^2 + fe(MID)))
# end



# Husband's income parameters


# Wife with probit adjustment
# NOTE: this code below shows how to do the probit adjustment and shows that the parameters can sure differ

# d = @chain D begin 
#     @select :MID :F_earn :M_AGE :exp2 :edW :M_hrs :M_wage
#     dropmissing()
# end
# pmod = glm(@formula(M_hrs>0 ~ edW*M_AGE*F_earn*exp2),d,Binomial(),ProbitLink())
# d[!,:pZ] = predict(pmod)

# @chain d begin
#     @subset :M_wage.>0 :edW.==0
#     reg(@formula(log(M_wage) ~ exp2 + exp2^2 + pZ + pZ^2),_)
# end