# add to data frame so we can use a different status variable
function kidmoms_data(d) #<-?
    m = @chain d begin
        @subset :AGE.<=17
        @subset :AGE.>=3
        #@transform :AGE = round.(:AGE ./ 4)
        groupby([:dgroup,:AGE])
        @combine :S = mean(skipmissing(:AP_raw))
    end
    return m
end

# TODO:
# finish and test code to calculate test score moments in data (above)
# write code to update function with age profiles
# write code to fetch the simulated moments
# write code to get the objective
# write code to minimize and see if it's enough
# in data patterns, look at test score moments relative to never divorced group


# TODO: fix this!!!
function kid_moments(S,θk,θ,F,kid_data)
    predict_k!(S,θk,θ,F,kid_data)
    return [mean(S[kid_data.AK.==a .&& kid_data.G.==g]) for g in 1:3 for a in 3:17]
end

function obj_stage5(S,θk,θ,F,kid_data,kmoms0)
    kmoms1 = kid_moments(S,θk,θ,F,kid_data)
    return sum((kmoms1 .- kmoms0).^2)
end

function updateθk(x,θk,θ)
    (;γ_ψ,δW,δH,Γa) = θk
    (;αΓ_τWa,αΓ_τHa,β) = θ
    γ_ψ0 = x[1]
    γ_ψ[:] .= x[2:5]
    δk = exp(x[6])/(1+exp(x[6]))
    get_Γa!(Γa,δk,β)
    δW[1] = exp(x[7])
    δH[1] = exp(x[8])
    for a in eachindex(δW)
        δW[a] = δW[1] * αΓ_τWa[a] * Γa[2] / (αΓ_τWa[1] * Γa[a+1])
        δH[a] = δH[1] * αΓ_τHa[a] * Γa[2] / (αΓ_τHa[1] * Γa[a+1])
    end
    return (;θk...,γ_ψ0,δk)
end

function get_Γa!(Γa,δk,β)
    Γ_next = 1/(1-β) 
    for a in reverse(eachindex(Γa))
        Γa[a] = 1 + β * δk * Γ_next
        Γ_next = Γa[a]
    end
end

# recall the recursive formula:
# Γ_{s,a} = δ_{s,a} * β Γ_{a+1}
# Γ_{a} = 1 + β * δk* Γ_{a+1}
# so:
# δ_{s,a+1} = δ_{s,a} * Γ_{s,a+1} * Γ_{a+1} / (Γ_{s,a} * Γ_{a+2})
# Γ_{a+1} / Γ_{a+2} = (1 + β * δ_{k} * Γ_{a+2}) / Γ_{a+2}