# tau_f = ϕ_{a} / (1+ρ(1-τbar)) 
# EE[PHI_{f}|divorced] = ϕ_{a} * EE[1/(1+ρ(1-τbar))]
# THEN ALSO: there is only a probability τbar of the father being seen on any given day, so this gives:
# EE[PHI_{f}|divorced] = ϕ_{a} * EE[τbar/(1+ρ(1-τbar))]


function stage3(d,τgrid,cprobs) # Ebar = E[1/(1+)]
    phi_m = @chain d begin
        @subset :DIV.==0 :AGE.<=17 :AGE.>=7
        groupby(:AGE)
        @combine :PHIm = mean(skipmissing(:tau_f)) :wght = sum(.!ismissing.(:tau_f))
    end
    phi_d = @chain d begin
        @subset :DIV.==1 :AGE.<=17 :AGE.>=7
        groupby(:AGE)
        @combine :PHId = mean(skipmissing(:tau_f)) 
    end
    m = @chain phi_d begin
        innerjoin(phi_m,on=:AGE)
        @transform :r = :PHId ./ :PHIm
        @combine :m=sum(:wght .* :r) / sum(:wght)
        _.m[1]
    end
    modratio(ρ,τgrid,p) = sum(((τgrid .+ τgrid.^2) ./ (1 .+ ρ.*(1 .- τgrid))).*p)
    res = optimize(ρ->(m - modratio(ρ,τgrid,cprobs))^2,0,10)
    ρ_est = res.minimizer
    return ρ_est
end

function update_Cτ(θ,τgrid,cprobs)
    (;Cτ,β,αΓ_τHa,ρ) = θ
    Cτ[:] .= 0.     
    Ec = dot(cprobs,log.(1 .+ ρ .* (1 .- τgrid)))
    for a in reverse(1:length(Cτ)-1)
        Cτ[a] = αΓ_τHa[a] *  Ec + β * Cτ[a+1]
    end
    return θ
end

function stage3(θ,d,τgrid,cprobs)
    ρ_est = stage3(d,τgrid,cprobs)
    θ = (;θ...,ρ = ρ_est)
    θ = update_Cτ(θ,τgrid,cprobs)
    return θ
end