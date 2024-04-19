function stage3(d,τgrid,cprobs)
    m = @chain d begin
        @subset .!ismissing.(:tau_f) :AGE.<=17
        groupby(:AGE)
        @transform :gamma_a = mean(:tau_f[.!:DIV])
        reg(@formula(tau_f ~ gamma_a + gamma_a&DIV + fe(KID)))
        coef(_)[2]
    end
    modratio(ρ,τgrid,p) = dot(log.(1 ./ (1 .+ ρ.*(1 .- τgrid))),p)
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