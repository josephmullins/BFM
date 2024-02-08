function stage2(d)
    phi = @chain d begin
        @subset :DIV.==0 :AGE.<=18
        groupby(:AGE)
        @combine begin
            :phi_f = mean(skipmissing(:PHI_F))
            :phi_m = mean(skipmissing(:PHI_M))
        end
        @orderby :AGE
        Matrix(_)[:,2:3]
    end
    αΓl = phi ./ (1 .- phi)
    return αΓl
end

function stage2(θ,d)
    αΓl = stage2(d)
    @views θ.αΓ_τHa[:] .= αΓl[:,1]
    @views θ.αΓ_τWa[:] .= αΓl[:,2]
    return θ
end