function data_moms(P)
    # E[PT / FT|married / divorced]
    m = @chain P begin
        @transform :F = :YEAR .>= :YBIRTH
        @subset :M_AGE.>=22 :M_AGE.<=40
        groupby(:M_AGE)
        @combine begin 
            :PT = mean(skipmissing((:M_hrs.>0) .& (:M_hrs.<1750)))
            :FT = mean(skipmissing(:M_hrs.>=1750))
            :DIV = mean(:DIV)
            :F = mean(:F)
            :DF = mean(:DIV .* :F)
        end
        stack(_)
    end
end

function model_moms(TD,TF,Lsim,tlength,age_mar)
    moms = zeros(19,5)
    denom = zeros(19)
    it = 1
    for i in eachindex(TD)
        am = age_mar[i]
        for t in 1:tlength[i]
            age_i = am+t - 22
            if age_i>=1 && age_i<=19
                d = t>=TD[i]+1
                f = t>=TF[i]+1
                ft = Lsim[it]==2
                pt = Lsim[it]==1
                denom[age_i] += 1
                moms[age_i,:] .+= (pt,ft,d,f,d*f)
            end
            it += 1
        end
    end
    moms ./= denom
     
    return moms[:]
end

function update(x,θ,F)
    π_ω = 1/(1+exp(x[12]))
    Π_ω = transmat_ω(π_ω, F.N_ω)
    return (;θ..., α_l=exp(x[1]), σ_L=exp(x[2]),α_F=x[3], σ_F=exp(x[4]), α_ω=[x[5],exp(x[6])], α_νH0=x[7:8], α_νH1=x[9:10], σ_ω=exp(x[11]), π_ω, Π_ω)
end


function get_x(θ)
    return [log(θ.α_l);log(θ.σ_L);θ.α_F;log(θ.σ_F);[θ.α_ω[1],log(θ.α_ω[2])];θ.α_νH0;θ.α_νH1;log(θ.σ_ω);log(1/θ.π_ω - 1)]
end

function update(x,ii,θ,F)
    xnew = get_x(θ)
    xnew[ii] .= x
    return update(xnew,θ,F)
end

# let's fix these?
function get_moments(θ,values,F,dat,legal)
    mod = (;θ,values,F)
    solve_all!(mod)
    TD,TF,L_sim, _, _ = data_gen(mod,dat,legal)
    return model_moms(TD,TF,L_sim,dat.tlength,dat.age_mar)
end

function ssq(θ,values,F,moms0,dat,legal)
    moms1 = get_moments(θ,values,F,dat,legal)
    ssq_ = sum((moms1 .- moms0).^2)
    if isnan(ssq_)
        #println(Inf)
        return Inf
    else
        #println(ssq_)
        return ssq_
    end
end

function ssq(θ,wght,values,F,moms0,dat,legal)
    moms1 = get_moments(θ,values,F,dat,legal)
    ssq_ = sum(wght .* (moms1 .- moms0).^2)
    if isnan(ssq_)
        #println(Inf)
        return Inf
    else
        #println(ssq_)
        return ssq_
    end
end

function estimate_blocks(θ,blocks,wght,values,F,moms0,dat,legal)
    bi = 1
    for b in blocks
        println("Doing block $bi ...")
        x0 = get_x(θ)[b]
        r = optimize(x->ssq(update(x,b,θ,F),wght,values,F,moms0,dat,legal),x0,Optim.Options(iterations=50))
        θ = update(r.minimizer,b,θ,F)
        bi += 1
    end
    return θ
end