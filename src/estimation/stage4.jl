function stage4(x0,θ,V,F,data,panel; R=10, num_iter=1000, show_trace = true)
    dat = prep_sim_data(data,panel;R = 10)
    moms0 = data_moms(data,panel)

    # this produces the best outcome from all the guesses
    res1 = optimize(x->ssq(update(x,θ,F),V,F,moms0,dat),x0,Optim.Options(iterations=num_iter,show_trace=show_trace))

    writedlm("output/stage4",res1.minimizer)

    θ = update(res1.minimizer,θ,F)
    return θ
end


function data_moms(M,P)
    # E[PT / FT|married / divorced]
    m0 = @chain P begin
        @subset .!ismissing.(:M_hrs)
        groupby(:DIV)
        @combine :PT = mean((:M_hrs.>0) .& (:M_hrs.<1500)) :FT = mean(skipmissing(:M_hrs.>=1500))
        stack(_) 
        @orderby :DIV :variable
    end
    
    # P[tD<x],for x = 5, 10, 15,
    # P[tF<x] for x = 2, 4, 6
    m1 = @chain M begin
        @transform :tD = :YDIV .- :YMAR :tF = :YBIRTH .- :YMAR
        @combine begin
            :tD5 = mean(:tD.<5)
            :tD10 = mean(:tD.<10)
            :tD15 = mean(:tD.<15)
            :tF2 = mean(:tF.<2)
            :tF4 = mean(:tF.<4)
            :tF6 = mean(:tF.<6)
        end
        stack(_)
    end

    # P[tD-tF<x|yD<=2010] for x in 0,5,10,15
    m2 = @chain M begin
        @subset :YDIV.<2010
        @transform :tDF = :YDIV .- :YBIRTH
        @combine begin
            :tDF0 = mean(:tDF.<0)
            :tDF5 = mean(:tDF.<5)
            :tDF10 = mean(:tDF.<10)
            :tDF15 = mean(:tDF.<15)
        end
        stack(_)
    end
    
    return [m0.value ; m1.value; m2.value] 
end

function model_moms(TD,TF,D,Lsim)
    moms = zeros(14) 
    moms[1:4] .= (mean(Lsim[.!D].==2), 
            mean(Lsim[.!D].==1), 
            mean(Lsim[D].==2), 
            mean(Lsim[D].==1))
    moms[5:7] .= (mean(TD.<x) for x in (5,10,15))
    moms[8:10] .= (mean(TF.<x) for x in (2,4,6))
    I = TD.<9998
    moms[11:14] .= (mean((TD[I] .- TF[I]).<x) for x in (0,5,10,15))
         
    return moms
end



function update(x,θ,F)
    π_ω = 1/(1+exp(x[12]))
    Π_ω = transmat_ω(π_ω, F.N_ω)
    return (;θ..., α_l=exp(x[1]), σ_L=exp(x[2]),α_νd=x[3], σ_F=exp(x[4]), α_ω=[x[5],exp(x[6])], α_νH0=x[7:8], α_νH1=x[9:10], σ_ω=exp(x[11]), π_ω, Π_ω)
end


function get_x(θ)
    return [log(θ.α_l);log(θ.σ_L);θ.α_νd;log(θ.σ_F);[θ.α_ω[1],log(θ.α_ω[2])];θ.α_νH0;θ.α_νH1;log(θ.σ_ω);log(1/θ.π_ω - 1)]
end

function update(x,ii,θ,F)
    xnew = get_x(θ)
    xnew[ii] .= x
    return update(xnew,θ,F)
end

# let's fix these?
function get_moments(θ,values,F,dat)
    mod = (;θ,values,F)
    solve_all!(mod)
    TD,TF,L_sim, _, D = data_gen(mod,dat)
    return model_moms(TD,TF,D,L_sim)
end

function ssq(θ,values,F,moms0,dat)
    moms1 = get_moments(θ,values,F,dat)
    ssq_ = sum((moms1 .- moms0).^2)
    if isnan(ssq_)
        #println(Inf)
        return Inf
    else
        #println(ssq_)
        return ssq_
    end
end

function ssq(θ,wght,values,F,moms0,dat)
    moms1 = get_moments(θ,values,F,dat)
    ssq_ = sum(wght .* (moms1 .- moms0).^2)
    if isnan(ssq_)
        #println(Inf)
        return Inf
    else
        #println(ssq_)
        return ssq_
    end
end

function estimate_blocks(θ,blocks,wght,values,F,moms0,dat)
    bi = 1
    for b in blocks
        println("Doing block $bi ...")
        x0 = get_x(θ)[b]
        r = optimize(x->ssq(update(x,b,θ,F),wght,values,F,moms0,dat),x0,Optim.Options(iterations=50))
        θ = update(r.minimizer,b,θ,F)
        bi += 1
    end
    return θ
end

function get_tsd(TD,dat)
    tsd = Int64[]
    for n in eachindex(TD)
        tlength = dat.tlength[n]
        for t in 1:tlength
            push!(tsd,t-1-TD[n])
        end
    end
    return tsd
end