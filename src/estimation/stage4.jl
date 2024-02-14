function data_moms(M,P)
    # E[PT / FT|married / divorced]
    m0 = @chain P begin
        @subset .!ismissing.(:M_hrs)
        groupby(:DIV)
        @combine :PT = mean((:M_hrs.>1000) .& (:M_hrs.<2000)) :FT = mean(skipmissing(:M_hrs.>=2000))
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

    m4 = @chain M begin
        @transform :tD = :YDIV .- :YMAR :tF = :YBIRTH .- :YMAR
        groupby(:edW)
        @combine :div = mean(:YDIV.<2010) :tD10 = mean(:tD.<10) :tF5 = mean(:tF.<5) :tF10=mean(:tF.<10)
        stack(_)
    end
    
    return [m0.value ; m1.value; m2.value; m4.value]
end

function model_moms(TD,TF,D,Lsim,edW)
    moms = zeros(22)
    moms[1:4] .= (mean(Lsim[.!D].==2), 
            mean(Lsim[.!D].==1), 
            mean(Lsim[D].==2), 
            mean(Lsim[D].==1))
    moms[5:7] .= (mean(TD.<x) for x in (5,10,15))
    moms[8:10] .= (mean(TF.<x) for x in (2,4,6))
    I = TD.<9998
    moms[11:14] .= (mean((TD[I] .- TF[I]).<x) for x in (0,5,10,15))
    
    @views moms[15:18] .= (mean(TD[.!edW].<9998), mean(TD[edW].<9998), mean(TD[.!edW].<10), mean(TD[edW].<10))
    @views moms[19:22] .= (mean(TF[.!edW].<5), mean(TF[edW].<5), mean(TF[.!edW].<10), mean(TF[edW].<10))
     
    return moms
end

# function model_moms(TD,TF,D,Lsim)
#     moms = zeros(14)
#     moms[1:4] .= (mean(Lsim[.!D].==2), 
#             mean(Lsim[.!D].==1), 
#             mean(Lsim[D].==2), 
#             mean(Lsim[D].==1))
#     moms[5:7] .= (mean(TD.<x) for x in (5,10,15))
#     moms[8:10] .= (mean(TF.<x) for x in (2,4,6))
#     I = TD.<9998
#     moms[11:14] .= (mean((TD[I] .- TF[I]).<x) for x in (0,5,10,15))
#     return moms
# end



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
    TD,TF,L_sim, Ω_sim, D = data_gen(mod,dat,legal)
    return model_moms(TD,TF,D,L_sim,dat.edW.==1)
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