function stage5(x0,θk,mod,dat,K;num_iter = 1000,show_trace = true, seed = 20240220)
    (;θ) = mod
    (;cprobs) = θ

    solve_all!(mod)
    sim_data = data_gen(mod,dat);

    kid_data = prep_child_data(sim_data,dat,cprobs);

    S = zeros(length(kid_data.ΩK))

    # the data moments:
    kmoms0 = kidmoms_data(K)

    #x0 = [-20.,0.75,0.75,0.05,0.00,2.5,log(0.2),log(0.7),log(20.)]

    res = optimize(x->obj_stage5(S,updateθk(x,θk,θ),θ,F,kid_data,kmoms0 ; seed),x0,Optim.Options(iterations = num_iter,show_trace = show_trace))
    θk =  updateθk(res.minimizer,θk,θ)
    
    return θk
end

function kidmoms_data(d)
    m = testscore_averages(d)
    ma = m.S
    wght = m.wght
    ma = testscore_averages(d).S
    sa = testscore_sd(d).sd
    k = construct_regression_data(d)
    mod = lm(@formula(AP_raw ~ log(tau_m) + AP_lag ),k)
    mb = coef(mod)[2:3]
    return ma,sa,mb,wght
end

function testscore_averages(d) #<-?
    d0 = DataFrame(AGE=repeat(4:15,3),dgroup = repeat(1:3,inner=12))
    m = @chain d begin
        @subset :AGE.<=15
        @subset :AGE.>=4
        #@transform :AGE = round.(:AGE ./ 4)
        groupby([:dgroup,:AGE])
        @combine :S = mean(skipmissing(:AP_raw))
        rightjoin(d0,on=[:dgroup,:AGE])
        @orderby :dgroup :AGE
        @transform :wght = 1. * .!ismissing.(:S) :S = coalesce.(:S,0.)
    end
    return m
end

function testscore_sd(d)
    m = @chain d begin
        @subset :AGE.<=15
        @subset :AGE.>=4
        #@transform :AGE = round.(:AGE ./ 4)
        groupby(:AGE)
        @combine :sd = std(skipmissing(:AP_raw))
    end
    return m
end


function construct_regression_data(d)
    # how does the regression look
    k1 = @chain d begin
        @subset :YEAR.==1997 .|| :YEAR.==2002
        @transform :g1 = :dgroup.==1 :g2 = :dgroup.==2
        select(:YEAR,:KID,:AP_raw,:LW_raw, :tau_m,:tau_f,:g1,:g2)
        @rename :AP_lag = :AP_raw :LW_lag = :LW_raw
    end
    k = @chain d begin
        @subset :YEAR.==2002 .|| :YEAR.==2007
        @transform :YEAR = :YEAR .- 5
        select(:KID,:YEAR,:AP_raw,:LW_raw,:AGE)
        innerjoin(k1,on=[:KID,:YEAR])
        @subset :tau_m.>0
        dropmissing()
    end
    return k
end


function kid_moments(S,θk,θ,F,kid_data ; seed = 20240220)
    (;Y,X) = kid_data
    predict_k!(S,θk,θ,F,kid_data ; seed)
    ma = [mean(S[kid_data.AK.==a .&& kid_data.G.==g]) for g in 1:3 for a in 4:15]
    sa = [std(S[kid_data.AK.==a]) for a in 4:15]
    β = inv(X' * X) * X' * Y
    mb = β[2:3]
    return ma,sa,mb
end

function obj_stage5(S,θk,θ,F,kid_data,kmoms0 ; seed = 20240220)
    ma0,sa0,mb0,wght = kmoms0
    ma1,sa1,mb1 = kid_moments(S,θk,θ,F,kid_data ; seed)
    return sum(wght .* (ma1 .- ma0).^2) + sum((sa1 .- sa0).^2) + 1000*sum((mb1 .- mb0).^2)
end

function updateθk(x,θk,θ)
    (;γ_ψ,δW,δH,Γa) = θk
    (;αΓ_τWa,αΓ_τHa,β) = θ
    γ_ψ0 = x[1]
    γ_ψ[:] .= x[2:5]
    δk = exp(x[6])/(1+exp(x[6]))
    get_Γa!(Γa,δk,β)
    δW[1] = exp(x[7])
    δH[1] = δW[1] * αΓ_τHa[1] / αΓ_τWa[1]
    #δH[1] = exp(x[8])
    for a in eachindex(δW)
        δW[a] = δW[1] * αΓ_τWa[a] * Γa[2] / (αΓ_τWa[1] * Γa[a+1])
        δH[a] = δH[1] * αΓ_τHa[a] * Γa[2] / (αΓ_τHa[1] * Γa[a+1])
    end
    σ_η = exp(x[8])
    γAP = exp(x[9])
    return (;θk...,γ_ψ0,δk,σ_η,γAP)
end

function get_xk(θk)
    (;γ_ψ0,γ_ψ,δk,δW,σ_η,γAP) = θk
    return [γ_ψ0 ; γ_ψ ; log(δk / (1-δk)) ; log(δW[1]) ; log(σ_η) ; log(γAP)]
end

function get_Γa!(Γa,δk,β)
    Γ_next = 1/(1-β) 
    for a in reverse(eachindex(Γa))
        Γa[a] = 1 + β * δk * Γ_next
        Γ_next = Γa[a]
    end
end

function input_decomposition(θk,θ,F,sim_data;seed=20240220)
    Random.seed!(seed)
    (;δW,δH,δk,γ_ψ0,γ_ψ) = θk
    (;τgrid,ω_grid) = F
    (;AK,DK,LK,TL,ΩK,Csim,reg_idx,Y,X,logτ_m) = sim_data
    (;αΓ_τWa,αΓ_τHa,ρ) = θ
    nt = 1
    ϕW = αΓ_τWa ./ (1 .+ αΓ_τWa)
    ϕH = αΓ_τHa ./ (1 .+ αΓ_τHa)
    inputs = zeros(3,length(TL))
    mstat = zeros(length(TL))
    for n in eachindex(TL)
        tfp = γ_ψ0
        mtime = 0.
        ftime = 0.
        for t in 1:TL[n]
            ai = AK[nt]+1 
            ωi = ΩK[nt]
            nonmarket_time = 112 - LK[nt]*40
            τW = ϕW[ai] * nonmarket_time
            τH = ϕH[ai] * 72
            if DK[nt]
                ψ = γ_ψ[1] + γ_ψ[4]*AK[nt]
                ci = Csim[n]
                time_penalty =  - δH[ai] * log(1 + ρ * (1-τgrid[ci]))
            else
                ψ = γ_ψ[2] + γ_ψ[3]*ω_grid[ωi] + γ_ψ[4]*AK[nt]
                time_penalty = 0.
            end
            tfp = ψ + δk*tfp
            mtime = δW[ai] * log(τW) + δk*mtime
            ftime = δH[ai]*log(τH) + time_penalty + δk*ftime
            nt += 1
        end
        inputs[:,n] .= (tfp,mtime,ftime)
        mstat[n] = DK[nt-1] 
    end
    return inputs,mstat
end
