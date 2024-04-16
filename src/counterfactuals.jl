function compare_stats(stats1,stats0,β)
    ΔwH = exp((1-β) * (stats1.welf_H - stats0.welf_H)) - 1
    ΔwW = exp((1-β) * (stats1.welf_W - stats0.welf_W)) - 1
    Δlogwage = stats1.log_wage - stats0.log_wage
    Δfert = stats1.fertility - stats0.fertility
    Δdiv = stats1.divorce - stats0.divorce
    Δskill = (stats1.decomp .- stats0.decomp) / stats0.se
    return (;ΔwH,ΔwW,Δlogwage,Δfert,Δdiv,Δskill)
end

function divorce_standard_counterfactual(dat,mod,θk,stats0)
    (;θ) = mod
    (;cprobs, β) = θ

    # all mutual consent
    dat.legal[:] .= 1
    sim_data,kid_data = full_simulation(dat,mod,cprobs)
    stats_mc = counterfactual_statistics(kid_data,dat,θ,θk,mod)
    r1 = compare_stats(stats_mc,stats0,β)

    # all unilateral
    dat.legal[:] .= 2
    sim_data,kid_data = full_simulation(dat,mod,cprobs)
    stats_ul = counterfactual_statistics(kid_data,dat,θ,θk,mod)
    r2 = compare_stats(stats_ul,stats0,β)

    return r1,r2
end

function custody_counterfactual(dat,mod,θk)
    (;θ,F) = mod
    (;cprobs, β) = θ

    # sole maternal custody
    cprobs0 = [1.,0.,0.,0.,0.]
    θ = update_Cτ(θ,F.τgrid,cprobs0)
    mod = (;mod...,θ)
    sim_data,kid_data = full_simulation(dat,mod,cprobs0)
    stats0 = counterfactual_statistics(kid_data,dat,θ,θk,mod)

    # complete split
    cprobs1 = [0.,0.,0.,0.,1.]
    θ = update_Cτ(θ,F.τgrid,cprobs1)
    mod = (;mod...,θ)
    sim_data,kid_data = full_simulation(dat,mod,cprobs1)
    stats1 = counterfactual_statistics(kid_data,dat,θ,θk,mod)

    ΔwH = exp((1-β) * (stats1.welf_H - stats0.welf_H)) - 1
    ΔwW = exp((1-β) * (stats1.welf_W - stats0.welf_W)) - 1
    Δlogwage = stats1.log_wage - stats0.log_wage
    Δfert = stats1.fertility - stats0.fertility
    Δdiv = stats1.divorce - stats0.divorce
    Δskill = (stats1.decomp .- stats0.decomp) / stats0.se
    return (;ΔwH,ΔwW,Δlogwage,Δfert,Δdiv,Δskill)

end

function child_support_counterfactual(dat,mod,θk)
    (;θ, F) = mod
    (;cprobs, β) = θ
    # baseline
    sim_data,kid_data = full_simulation(dat,mod,cprobs)
    stats0 = counterfactual_statistics(kid_data,dat,θ,θk,mod)

    # increase child support to 50% (??) did I set this too high to begin with??
    F = (;F...,π_H = 0.5)
    mod = (;mod...,F)
    sim_data,kid_data = full_simulation(dat,mod,cprobs)
    stats1 = counterfactual_statistics(kid_data,dat,θ,θk,mod)

    ΔwH = exp((1-β) * (stats1.welf_H - stats0.welf_H)) - 1
    ΔwW = exp((1-β) * (stats1.welf_W - stats0.welf_W)) - 1
    Δlogwage = stats1.log_wage - stats0.log_wage
    Δfert = stats1.fertility - stats0.fertility
    Δdiv = stats1.divorce - stats0.divorce
    Δskill = (stats1.decomp .- stats0.decomp) / stats0.se
    return (;ΔwH,ΔwW,Δlogwage,Δfert,Δdiv,Δskill)
end


function counterfactual_statistics(kid_data,dat,θ,θk,mod)
    # Ex-Ante Welfare and 
    welf_H, welf_W, log_wage, fertility, divorce = model_stats(mod,dat)

    # Child Skills
    d, se = child_skill_outcomes(θk,θ,F,kid_data)
    (; welf_H,welf_W,log_wage,fertility,divorce,decomp = d, se)
end

function child_skill_outcomes(θk,θ,F,kid_data)
    inputs,_ = input_decomposition(θk,θ,F,kid_data)
    d = mean(inputs,dims=2)
    d_sum = sum(d,dims=1)
    sd_skills = std(sum(inputs,dims=1))
    d = [d_sum ; d]
    return d, sd_skills
end


# this function calculates:
# - initial welfare of both agents
# - final log-wages
# - final divorce rates
# - final fertility rates
function model_stats(mod,dat;seed=1234)
    (;values,θ,F) = mod
    (;Π_ϵ,Π_ω) = θ
    (;N_ϵ,N_t,N_ω,A_bar,N_d,T_f) = F
    (;VW1,VH1) = values
    (;γ_YW) = θ

    # = MersenneTwister(seed)
    Random.seed!(seed)

    # interpolate choice probabilities
    pL1,pL2,pL3,pL4,pL5,pD1,pD2,pD3,pF = interpolate_probs(values,F)
    
    # interpolate stage 1 values
    etpW = [interpolate3(VW1,F,t) for t in 1:24]
    etpH = [interpolate3(VH1,F,t) for t in 1:24]

    # create an indexing rule
    idx = LinearIndices((N_ϵ,N_ω,N_d,2,2,2))

    # create distribution objects (including initial)
    πϵ = [Categorical(Π_ϵ[ϵi,:]) for ϵi=1:N_ϵ]
    πϵ0 = eigvecs(I(N_ϵ) .- Π_ϵ')[:,1]
    πϵ0 ./= sum(πϵ0)
    πϵ_init = Categorical(πϵ0)
    πω = [Categorical(Π_ω[ωi,:]) for ωi=1:N_ω]
    πω0 = Categorical(fill(1/N_ω,N_ω))

    # 
    N = size(dat.edW,1)

    # welfare
    welf_H = 0.
    welf_W = 0.

    # mother's log-wage
    log_wage = 0.

    # fertility and divorce rate
    fertility = 0.
    divorce = 0.

    nt = 1
    for n in 1:N
        eW = dat.edW[n]
        eH = dat.edH[n]
        ad = dat.AD[n]
        adi = 1 + (ad>-5) + (ad>5)
        aw0 = dat.AW0[n]
        maxT = dat.tlength[n]
        κ = dat.κ[n]
        #ωt = N_ω #<- everyone starts at the highest draw (an alternative assumption)
        ωt = rand(πω0) 
        ϵt = rand(πϵ_init) 
        tM = 1
        ttF = 9999
        ttD = 9999
        stage = 1
        AK = -1

        # add welfare of wife and husband at beginning
        l = dat.legal[nt]
        i = idx[ϵt,ωt,adi,eH,eW,l]
        tt = aw0+1-18
        welf_H += etpH[tt](κ,i)
        welf_W += etpW[tt](κ,i)

        for t in 1:maxT
            l = dat.legal[nt]
            tt = aw0+t-18
            #println(n," ",stage," ",tt)
            if stage==1
                # birth decisions
                #println("hi")
                i = idx[ϵt,ωt,adi,eH,eW,l]
                if rand()<pF[tt](κ,i)
                    stage=2
                    AK = 0
                    ttF = t
                else
                    pD = pD1[tt](κ,i)
                    if rand()<pD
                        stage=5
                        ttD = t
                    else
                        p1 = pL1[tt](1,κ,i)
                        p2 = pL1[tt](2,κ,i)
                        L = draw_L(rand(),p1,p2)
                        κ += L==2
                    end
                end
            end
            if stage==2
                i = idx[ϵt,ωt,adi,eH,eW,l]
                pD = pD2[tt](κ,AK,i)
                if rand()<pD
                    stage=4
                    ttD = t
                else
                    p1 = pL2[tt](1,κ,AK,i)
                    p2 = pL2[tt](2,κ,AK,i)
                    L = draw_L(rand(),p1,p2)
                    κ += L==2
                end
            end
            if stage==3
                i = idx[ϵt,ωt,adi,eH,eW,l]
                pD = pD3[tt](κ,i)
                if rand()<pD
                    stage=5
                    ttD = t
                else
                    p1 = pL3[tt](1,κ,i)
                    p2 = pL3[tt](2,κ,i)
                    L = draw_L(rand(),p1,p2)
                    κ += L==2
                end
            end
            if stage==4
                p1 = pL4[tt](1,eW,eH,ad,ϵt,κ,AK)
                p2 = pL4[tt](2,eW,eH,ad,ϵt,κ,AK)
                L = draw_L(rand(),p1,p2)
                κ += L==2
            end
            if stage==5
                p1 = pL5[tt](1,eW,κ)
                p2 = pL5[tt](2,eW,κ)
                L = draw_L(rand(),p1,p2)            
                κ += L==2
            end
            # final updates on stages
            if (stage==1) & (tt==T_f-1)
                stage = 3
            end
            if stage==2
                if AK==A_bar-1
                    stage = 3
                else
                    AK += 1
                end
            end
            if stage==4
                if AK==A_bar-1
                    stage = 5
                else
                    AK += 1
                end
            end

            # draw new shocks
            ϵt = rand(πϵ[ϵt])
            ωt = rand(πω[ωt])
            nt += 1
        end
        log_wage += γ_YW[eW,1] + γ_YW[eW,2]*κ + γ_YW[eW,3]*κ^2
        fertility += ttF<9999
        divorce += ttD<9999
        #println(tF," ",tD)
    end
    return welf_H / N, welf_W / N, log_wage / N, fertility / N, divorce / N
end
