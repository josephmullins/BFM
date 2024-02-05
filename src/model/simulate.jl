function draw_L(u,p1,p2)
    if u<p1
        return 1
    elseif u<(p1+p2)
        return 2
    else
        return 0
    end
end

function data_gen(mod,dat;seed=1234)
    (;values,θ,F) = mod
    (;Π_ϵ,Π_ω) = θ
    (;N_ϵ,N_t,N_ω,A_bar,N_d,T_f) = F
    # = MersenneTwister(seed)
    Random.seed!(seed)
    # interpolate choice probabilities
    pL1,pL2,pL3,pL4,pL5,pD1,pD2,pD3,pF = interpolate_probs(values,F)
    idx = LinearIndices((N_ϵ,N_ω,N_d,2,2,2))

    πϵ = [Categorical(Π_ϵ[ϵi,:]) for ϵi=1:N_ϵ]
    πω = [Categorical(Π_ω[ωi,:]) for ωi=1:N_ω]
    Lm = 0
    Nm = 0
    Ld = 0
    Nd = 0
    N = length(dat.eW)
    TD = zeros(N)
    TF = zeros(N)
    L_sim = zeros(N,N_t)
    Ω_sim = zeros(Int64,N,N_t)

    for n in 1:N
        eW = dat.eW[n]
        eH = dat.eW[n]
        ad = dat.AD[n]
        AW0 = max(dat.AW0[n],18)
        maxT = min(dat.maxT[n],60-AW0)
        κ = dat.κ0[n]
        l = dat.L[n]
        ωt = N_ω
        ϵt = 3 #<- TODO: we need to make this a draw from the steady state
        #ϵt = rand(πϵ[3]) # a random draw of ϵt given ϵ(t-1)=0
        tM = 1
        tF = 9999
        tD = 9999
        stage = 1
        AK = -1
        for t in 1:maxT
            Ω_sim[n,t] = ωt
            tt = AW0+t-18
            #println(n," ",stage," ",tt)
            if stage==1
                # birth decisions
                #println("hi")
                i = idx[ϵt,ωt,ad,eH,eW,l]
                if rand()<pF[tt](κ,i)
                    stage=2
                    AK = 0
                    tF = t
                else
                    pD = pD1[tt](κ,i)
                    if rand()<pD
                        stage=5
                        tD = t
                    else
                        p1 = pL1[tt](1,κ,i)
                        p2 = pL1[tt](2,κ,i)
                        L = draw_L(rand(),p1,p2)
                        L_sim[n,t] = L
                        Nm += 1
                        Lm += L
                        κ = κ + 0.5L
                    end
                end
            end
            if stage==2
                i = idx[ϵt,ωt,ad,eH,eW,l]
                pD = pD2[tt](κ,AK,i)
                if rand()<pD
                    stage=4
                    tD = t
                else
                    p1 = pL2[tt](1,κ,AK,i)
                    p2 = pL2[tt](2,κ,AK,i)
                    L = draw_L(rand(),p1,p2)
                    L_sim[n,t] = L
                    Nm += 1
                    Lm += L
                    κ = κ + 0.5L
                end
            end
            if stage==3
                i = idx[ϵt,ωt,ad,eH,eW,l]
                pD = pD3[tt](κ,i)
                if rand()<pD
                    stage=5
                    tD = t
                else
                    p1 = pL3[tt](1,κ,i)
                    p2 = pL3[tt](2,κ,i)
                    L = draw_L(rand(),p1,p2)
                    L_sim[n,t] = L
                    Nm += 1
                    Lm += L
                    κ = κ + 0.5L
                end
            end
            if stage==4
                p1 = pL4[tt](1,eW,eH,ad,ϵt,κ,AK)
                p2 = pL4[tt](2,eW,eH,ad,ϵt,κ,AK)
                L = draw_L(rand(),p1,p2)
                L_sim[n,t] = L
                Nm += 1
                Lm += L
                κ = κ + 0.5L
            end
            if stage==5
                p1 = pL4[tt](1,eW,κ)
                p2 = pL4[tt](2,eW,κ)
                L = draw_L(rand(),p1,p2)
                L_sim[n,t] = L
                Nm += 1
                Lm += L
                κ = κ + 0.5L
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
        end
        TD[n] = tD-1
        TF[n] = tF-1
        #println(tF," ",tD)
    end
    return Ld,Lm,TD,TF,L_sim, Ω_sim
end
