function draw_L(u,p1,p2)
    if u<p1
        return 1
    elseif u<(p1+p2)
        return 2
    else
        return 0
    end
end

function data_gen(mod,dat,legal::Vector{Int64};seed=1234)
    (;values,θ,F) = mod
    (;Π_ϵ,Π_ω) = θ
    (;N_ϵ,N_t,N_ω,A_bar,N_d,T_f) = F
    # = MersenneTwister(seed)
    Random.seed!(seed)
    # interpolate choice probabilities
    pL1,pL2,pL3,pL4,pL5,pD1,pD2,pD3,pF = interpolate_probs(values,F)
    idx = LinearIndices((N_ϵ,N_ω,N_d,2,2,2))
    πϵ = [Categorical(Π_ϵ[ϵi,:]) for ϵi=1:N_ϵ]
    πϵ0 = eigvecs(I(N_ϵ) .- Π_ϵ')[:,1]
    πϵ0 ./= sum(πϵ0)
    πϵ_init = Categorical(πϵ0)

    πω = [Categorical(Π_ω[ωi,:]) for ωi=1:N_ω]
    πω0 = Categorical(fill(1/N_ω,N_ω))
    N = length(dat.edW)
    TD = zeros(Int64,N)
    TF = zeros(Int64,N)
    NT = size(legal,1)
    L_sim = zeros(NT)
    Ω_sim = zeros(NT)
    D = zeros(Bool,NT)
    nt = 1
    for n in 1:N
        eW = 1+dat.edW[n]
        eH = 1+dat.edW[n]
        ad = dat.YBIRTH_M[n] - dat.YBIRTH_F[n]
        adi = 1 + (ad>-5) + (ad>5)
        AW0 = max(dat.YMAR[n]-dat.YBIRTH_M[n],18)
        maxT = dat.tlength[n]
        κ = dat.exp0[n]
        #ωt = N_ω #<- everyone starts at the highest draw
        ωt = rand(πω0) #<- try this also
        ϵt = rand(πϵ_init) #<- TODO: we need to make this a draw from the steady state
        #ϵt = rand(πϵ[3]) # a random draw of ϵt given ϵ(t-1)=0
        tM = 1
        tF = 9999
        tD = 9999
        stage = 1
        AK = -1
        for t in 1:maxT
            l = legal[nt]
            Ω_sim[nt] = ωt
            tt = AW0+t-18
            #println(n," ",stage," ",tt)
            if stage==1
                # birth decisions
                #println("hi")
                i = idx[ϵt,ωt,adi,eH,eW,l]
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
                        L_sim[nt] = L
                        κ = κ + L==2
                    end
                end
            end
            if stage==2
                i = idx[ϵt,ωt,adi,eH,eW,l]
                pD = pD2[tt](κ,AK,i)
                if rand()<pD
                    stage=4
                    tD = t
                else
                    p1 = pL2[tt](1,κ,AK,i)
                    p2 = pL2[tt](2,κ,AK,i)
                    L = draw_L(rand(),p1,p2)
                    L_sim[nt] = L
                    
                
                    κ = κ + L==2
                end
            end
            if stage==3
                i = idx[ϵt,ωt,adi,eH,eW,l]
                pD = pD3[tt](κ,i)
                if rand()<pD
                    stage=5
                    tD = t
                else
                    p1 = pL3[tt](1,κ,i)
                    p2 = pL3[tt](2,κ,i)
                    L = draw_L(rand(),p1,p2)
                    L_sim[nt] = L
                    
                
                    κ = κ + L==2
                end
            end
            if stage==4
                p1 = pL4[tt](1,eW,eH,ad,ϵt,κ,AK)
                p2 = pL4[tt](2,eW,eH,ad,ϵt,κ,AK)
                L = draw_L(rand(),p1,p2)
                L_sim[nt] = L
                
            
                κ = κ + L==2
            end
            if stage==5
                p1 = pL5[tt](1,eW,κ)
                p2 = pL5[tt](2,eW,κ)
                L = draw_L(rand(),p1,p2)
                L_sim[nt] = L
                
            
                κ = κ + L==2
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
            D[nt] = stage>=4
            ϵt = rand(πϵ[ϵt])
            ωt = rand(πω[ωt])
            nt += 1
        end
        TD[n] = tD-1
        TF[n] = tF-1
        #println(tF," ",tD)
    end
    return TD,TF,L_sim, Ω_sim, D
end
