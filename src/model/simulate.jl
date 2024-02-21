function draw_L(u,p1,p2)
    if u<p1
        return 1
    elseif u<(p1+p2)
        return 2
    else
        return 0
    end
end

# this function takes data from the dataframe and wraps it as an immutable tuple, repeating it R times where R is the number of simulations per observation
# - this supports type stability and avoids unnecessary allocations
function prep_sim_data(dat,panel;R = 1)
    dat = repeat(dat,R)
    legal = repeat(panel.L,R)

    (;edW = 1 .+ dat.edW,
    edH = 1 .+ dat.edH,
    tlength = dat.tlength,
    AW0 = max.(dat.YMAR .- dat.YBIRTH_M,18),
    AD = dat.YBIRTH_M .- dat.YBIRTH_F,
    κ = dat.exp0,
    legal)
end

# TODO: break this up into functions?
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
    πϵ0 = eigvecs(I(N_ϵ) .- Π_ϵ')[:,1]
    πϵ0 ./= sum(πϵ0)
    πϵ_init = Categorical(πϵ0)

    πω = [Categorical(Π_ω[ωi,:]) for ωi=1:N_ω]
    πω0 = Categorical(fill(1/N_ω,N_ω))
    N = size(dat.edW,1)
    TD = zeros(Int64,N)
    TF = zeros(Int64,N)
    NT = length(dat.legal)
    L_sim = zeros(NT)
    Ω_sim = zeros(NT)
    D = zeros(Bool,NT)
    nt = 1
    for n in 1:N
        eW = dat.edW[n]
        eH = dat.edW[n]
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
        for t in 1:maxT
            l = dat.legal[nt]
            Ω_sim[nt] = ωt
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
                        L_sim[nt] = L
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
                    L_sim[nt] = L
                    
                
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
                    L_sim[nt] = L
                    
                
                    κ += L==2
                end
            end
            if stage==4
                p1 = pL4[tt](1,eW,eH,ad,ϵt,κ,AK)
                p2 = pL4[tt](2,eW,eH,ad,ϵt,κ,AK)
                L = draw_L(rand(),p1,p2)
                L_sim[nt] = L                
                κ += L==2
            end
            if stage==5
                p1 = pL5[tt](1,eW,κ)
                p2 = pL5[tt](2,eW,κ)
                L = draw_L(rand(),p1,p2)
                L_sim[nt] = L
                
            
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
            D[nt] = stage>=4
            ϵt = rand(πϵ[ϵt])
            ωt = rand(πω[ωt])
            nt += 1
        end
        TD[n] = ttD-1
        TF[n] = ttF-1
        #println(tF," ",tD)
    end
    return TD,TF,L_sim, Ω_sim, D
end

# a function to extract the child sample
function prep_child_data(sim_data,dat,cprobs;seed=131332)
    (_,TF,L_sim,Ω_sim,D) = sim_data
    Random.seed!(seed)
    AK = Int64[]
    ΩK = Int64[]
    DK = Int64[]
    LK = Int64[]
    TL = Int64[]
    Csim = Int64[]
    nt = 1
    cdist = Categorical(cprobs)
    for n in eachindex(TF)
        maxT = dat.tlength[n]
        if TF[n]<9998
            t0 = nt+TF[n]
            t1 = min(t0+17,nt+maxT-1)
            tlength = min(18,t1-t0+1)
            ak = collect(0:(t1-t0))
            # now push the relevant data to the new vectors
            push!(TL,tlength)
            push!(AK,ak...)
            push!(ΩK,Ω_sim[t0:t1]...)
            push!(DK,D[t0:t1]...)
            push!(LK,L_sim[t0:t1]...)
            if TD[n]-TF[n]<18
                push!(Csim,rand(cdist))
            else
                push!(Csim,0)
            end
        end
        nt += maxT
    end
    return (;AK,DK,LK,TL,ΩK,Csim)
end


# a function to predict child outcomes (what's the initial condition?)
function predict_k!(TH,θk,θ,sim_data;seed=20240220)
    (;AK,DK,LK,TL,ΩK) = sim_data
    (;αΓ_τWa,αΓ_τHa) = θ
    nt = 1
    ϕW = αΓ_τWa ./ (1 .+ αΓ_τWa)
    ϕH = αΓ_τWa ./ (1 .+ αΓ_τWa)
    for n in eachindex(TL)
        k = γ_ψ0
        for t in 1:TL[n]    
            nonmarket_time = 112 - L_sim[nt]*40
            τW = ϕH[]
            if DK[nt]
                ψ = γ_ψ[1] + γ_ψ[4]*AK[nt]
            else
                ψ = γ_ψ[2] + γ_ψ[3]*ω_grid[ω] + γ_ψ*AK[nt]
            end
        end
    end
end

# a function to calculate moments (easy but let's check)

# a function to calculate the