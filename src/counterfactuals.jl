

function counterfactual_statistics(kid_data,sim_data,dat,θ,θk,F,mod)
    # Ex-Ante Welfare

    # Child Skills

    # Mother's Log Potential Wages
end

function average_welfare()
    (;values,θ,F)
    (;VW1,VH1,VW2,VH2,VW3,VH3,pL1,pD1,pF,VW5,VH5) = values

    #VW & VH : T_f-1 x N_d x N_ϵ x N_κ x N_ω
    etpW = interpolate3(VW1,F,t)
    etpH = interpolate3(VH1,F,t)

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
            V = 
            if stage==1
                # birth decisions
                #println("hi")
                i = idx[ϵt,ωt,adi,eH,eW,l]
            end
        end
    end
end