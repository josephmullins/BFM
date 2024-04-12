using Printf

function stack_ests(θ,θk)
    # stage1
    x1 = [θ.γ_YH[:];θ.γ_YW[:];θ.ρ_ϵ;θ.σ_η] #(+ rouwenhorst)

    # stage 2
    x2 = [θ.αΓ_τHa;θ.αΓ_τWa]

    # stage 3
    x3 = [θ.ρ;θ.Cτ]

    # stage 4
    x4 = [θ.α_l;θ.σ_L;θ.α_νd;θ.σ_F;θ.α_ω;θ.α_νH0;θ.α_νH1;θ.σ_ω;θ.π_ω] # (+ transmat)

    # stage 5 (θk)
    x5 = [θk.γ_ψ0;θk.γ_ψ;θk.δH;θk.δW;θk.δk;θk.σ_η;θk.γAP] #(+ get_Γa!)

    return x1,x2,x3,x4,x5
end

function update_all(x,θ,θk,F)
    (x1,x2,x3,x4,x5) = x
    (;N_ϵ) = F
    # stage 1
    θ.γ_YH[:] .= x1[1:6] #<- this will edit the original θ?
    θ.γ_YW[:] .= x1[7:12]
    ρ_ϵ = x1[13]
    σ_η = x1[14]
    Λ_ϵ, Π_ϵ = rouwenhorst(N_ϵ,σ_η,ρ_ϵ)
    θ = (;θ...,ρ_ϵ,σ_η,Λ_ϵ,Π_ϵ)

    # stage 2
    αΓ_τHa = x2[1:19]
    αΓ_τWa = x2[20:38]
    θ = (;θ...,αΓ_τHa,αΓ_τWa)

    # stage 3
    ρ = x3[1]
    Cτ = x3[2:20]
    θ = (;θ...,ρ,Cτ)

    # stage 4:
    α_l,σ_L,α_νd,σ_F = x4[1:4]
    α_ω = x4[5:6]
    α_νH0 = x4[7:8]
    α_νH1 = x4[9:10]
    σ_ω,π_ω = x4[11:12]
    Π_ω = transmat_ω(π_ω, F.N_ω)
    θ = (;θ...,α_l,σ_L,α_νd,σ_F,α_ω,α_νH0,α_νH1,σ_ω,π_ω,Π_ω)

    # stage 5
    γ_ψ0 = x5[1]
    γ_ψ = x5[2:5]
    δH = x5[6:23]
    δW = x5[24:41]
    δk = x5[42]
    σ_η,γAP = x5[43:44] #(+ get_Γa!)
    get_Γa!(θk.Γa,δk,θ.β)
    θk = (;θk...,γ_ψ0,γ_ψ,δH,δW,δk,σ_η,γAP)

    return θ,θk
end

form(x) = @sprintf("%0.2f",x) #<- print to 2 digits
formse(x) = string("(",@sprintf("%0.2f",x),")")
form3(x) = @sprintf("%0.3f",x) #<- print to 2 digits
form3se(x) = string("(",@sprintf("%0.3f",x),")")
# a helper function to write a collection of strings into separate columns
function tex_delimit(x)
    str = x[1]
    num_col = length(x)
    for i in 2:num_col
        str *=  "&" * x[i]
    end
    return str
end

function write_inc_table(θ,θse)
    file = open("output/tables/income_ests.tex", "w")
    write(file," & \\multicolumn{2}{c}{Wife} & \\multicolumn{2}{c}{Husband} \\\\ \n ")
    write(file," & Non-College & College & Non-College & College \\\\ \\cmidrule(r){2-3}\\cmidrule(r){4-5} \n")
    write(file, " Const. &",tex_delimit(form.([θ.γ_YW[:,1] ; θ.γ_YH[:,1]]))," \\\\ \n")
    write(file," & ",tex_delimit(formse.([θse.γ_YW[:,1] ; θse.γ_YH[:,1]]))," \\\\ \n")
    write(file, "FT Exp (\$\\kappa\$) & ",tex_delimit(form3.(θ.γ_YW[:,2]))," & & \\\\ \n")
    write(file," & ",tex_delimit(form3se.(θse.γ_YW[:,2]))," & & \\\\ \n")
    write(file, "FT Exp \$^2\$ (\$\\kappa^2\$) & ",tex_delimit(form3.(θ.γ_YW[:,3]))," & \\\\ \n")
    write(file," & ",tex_delimit(form3se.(θse.γ_YW[:,3]))," & &  \\\\ \n")
    write(file, " Age & & & ",tex_delimit(form3.(θ.γ_YH[:,2]))," \\\\ \n")
    write(file," & & & ",tex_delimit(form3se.(θse.γ_YH[:,2]))," \\\\ \n")
    write(file, " Age \$^2\$ & & & ",tex_delimit(form3.(θ.γ_YH[:,3]))," \\\\ \n")
    write(file," & & & ",tex_delimit(form3se.(θse.γ_YH[:,3]))," \\\\ \n")
    write(file, "\$\\rho_{\\epsilon}\$ & & & ",form(θ.ρ_ϵ)," & ",form(θ.ρ_ϵ),"\\\\ \n")
    write(file, " & & & ",formse(θse.ρ_ϵ)," & ",formse(θse.ρ_ϵ),"\\\\ \n")
    write(file, "\$\\sigma_{\\eta}\$ & & & ",form(θ.σ_η)," & ",form(θ.σ_η),"\\\\ \n")
    write(file, " & & & ",formse(θse.σ_η)," & ",formse(θse.σ_η),"\\\\ \n")
    close(file)
end

# tables remaining:
# - everything from stage 4 
function write_prefs(θ,θse)
    # x4 = [θ.α_l;θ.σ_L;θ.α_νd;θ.σ_F;θ.α_ω;θ.α_νH0;θ.α_νH1;θ.σ_ω;θ.π_ω]
    str = ["\$\\alpha_l\$","\$\\sigma_{L}\$","\$\\alpha_{\\omega,0}\$","\$\\alpha_{\\omega,1}\$",
        "\$\\sigma_{\\omega}\$","\$\\pi_{\\omega}\$","\$\\alpha_{\\nu,0,0}\$","\$\\alpha_{\\nu,0,1}\$",
        "\$\\alpha_{\\nu,1,0}\$","\$\\alpha_{\\nu,1,1}\$","\$\\alpha_{\\nu,D}\$","\$\\sigma_{F}\$"]
    par = [θ.α_l ; θ.σ_L ; θ.α_ω ; θ.σ_ω ; θ.π_ω ; θ.α_νH0 ; θ.α_νH1 ; θ.α_νd ; θ.σ_F]
    se = [θse.α_l ; θse.σ_L ; θse.α_ω ; θse.σ_ω ; θse.π_ω ; θse.α_νH0 ; θse.α_νH1 ; θse.α_νd ; θse.σ_F]
    file = open("output/tables/stage_4.tex", "w")
    write(file,"Parameter & Estimate & (Std Error) \\\\ \\cmidrule(r){1-3} \n")
    for r in eachindex(str)
        write(file,str[r]," & ",form(par[r])," & ",formse(se[r])," \\\\ \n")
    end
    close(file)
end
# - model fit from stage 4
function write_modelfit_s4(moms0,moms1)
    mstring = ["FT - Married","PT - Married","FT - Divorced","PT - Divorced"]
    mstring = [mstring;["Divorce \$<\$ $x years" for x in (5,10,15)]]
    mstring = [mstring;["Birth \$<\$ $x years" for x in (2,4,6)]]
    mstring = [mstring;["Time to Divorce - Time to Birth \$<\$ $x years" for x in (0,5,10,15)]]
    file = open("output/tables/model_fit_s4.tex", "w")
    write(file,"Moment & Data & Model \\\\ \\cmidrule(r){1-3} \n")
    for r in eachindex(mstring)
        write(file,mstring[r]," & ",form(moms0[r])," & ",form(moms1[r])," \\\\ \n")
    end
    close(file)
end

function write_modelfit_s5(moms0,moms1)
    coefficient = ["\$\\log(\\text{mother's time})\$","Lagged Test Score"]
    file = open("output/tables/model_fit_s5.tex", "w")
    write(file,"Coefficient & Data & Model \\\\ \\cmidrule(r){1-3} \n")
    for r in eachindex(coefficient)
        write(file,coefficient[r]," & ",form(moms0[r])," & ",form(moms1[r])," \\\\ \n")
    end
    close(file)
end

# - production parameters with ρ as well.
function write_production(θk,θkse)
    #x5 = [θk.γ_ψ0;θk.γ_ψ;θk.δH;θk.δW;θk.δk;θk.σ_η;θk.γAP] #(+ get_Γa!)
    str = [["\$\\gamma_{\\psi,$x}\$" for x in 0:4] ; "\$\\delta_{k}\$" ; "\$\\sigma_{\\xi}\$" ; "\$\\gamma_{AP}\$" ; "\$\\rho\$"]
    par = [θk.γ_ψ0; θk.γ_ψ; θk.δk; θk.σ_η; θk.γAP ; θk.ρ]
    se = [θkse.γ_ψ0; θkse.γ_ψ; θkse.δk; θkse.σ_η; θkse.γAP ; θkse.ρ]
    file = open("output/tables/stage_5.tex", "w")
    write(file,"Parameter & Estimate & (Std Error) \\\\ \\cmidrule(r){1-3} \n")
    for r in eachindex(str)
        write(file,str[r]," & ",form(par[r])," & ",formse(se[r])," \\\\ \n")
    end
    close(file)
end
# - model fit from stage 5 plus picture of test scores.
# - also look at averages of each group relative to never divorced

function write_decomposition(d,dse)
    str = ["Total","TFP","Mother's Time","Father's Time"]
    file = open("output/tables/skill_decomposition.tex","w")
    write(file,"Input & Estimate & (Std Error) \\\\ \\cmidrule(r){1-3} \n")
    for r in eachindex(str)
        write(file,str[r]," & ",form(d[r])," & ",formse(dse[r])," \\\\ \n")
    end
    close(file)
end
