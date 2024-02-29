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
    write(file," & \\multicolumn{2}{c}{Wife} & \\multicolumn{2}{c}{Wife} \\\\ \n ")
    write(file," & Non-College & College & Non-College & College \\\\ \\cmidrule(r){2-5} \n")
    write(file, " Const. &",tex_delimit(form.([θ.γ_YW[:,1] ; θ.γ_YH[:,1]]))," \\\\ \n")
    write(file," & ",tex_delimit(formse.([θse.γ_YW[:,1] ; θse.γ_YH[:,1]]))," \\\\ \n")
    write(file, " Age & & & ",tex_delimit(form3.(θ.γ_YH[:,2]))," \\\\ \n")
    write(file," & & & ",tex_delimit(form3se.(θse.γ_YH[:,2]))," \\\\ \n")
    write(file, " Age \$^2\$ & & & ",tex_delimit(form3.(θ.γ_YH[:,3]))," \\\\ \n")
    write(file," & & & ",tex_delimit(form3se.(θse.γ_YH[:,3]))," \\\\ \n")
    write(file, "FT Exp (\$\\kappa\$) & ",tex_delimit(form3.(θ.γ_YW[:,2]))," & & \\\\ \n")
    write(file," & ",tex_delimit(form3se.(θse.γ_YW[:,2]))," & & \\\\ \n")
    write(file, "FT Exp \$^2\$ (\$\\kappa^2\$) & ",tex_delimit(form3.(θ.γ_YW[:,3]))," & \\\\ \n")
    write(file," & ",tex_delimit(form3se.(θse.γ_YW[:,3]))," & &  \\\\ \n")
    close(file)
end

# tables remaining:
# - everything from stage 4 
# - model fit from stage 4
# - production parameters with ρ as well.
# - model fit from stage 4 plus picture of test scores.
# - also look at averages of each group relative to never married