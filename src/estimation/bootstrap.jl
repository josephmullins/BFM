function draw_boot_sample(M,P,K)
    N = nrow(M)
    Mb = M[rand(1:N,N),:]
    Mb[!,:MIDb] = 1:N

    Pb = @chain Mb begin
        @select :MID :MIDb
        innerjoin(P,on=:MID)
        @orderby :MIDb :YEAR
    end
    Kb = @chain Mb begin
        @select :MID :MIDb
        innerjoin(K,on=:MID)
        @orderby :MIDb :KID :YEAR
    end
    return Mb,Pb,Kb
end

function bootstrap_pars(θ,θk,V,F,M,P,K,num_boot ; R = 10, num_iter = 500, show_trace = false)
    seed0 = 1010
    seed1 = 2020
    seed2 = 3030
    x4_0 = get_x(θ)
    x5_0 = get_xk(θk)
    x1,x2,x3,x4,x5 = stack_ests(θ,θk)
    X1b = zeros(length(x1),num_boot)
    X2b = zeros(length(x2),num_boot)
    X3b = zeros(length(x3),num_boot)
    X4b = zeros(length(x4),num_boot)
    X5b = zeros(length(x5),num_boot)
    for b in 1:num_boot
        println("Performing bootstrap trial $b")
        Random.seed!(seed0+b)
        Mb,Pb,Kb = draw_boot_sample(M,P,K)
        seed4 = seed1+b
        seed5 = seed2+b
        θ,θk = estimate_model(x4_0, x5_0, θ, θk, V, F, Mb, Pb, Kb ; R, num_iter, show_trace, seed4, seed5)
        X1b[:,b],X2b[:,b],X3b[:,b],X4b[:,b],X5b[:,b] = stack_ests(θ,θk)
    end
    return X1b,X2b,X3b,X4b,X5b    
end
