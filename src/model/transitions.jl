
# Discretize AR(1)
# Grid and Markov transition matrix for husband income process shock
function rouwenhorst(n,σ,ρ)
    Λ = range(-σ*sqrt(n-1), σ*sqrt(n-1), length=n)
    p = (1+ρ)/2
    P_old = [p 1-p; 1-p p]
    for i in 3:n
        P_new = p.*[P_old zeros(i-1); zeros(i-1)' 0.] + (1-p).*[zeros(i-1) P_old; 0. zeros(i-1)'] + (1-p).*[zeros(i-1)' 0.; P_old zeros(i-1)] + p.*[0. zeros(i-1)'; zeros(i-1) P_old]
        P_new[2:i-1,:] = P_new[2:i-1,:]./2
        P_old = P_new
    end
    return (Λ=Λ,Π=P_old)
end

# might want to improve this one?
# or we might not use at all?
transmat_ω(π_ω,N_ω) = sparse(Matrix(π_ω*I, N_ω, N_ω)) + sparse([collect(1:N_ω-1)' N_ω]',[collect(2:N_ω)' N_ω]',[fill((1-π_ω)/2,N_ω-1)' 0]') + sparse([collect(2:N_ω)' N_ω]',[collect(1:N_ω-1)' N_ω]',[fill((1-π_ω)/2,N_ω-1)' 0]') + sparse([1,N_ω],[1,N_ω],[(1-π_ω)/2,(1-π_ω)/2]);
