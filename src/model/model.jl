using LinearAlgebra, Interpolations, SparseArrays, Distributions
include("transitions.jl")

#-- define the parameters that won't change through the problem
function FixedParams(;
    β = 0.95,  # Discounting factor
    b = 16.75,         # Lower bound of consumption (consumption floor)
    A_W = 18:59,        # Wife age grid
    N_t = length(A_W),   # Length of wife age grid
    A_d = -10:10:10,     # Husband-wife age difference grid
    N_d = length(A_d),   # Length of grid of husband-wife age difference
    N_κ = 5 ,            # Length of wife human capital (experience) grid for each age
    κ_W  = hcat(collect.(range.(0, A_W .- 18 , length=N_κ))...)',
    κ_W_grid = push!(range.(0, A_W .- 18 .+10e-6 , length=N_κ),range(0,10e-6,length=N_κ)), # Wife human capital grid (N_t x N_κ)
    A_bar = 19 ,  # The age at which the child is fully developed
    A = collect(0:A_bar-1), # Child age
    A_grid = LinRange(0, 18, 4), # Child age grid
    N_a = length(A_grid), # Length of child age grid
    N_ϵ = 5,  # Length of husband income shock grid
    π_H = 1/3,# When divorced with a developing child, husband's income is taxed at rate π_H
    N_ω = 5,   # Length of marriage quality grid
    ω_grid = LinRange(-1,1,N_ω), # Marriage quality grid
    τ_R = 1/2, # Default legal custody arragement
    T_f = 25,  # The perod at which the couple can no longer bear chidren
    τgrid = [0.1,0.3,0.5,0.7,0.9], #<- the grid of potential allocations
    N_τ = length(τgrid),
    N_c = 1) #<- number of custody regimes
    return (;β,b,A_W,N_t,A_d,N_d,N_κ,κ_W,κ_W_grid,A_bar,A,A_grid,N_a,N_ϵ,π_H,N_ω,ω_grid,τ_R,T_f,τgrid,N_τ,N_c)
end

# - return a named tuple of parameters with default values
function Params(F)
    (;N_ϵ, N_ω, A, A_bar) = F
    ρ_ϵ = 0.9
    σ_η = 1.
    Λ,Π = rouwenhorst(N_ϵ ,σ_η, ρ_ϵ)
    π_ω = 0.8
    Π_ω = transmat_ω(π_ω,N_ω)
    return (;
    β = F.β , # inherit the discount factor
    α_C = 1. ,    # Consumption coefficient
    α_l = 0.5 ,   # Leisure coefficient
    γ_YW = fill(0.5,2,3), # Coefficients in wife's income process, Y_Wt = γ_YW[1] + γ_YW[2]*A_W + γ_YW[3]*κ_W
    σ_L = 1.          ,  # Scale parameter of wife labor shock, ϵ_L ∼ Gumbel(0,σ_L)
    γ_YH = fill(0.5,2,3) , # Coefficients in husband's income process, Y_Ht = γ_YH[1] + γ_YH[2]*A_H + ϵ
    # Husband income shock ϵ_t+1 = ρ_ϵ*ϵ_t + η_t+1, η ∼ N(0,σ_η^2)
    ρ_ϵ = 0.9 ,
    σ_η = 1. ,
    αΓ_τWa = zeros(A_bar) ,
    αΓ_τHa = zeros(A_bar) ,
    α_ω = [1., 1.5] ,  # Coefficients in marriage utility, α_ω[1] + α_ω[2]*ω
    σ_ω = 1.2 ,    # Scale parameter of marriage shock, ϵ_ω ∼ Gumbel(0,σ_ω)
    π_ω = π_ω ,    # Probability of staying in the same marriage quality

    # Residual terms coming from child skill production coefficients: ν_0 + ν_1 * ω
    α_νH0 = fill(0.3,2) ,  # Assume α_kH * ν_0 = α_νH0[1] + α_νH0[2]*a & α_kS * ν_0 = (α_kW/α_kH) * (α_νH0[1] + α_νH0[2]*a)
    α_νH1 = fill(0.3,2) , # Assume α_kH * ν_1 = α_νH1[1] + α_νH1[2]*a & α_kS * ν_1 = (α_kW/α_kH) * (α_νH1[1] + α_νH1[2]*a)
    α_τ = 5. , # First term in one-time payoff of custody decision
    π_lc = [0.13, 0.29] ,

    Cτ = fill(0.5,A_bar) ,
    ρ = 0.7   ,    # Additional marginal time cost of investing in the child.
    σ_τ = 1.   ,   # Scale parameter of the custody shock, ϵ_τ ∼ Gumbel(0,σ_τ)
    α_F = 120. ,
    σ_F = 2.    ,  # Scale parameter of the fertility shock, ϵ_F ∼ Gumbel(0,σ_F)
    Λ_ϵ = Λ, 
    Π_ϵ = Π, 
    Π_ω = Π_ω
    )
end
# -- set up storage for value functions
function values(F);
    (;N_t,N_κ,N_ϵ,N_d,N_a,N_ω,T_f) = F
    return (;
    # pre-allocated value function arrays
    # Stage 5: education group x ...
    VW5   = zeros(2,N_κ, N_t+1),
    pL5   = zeros(2,2,N_κ,N_t),
    VH5   = zeros(2,N_d, N_ϵ,N_t + 1),
    # Stage 4: education group x ...
    VW4   = zeros(2,2,N_d, N_ϵ, N_κ, N_a, N_t+1),
    VH4   = zeros(2,2,N_d, N_ϵ, N_κ, N_a, N_t+1),
    pL4   = zeros(2,2,2, N_d, N_ϵ, N_κ, N_a, N_t),
    # Stage 3: divorce law (mutual/unilateral) x education group W x education group H ...
    VW3   = zeros(N_κ, 2 * 2 * 2 *  N_d * N_ϵ * N_ω, N_t+1),
    VH3   = zeros(N_κ, 2 * 2 * 2 *  N_d * N_ϵ * N_ω, N_t+1),
    pL3  = zeros(2,N_κ, 2 * 2 * 2 *  N_d * N_ϵ * N_ω, N_t),
    pD3 = zeros(N_κ, 2 * 2 * 2 *  N_d * N_ϵ * N_ω, N_t),
    # Stage 2: divorce law (mutual/unilateral) x education group W x education group H ...
    VW2  = zeros(N_κ,N_a, 2 * 2 * 2 *  N_d * N_ϵ * N_ω , N_t+1),
    VH2  = zeros(N_κ,N_a, 2 * 2 * 2 *  N_d * N_ϵ * N_ω , N_t+1),
    pL2 = zeros(2,N_κ,N_a, 2 * 2 * 2 *  N_d * N_ϵ * N_ω , N_t+1),
    pD2 = zeros(N_κ,N_a, 2 * 2 * 2 *  N_d * N_ϵ * N_ω , N_t),
    # Stage 1: divorce law (mutual/unilateral) x education group W x education group H ...
    VW1   = zeros(N_κ, 2 * 2 * 2 *  N_d * N_ϵ * N_ω, T_f),
    VH1   = zeros(N_κ, 2 * 2 * 2 *  N_d * N_ϵ * N_ω, T_f),
    pL1  = zeros(2,N_κ, 2 * 2 * 2 *  N_d * N_ϵ * N_ω, T_f-1),
    pD1 = zeros(N_κ, 2 * 2 * 2 *  N_d * N_ϵ * N_ω, T_f-1),
    pF = zeros(N_κ, 2 * 2 * 2 *  N_d * N_ϵ * N_ω, T_f-1),

    # Expected custody cost
    Cτ = zeros(N_a) )
end

include("utility.jl")
include("interpolate.jl")
include("solve.jl")
include("simulate.jl")