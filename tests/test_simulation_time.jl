using Random

# test 1: drawing 1000 binary variables
p = rand(10000)
function sim1(p)
    for n in eachindex(p)
        r = rand()<p[n]
    end
end

u = rand(10000)
function sim2(p,u)
    for n in eachindex(p)
        r = u[n]<p[n]
    end
end

@show "Test for the simple case"

sim1(p)
sim2(p,u)

@time sim1(p)
@time sim2(p,u)

# test 2: testing a choice probability simulation for a thousand agents
P = rand(5,20)
U = rand(2,20,100)
S = zeros(20,100)
function sim1!(S,P)
    s = 1
    for r in axes(S,2)
        for n in axes(S,1)
            S[n,r] = s
            choice = rand()<P[s,n]
            if choice
                up = rand()<0.1
            else
                up = false
            end
            s = min(5,s + up)
        end
    end
end

function sim1!(S,P,U)
    s = 1
    for r in axes(S,2)
        for n in axes(S,1)
            S[n,r] = s
            choice = U[1,n,r]<P[s,n]
            if choice
                up = U[2,n,r]<0.1
            else
                up = false
            end
            s = min(5,s + up)
        end
    end
end

@show "Test for the less simple case"
sim1!(S,P)
@time sim1!(S,P)
sim1!(S,P,U)
@time sim1!(S,P,U)