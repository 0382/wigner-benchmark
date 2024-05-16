using Primes

prime32 = primes(136)

for i = 2:136
    f = Dict(factor(i).pe)
    maxp = maximum(keys(f))
    v = Vector{Int}()
    for p in prime32
        if p > maxp
            break
        end
        push!(v, get(f, p, 0))
    end
    print("{", join(v, ", "), "}, // $i = ")
    display(factor(i))
    println()
end