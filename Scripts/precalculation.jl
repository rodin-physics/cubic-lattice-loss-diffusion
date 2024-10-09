include("../src/main.jl")

# Precalculate W
function W_precalc(lim)
    displacements = Iterators.product([-lim:lim, -lim:lim, -lim:lim]...) |> collect |> vec
    Ws = Vector{Matrix{Float64}}(undef, length(displacements))
    pr = Progress(length(Ws))
    Threads.@threads for ii in eachindex(Ws)
        Ws[ii] = W(m, dyn_mat, 0, displacements[ii])[1]
        next!(pr)
    end
    return Dict(zip(displacements, Ws))
end

if !isfile("Data/W_precalc.jld2")
    println("Calculating W's")
    lim = 7
    W_dict = W_precalc(lim)
    save_object("Data/W_precalc.jld2", W_dict)
end
