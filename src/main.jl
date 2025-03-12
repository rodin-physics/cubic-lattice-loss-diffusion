using Distributions
using ForwardDiff
using HCubature
using JLD2
using LinearAlgebra
using ProgressMeter
using QuadGK
using SparseArrays
using SpecialFunctions

include("full_approach.jl")
include("time_local_approach.jl")

dirs = [
    "Data",
    "Data/System",
    "Data/Loss",
    "Data/Diffusion",
    "Data/Diffusion_Full",
    "Data/Reflection",
]
[isdir(d) ? nothing : mkdir(d) for d in dirs]

function relaxation(sys, R, λ, U0, nStep)
    m = 1
    M = 1
    atoms = Iterators.product([1:size_x, 1:size_y, 1:size_z]...) |> collect |> vec
    atom_indices = [
        [
            get(sys[1], c, ErrorException("Coordinate not found")) for
            c in [(a..., d) for d in [1, 2, 3]]
        ] for a in atoms
    ]
    @inline function U(r)
        res = U0 * exp(-norm(r) / λ) / norm(r)
        return res
    end

    r = zeros(length(sys[1]))
    r_dot = zeros(length(sys[1]))
    if nStep > 0

        @showprogress for _ = 1:nStep
            d = derivatives(m, M, a, sys..., atoms, U, r, r_dot, R, [0, 0, 0])
            step = min(1e-1, maximum(abs.(d[2])))
            r += step .* normalize(d[2])
        end

        @showprogress for _ = 1:nStep
            d = derivatives(m, M, a, sys..., atoms, U, r, r_dot, R, [0, 0, 0])
            step = min(1e-2, maximum(abs.(d[2])))
            r += step .* normalize(d[2])
        end

        @showprogress for _ = 1:nStep
            d = derivatives(m, M, a, sys..., atoms, U, r, r_dot, R, [0, 0, 0])
            step = min(1e-3, maximum(abs.(d[2])))
            r += step .* normalize(d[2])
        end

        @showprogress for _ = 1:nStep
            d = derivatives(m, M, a, sys..., atoms, U, r, r_dot, R, [0, 0, 0])
            step = min(1e-4, maximum(abs.(d[2])))
            r += step .* normalize(d[2])
        end

        @showprogress for _ = 1:nStep
            d = derivatives(m, M, a, sys..., atoms, U, r, r_dot, R, [0, 0, 0])
            step = min(1e-5, maximum(abs.(d[2])))
            r += step .* normalize(d[2])
        end

    end
    atom_positions = [r[atom_indices[n]] .+ (atoms[n] .- 1) .* a for n in eachindex(atoms)]
    interaction = [U(pos - R) for pos in atom_positions] |> sum
    elastic = r' * sys[2] * r / 2
    max_disp = maximum(abs.(r))
    d = derivatives(m, M, a, sys..., atoms, U, r, r_dot, R, [0, 0, 0])
    return (interaction, elastic, max_disp, d[2], r)
end
