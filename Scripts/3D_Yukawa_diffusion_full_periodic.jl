using Random
include("../src/main.jl")
include("../src/plotting.jl")
Random.seed!(150)
## PARAMETERS
M = 0.7                     # Mass of the particle in meV * (ps / Å)²
a = 3                       # Lattice constant in Å
k1 = 520                    # Nearest neighbor force constant in meV / Å²
k2 = 170                    # Next-nearest neighbor force constant in meV / Å²
m = 3.5                     # Lattice mass in meV * (ps / Å)²
ħ = 0.6582119569            # Planck constant in meV * ps

# Simulation parameters
δt = 5e-3                   # Time step in ps
t_max = 30                  # Time in ps
nPts = floor(t_max / δt) |> Int
extent = 1                  # Interaction length
# LATTICE
size_x = size_y = size_z = 20

momenta =
    Iterators.product(
        [
            2 * π .* (0:(size_x-1)) ./ size_x,
            2 * π .* (0:(size_y-1)) ./ size_y,
            2 * π .* (0:(size_z-1)) ./ size_z,
        ]...,
    ) |>
    collect |>
    vec
# Drop the (0,0,0) momentum
momenta = momenta[2:end]

couplings = cubic_lattice(k1, k2)
sys = system(size_x, size_y, size_z, couplings)

dyn_mat = dynamical_matrix(m, couplings)

eigs = [dyn_mat(q) |> eigen for q in momenta]
ηs = hcat([ei.vectors for ei in eigs]...)
Ωs = vcat([sqrt.(ei.values |> real) for ei in eigs]...)
qs = vcat([[collect(q), collect(q), collect(q)] for q in momenta]...)

# Initial conditions
R = [size_x + 1, size_y + 1, size_z + 1] .* a / 2
R_dot = [0, 0, 0]

# CALCULATIONS
ħΩTs = [25, 50]
U0s = [4000, 8000]

λ = 1 / 2

params = [(x, 4000) for x in [50]]
# params = [(x, 4000) for x in [10, 15, 20, 25, 30, 40, 50]]

for par in params
    ħΩT = par[1]
    U0 = par[2]

    if !isfile(
        "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U0$(U0)_λ$(λ)_ħΩT$(ħΩT)_size$(size_x).jld2",
    )

        @inline function U(r)
            res = U0 * exp(-norm(r) / λ) / norm(r)
            return res
        end

        ΩT = ħΩT / ħ

        # Calculate n_traj trajectories

        pos_particle = Vector{Vector{Float64}}(undef, nPts)
        speed_particle = Vector{Vector{Float64}}(undef, nPts)
        pos_particle[1] = R
        speed_particle[1] = R_dot

        println(ΩT)
        println(U0)
        r_init, r_dot_init = homogeneous_init(ħ, m, dyn_mat, ΩT, size_x, size_y, size_z)

        current_state = (r_init, r_dot_init, R, R_dot)

        # pos = ceil.(Int, (pos_particle[1] ./ a))
        # atoms =
        #     Iterators.product(
        #         [
        #             pos[1]+1-extent:pos[1]+extent,
        #             pos[2]+1-extent:pos[2]+extent,
        #             pos[3]+1-extent:pos[3]+extent,
        #         ]...,
        #     ) |>
        #     collect |>
        #     vec
        # atoms = [mod1.(atom, size_x) for atom in atoms]
        @showprogress for ii = 2:nPts
            # RKstep_periodic(m, M, a, sys, U, current_state, δt)
            current_state = RKstep_periodic(m, M, a, sys, U, current_state, δt)
            # current_state = RKstep_periodic(m, M, a, sys, atoms, U, current_state, δt)
            pos_particle[ii] = current_state[3]
            speed_particle[ii] = current_state[4]
            # Get the unit cell in which the particle currently is
            # pos = ceil.(Int, (pos_particle[ii] ./ a))
            # # The the framework atoms bounding the unit cell
            # atoms =
            #     Iterators.product(
            #         [
            #             pos[1]+1-extent:pos[1]+extent,
            #             pos[2]+1-extent:pos[2]+extent,
            #             pos[3]+1-extent:pos[3]+extent,
            #         ]...,
            #     ) |>
            #     collect |>
            #     vec
            # # Fold the atoms back into the finite size system
            # atoms = [mod1.(atom, size_x) for atom in atoms]
            GC.safepoint()
        end

        save_object(
            "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U0$(U0)_λ$(λ)_ħΩT$(ħΩT)_size$(size_x).jld2",
            (
                pos_particle,
                speed_particle,
                r_init,
                r_dot_init,
                current_state[1],
                current_state[2],
            ),
        )
    end

end
