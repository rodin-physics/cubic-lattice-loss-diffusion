using Random
include("../src/main.jl")
include("../src/plotting.jl")

## PARAMETERS
M = 0.7                     # Mass of the particle in meV * (ps / Å)²
a = 3                       # Lattice constant in Å
k1 = 520                    # Nearest neighbor force constant in meV / Å²
k2 = 170                    # Next-nearest neighbor force constant in meV / Å²
m = 3.5                     # Lattice mass in meV * (ps / Å)²
ħ = 0.6582119569            # Planck constant in meV * ps

# Initial conditions
R = [a, a, a] / 2
R_dot = [0, 0, 0]
# Simulation parameters
δt = 5e-3                   # Time step in ps
t_max = 3000                # Time in ps
nPts = floor(t_max / δt) |> Int

# LATTICE
size_x = size_y = size_z = 20         # Used ONLY for harmonics

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
dyn_mat = dynamical_matrix(m, couplings)
dyn_mat_small = dynamical_matrix_small(a, m, couplings)
loss_mat = Loss_Matrix(a, m, dyn_mat_small)[1] |> real

eigs = [dyn_mat(q) |> eigen for q in momenta]
ηs = hcat([ei.vectors for ei in eigs]...)
Ωs = vcat([sqrt.(ei.values |> real) for ei in eigs]...)
qs = vcat([[collect(q), collect(q), collect(q)] for q in momenta]...)

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

if !isfile("Data/Diffusion/W_precalc.jld2")
    println("Calculating W's")
    lim = 2 * extent + 1
    W_dict = W_precalc(lim)
    save_object("Data/Diffusion/W_precalc.jld2", W_dict)
end

W_dict = load_object("Data/W_precalc.jld2")

# CALCULATIONS
ħΩTs = [25, 50]
U0s = [1000, 2000, 4000]

λ = a / 4

extent = 1                  # Interaction length
params = [(ħΩT, U0) for ħΩT in ħΩTs, U0 in U0s] |> vec

Threads.@threads for par in params
    Random.seed!(150)
    ħΩT = par[1]
    U0 = par[2]

    @inline function U(r)
        res = U0 * exp(-dot(r, r) / 2 / λ^2)
        return res
    end

    ΩT = ħΩT / ħ

    ζs = [Amplitude(ħ, Ωq, ΩT) * exp(-2im * π * rand()) for Ωq in Ωs]
    ζsηs =
        [ζs[j] .* ηs[:, j] for j in eachindex(ζs)] ./ √(size_x * size_y * size_z) ./ sqrt(m)
    ζsηs_dot = -1im .* Ωs .* ζsηs

    # Statistics check: get a trajectory of one of the masses for an extended period of time
    if !isfile("Data/Diffusion/Homogeneous_SingleMass_ħΩT$(ħΩT)_size$(size_x).jld2")
        # t_extended = 3000
        δt_extended = 5e-2
        ts = range(0, t_max, step = δt_extended)

        Ds = Vector{Vector{Float64}}(undef, length(ts))
        Vs = Vector{Vector{Float64}}(undef, length(ts))
        pr = Progress(length(ts))
        for ii in eachindex(ts)
            # Threads.@threads for ii in eachindex(ts)

            time_phase = exp.(-1im * ts[ii] .* Ωs)
            harmonics_pos = ζsηs .* time_phase
            harmonics_vel = ζsηs_dot .* time_phase

            Ds[ii] = real(sum(harmonics_pos))
            Vs[ii] = real(sum(harmonics_vel))
            next!(pr)
        end
        disp = [norm(d) for d in Ds]
        speed = [norm(v) for v in Vs]
        save_object(
            "Data/Diffusion/Homogeneous_SingleMass_ħΩT$(ħΩT)_size$(size_x).jld2",
            (disp, speed),
        )
    end

    pos_particle = Vector{Vector{Float64}}(undef, nPts)
    speed_particle = Vector{Vector{Float64}}(undef, nPts)
    pos_particle[1] = R
    speed_particle[1] = R_dot

    current_state = (R, R_dot)

    if !isfile("Data/Diffusion/Diffusion_U0$(U0)_λ$(λ)_ħΩT$(ħΩT)_size$(size_x).jld2")
        @showprogress for ii = 2:nPts
            t = (ii - 1) * δt

            # Get the lattice atoms with which the particle interacts
            # RECALL (1,1,1) ATOM is at (0,0,0) COORDINATE
            pos = ceil.(Int, (current_state[1] ./ a))

            atoms =
                Iterators.product(
                    [
                        pos[1]+1-extent:pos[1]+extent,
                        pos[2]+1-extent:pos[2]+extent,
                        pos[3]+1-extent:pos[3]+extent,
                    ]...,
                ) |>
                collect |>
                vec

            current_state = RKstep_local(
                current_state,
                ζsηs,
                ζsηs_dot,
                qs,
                Ωs,
                atoms,
                U,
                t,
                δt,
                loss_mat,
            )
            pos_particle[ii] = current_state[1]
            speed_particle[ii] = current_state[2]
            GC.safepoint()
        end
        save_object(
            "Data/Diffusion/Diffusion_U0$(U0)_λ$(λ)_ħΩT$(ħΩT)_size$(size_x).jld2",
            (pos_particle, speed_particle),
        )
    end

end






Random.seed!(150)
ħΩT = 25
U0 = 4000

@inline function U(r)
    res = U0 * exp(-dot(r, r) / 2 / λ^2)
    return res
end

ΩT = ħΩT / ħ

ζs = [Amplitude(ħ, Ωq, ΩT) * exp(-2im * π * rand()) for Ωq in Ωs]
ζsηs = [ζs[j] .* ηs[:, j] for j in eachindex(ζs)] ./ √(size_x * size_y * size_z) ./ sqrt(m)
ζsηs_dot = -1im .* Ωs .* ζsηs

# Statistics check: get a trajectory of one of the masses for an extended period of time


pos_particle = Vector{Vector{Float64}}(undef, nPts)
speed_particle = Vector{Vector{Float64}}(undef, nPts)
pos_particle[1] = R
speed_particle[1] = R_dot

current_state = (R, R_dot)
ii = 2
t = (ii - 1) * δt

# Get the lattice atoms with which the particle interacts
# RECALL (1,1,1) ATOM is at (0,0,0) COORDINATE
pos = ceil.(Int, (current_state[1] ./ a))
extent = 1
atoms =
    Iterators.product(
        [
            pos[1]+1-extent:pos[1]+extent,
            pos[2]+1-extent:pos[2]+extent,
            pos[3]+1-extent:pos[3]+extent,
        ]...,
    ) |>
    collect |>
    vec



# Times used in RK routine for the time-dependent portion
ts = [t, t + δt / 4, t + δt / 4, t + δt / 2, t + 3 * δt / 4, t + δt]

# Compute rH for the times
rHs = Vector{Tuple{Vector{Vector{Float64}},Vector{Vector{Float64}}}}(undef, 6)

if isnothing(ζsηs) | isnothing(ζsηs_dot)
    for ii in eachindex(ts)
        res_pos = [atom .- [1, 1, 1] for atom in atoms] .* a |> real
        res_vel = [atom .- [1, 1, 1] for atom in atoms] .* 0 |> real

        rHs[ii] = (res_pos, res_vel)
    end
else
    phases = [exp(1im * dot(q, atom)) for atom in atoms, q in qs]
    for ii in eachindex(ts)
        # Threads.@threads for ii in eachindex(ts)
        time_phase = exp.(-1im * ts[ii] .* Ωs)
        harmonics_pos = ζsηs .* time_phase
        harmonics_vel = ζsηs_dot .* time_phase
        res_pos =
            eachrow(phases * transpose(hcat(harmonics_pos...))) .+
            [atom .- [1, 1, 1] for atom in atoms] .* a |> real
        res_vel =
            eachrow(phases * transpose(hcat(harmonics_vel...))) .+
            [atom .- [1, 1, 1] for atom in atoms] .* 0 |> real

        rHs[ii] = (res_pos, res_vel)
    end
end



# rHs[1]


derivatives_local(atoms, rHs[1]..., s1..., U, loss_mat)
rHs[1][1]

@benchmark r_eff, R_eff = pos_eff(atoms, rHs[1][1], rHs[1][2], R, R_dot, U, loss_mat)
r_eff
# rHs[1][1]
# r_eff
# current_state =
#     RKstep_local(current_state, ζsηs, ζsηs_dot, qs, Ωs, atoms, U, t, δt, loss_mat)
# pos_particle[ii] = current_state[1]
# speed_particle[ii] = current_state[2]
# GC.safepoint()
r = rHs[1][1]
r_dot = rHs[1][2]
loss_mat * (
    ForwardDiff.gradient(
        R_ ->
            sum([dot(ForwardDiff.gradient(U, r[n] .- R_), r_dot[n]) for n in eachindex(r)]),
        R,
    ) + ForwardDiff.hessian(R -> sum([U(R - x) for x in r]), R) * R_dot
) |> real
ForwardDiff.gradient(
    R_ -> sum([dot(ForwardDiff.gradient(U, r[n] .- R_), r_dot[n]) for n in eachindex(r)]),
    R,
)
# W_dd(m, dyn_mat, 0.1, [0,0,0])[1]

# ts = range(0, .5, length = 20)
# rr = [W_dd(m, dyn_mat, t, [0,0,0])[1][1,1] for t in ts]
# ww = [W(m, dyn_mat, t, [0,0,0])[1][1,1] for t in ts]
# W(m, dyn_mat, 0.2, [0,0,0])[1]
# lines(ts, rr)
# lines(ts, ww)


# W_dict
