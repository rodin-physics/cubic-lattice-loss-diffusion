include("../src/main.jl")
include("../src/plotting.jl")

## PARAMETERS
M = [Inf, Inf, 0.7]         # Mass of the particle in meV * (ps / Å)²
a = 3                       # Lattice constant in Å
k1 = 520                    # Nearest neighbor force constant in meV / Å²
k2 = 170                    # Next-nearest neighbor force constant in meV / Å²
m = 3.5                     # Lattice mass in meV * (ps / Å)²
ħ = 0.6582119569            # Planck constant in meV * ps

# System size
size_x = size_y = size_z = 20

# Momentum grid
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

# Create the lattice
couplings = cubic_lattice(k1, k2)
dyn_mat = dynamical_matrix(m, couplings)
dyn_mat_small = dynamical_matrix_small(a, m, couplings)
loss_mat = Loss_Matrix(a, m, dyn_mat_small)[1] |> real

# Solve the eigenproblem for all the momenta on the grid
eigs = [dyn_mat(q) |> eigen for q in momenta]
ηs = hcat([ei.vectors for ei in eigs]...)
Ωs = vcat([sqrt.(ei.values |> real) for ei in eigs]...)
qs = vcat([[collect(q), collect(q), collect(q)] for q in momenta]...)

# Interaction properties
U0 = 500
λ = 3 / 4
@inline function U(r)
    res = U0 * exp(-dot(r, r) / 2 / λ^2)
    return res
end
# The particle interacts with a fixed collection of atoms
atoms = [(1, 1, n) for n = 1:8]

# Simulation parameters
δt = 5e-3                   # Time step in ps
t_max = 30                  # Maximum time
nPts = floor(t_max / δt) |> Int

# Initial conditions
init_pos = 3.5 * a
init_speed = 7.5
R = [0, 0, init_pos]
R_dot = [0, 0, init_speed]

# ## Full solution
# sys = system(size_x, size_y, size_z, couplings)

# pos = Vector{Vector{Float64}}(undef, nPts)
# speed = Vector{Vector{Float64}}(undef, nPts)
# pos[1] = R
# speed[1] = R_dot
# # No homogeneous motion
# r_init = zeros(length(sys[1]))
# r_dot_init = zeros(length(sys[1]))

# function solve_full()
#     current_state = (r_init, r_dot_init, R, R_dot)

#     @showprogress for ii = 2:nPts
#         current_state = RKstep(m, M, a, sys, atoms, U, current_state, δt)
#         pos[ii] = current_state[3]
#         speed[ii] = current_state[4]
#     end
#     return (pos, speed)
# end

## Single-particle solution
# Get the recoil terms
include("precalculation.jl")
W_dict = load_object("Data/W_precalc.jld2")

pos_particle = Vector{Vector{Float64}}(undef, nPts)
speed_particle = Vector{Vector{Float64}}(undef, nPts)
pos_particle[1] = R
speed_particle[1] = R_dot

rH = [atom .- [1, 1, 1] for atom in atoms] .* a |> real
rH_dot = [atom .- [1, 1, 1] for atom in atoms] .* 0 |> real

for ii = 1:20
    rH, R = pos_eff(atoms, rH, rH_dot, R, R_dot, U, loss_mat)
end
R
function pos_eff(atoms, r, r_dot, R, R_dot, U, loss_mat)

    # Get ∇ᵣ assuming a pairwise interaction for a 
    # set of framework positions r and the particle position R
    gradient_r = [ForwardDiff.gradient(U, x .- R) for x in r]
    # The correction to nth position is given by ∑ₖ Wₙ₋ₖ grad_r[k]
    r_corr = [
        sum([
            get(W_dict, atoms[n] .- atoms[k], nothing) * gradient_r[k] for
            k in eachindex(atoms)
        ]) for n in eachindex(atoms)
    ]

    # The correction to R is given by loss_mat * ∇_R dU/dt = loss_mat * (∇_R (∇_r U ⋅ ṙ) + Hess(U)⋅Ṙ)
    R_corr =
        loss_mat * (
            ForwardDiff.gradient(
                R_ -> sum([
                    dot(ForwardDiff.gradient(U, r[n] .- R_), r_dot[n]) for n in eachindex(r)
                ]),
                R,
            ) + ForwardDiff.hessian(R -> sum([U(R - x) for x in r]), R) * R_dot
        ) |> real

    r_eff = r .- r_corr
    R_eff = R + R_corr

    return (r_eff, R_eff)
end

# rHs[ii] = (res_pos, res_vel)

function init_state(atoms, r, r_dot, R, R_dot, U, loss_mat)

    # Get ∇ᵣ assuming a pairwise interaction for a 
    # set of framework positions r and the particle position R
    gradient_r = [ForwardDiff.gradient(U, x .- R) for x in r]
    # The correction to nth position is given by ∑ₖ Wₙ₋ₖ grad_r[k]

    δ =
        -[
            sum([
                get(W_dict, atoms[n] .- atoms[k], nothing) * gradient_r[k] for
                k in eachindex(atoms)
            ]) for n in eachindex(atoms)
        ]


    Hess_r = [ForwardDiff.hessian(U, x .- R) for x in r]
    Hess_r .* δ

    δ1 =
        -[
            sum([
                get(W_dict, atoms[n] .- atoms[k], nothing) * (Hess_r.*δ)[k] for
                k in eachindex(atoms)
            ]) for n in eachindex(atoms)
        ]
    # Hess_R = sum(Hess_r)
    # Grad_r_R = -hcat(Hess_r...)
    # Grad_R_r = Grad_r_R'

    # Hess_r = Matrix(blockdiag(sparse.(Hess_r)...))

    # V_inv = [
    #     get(W_dict, atoms[n] .- atoms[k], nothing) for k in eachindex(atoms),
    #     n in eachindex(atoms)
    # ]

    # V_inv = hcat([vcat(row...) for row in eachrow(V_inv)]...)

    # V_inv = (V_inv + V_inv') / 2

    # Δ1 = loss_mat * (Hess_R * R_dot)

    # Δ2 =
    #     loss_mat *
    #     Grad_r_R *
    #     inv((I(length(atoms) * 3) + V_inv * Hess_r)) *
    #     (-V_inv * Grad_R_r * R_dot)

    # # println(Δ1)
    # # println(Δ2)

    # # return (r .+ δ, R + Δ1 + Δ2)
    # return (r .+ δ)
    # return (R + Δ1 + Δ2)
    return (δ, δ1)
end

for ii = 1:20
    rH, R = init_state(atoms, rH, rH_dot, R, R_dot, U, loss_mat)
end
R
R = init_state(atoms, rH, rH_dot, R, R_dot, U, loss_mat)
for ii = 1:100
    rH = init_state(atoms, rH, rH_dot, R, R_dot, U, loss_mat)
end

init_state(atoms, rH, rH_dot, R, R_dot, U, loss_mat)[1]
init_state(atoms, rH, rH_dot, R, R_dot, U, loss_mat)[2]
vcat(rH...)
U(1.5)

R

Hess_r = [ForwardDiff.hessian(U, x .- R) for x in rH]
Hess_R = sum(Hess_r)
R = R + loss_mat * Hess_R * R_dot

rH = rH .+ init_state(atoms, rH, rH_dot, R, R_dot, U, loss_mat)[1]


res1 = load_object("Data/Loss/Local_Loss_U50_λ0.75_TEST.jld2")
res2 = load_object("Data/Loss/Local_Loss_U50_λ0.75_TEST2.jld2")
lines([x[3] for x in (res1[1] .- res2[1])])
