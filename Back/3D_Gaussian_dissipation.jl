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
size_x = size_y = size_z = 100

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
U0 = 50
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

## Full solution
sys = system(size_x, size_y, size_z, couplings)

pos = Vector{Vector{Float64}}(undef, nPts)
speed = Vector{Vector{Float64}}(undef, nPts)
pos[1] = R
speed[1] = R_dot
# No homogeneous motion
r_init = zeros(length(sys[1]))
r_dot_init = zeros(length(sys[1]))

function solve_full()
    current_state = (r_init, r_dot_init, R, R_dot)

    @showprogress for ii = 2:nPts
        current_state = RKstep(m, M, a, sys, atoms, U, current_state, δt)
        pos[ii] = current_state[3]
        speed[ii] = current_state[4]
    end
    return (pos, speed)
end

## Single-particle solution
# Get the recoil terms
include("precalculation.jl")
W_dict = load_object("Data/W_precalc.jld2")

pos_particle = Vector{Vector{Float64}}(undef, nPts)
speed_particle = Vector{Vector{Float64}}(undef, nPts)
pos_particle[1] = R
speed_particle[1] = R_dot

function solve_particle()
    current_state = (R, R_dot)
    @showprogress for ii = 2:nPts
        current_state = RKstep_local(
            current_state,
            nothing,
            nothing,
            qs,
            Ωs,
            atoms,
            U,
            δt * ii,
            δt,
            loss_mat,
        )
        pos_particle[ii] = current_state[1]
        speed_particle[ii] = current_state[2]
    end
    return (pos_particle, speed_particle)

end

if !isfile("Data/Loss/Local_Loss_U$(U0)_λ$(λ)_TEST2.jld2")
    pos_particle, speed_particle = solve_particle()
    save_object(
        "Data/Loss/Local_Loss_U$(U0)_λ$(λ)_TEST2.jld2",
        (pos_particle, speed_particle),
    )
end

if !isfile("Data/Loss/Full_Loss_U$(U0)_λ$(λ).jld2")
    pos, speed = solve_full()
    save_object("Data/Loss/Full_Loss_U$(U0)_λ$(λ).jld2", (pos, speed))
end

## PLOTTING

pos_local, speed_local = load_object("Data/Loss/Local_Loss_U$(U0)_λ$(λ)_TEST2.jld2")
pos_full, speed_full = load_object("Data/Loss/Full_Loss_U$(U0)_λ$(λ).jld2")

## FIGURES
set_theme!(CF_theme)
colors = [CF_vermillion, CF_orange, CF_green, CF_sky]
fig = Figure(size = (1200, 1000))

supertitle = fig[1, 1]
Label(
    supertitle,
    "Trapped dissipation",
    tellwidth = false,
    tellheight = false,
    font = :latex,
    fontsize = 36,
    valign = :center,
)

main_grid = fig[2:30, 1] = GridLayout()
legend_grid = fig[31, 1] = GridLayout()

ax_pos = Axis(
    main_grid[1, 1],
    # title = "Position",
    ylabel = L"Position $R$ (Å)",
    xlabel = L"Time $t$ (ps)",
)

ax_speed = Axis(
    main_grid[2, 1],
    # title = "Velocity",
    ylabel = L"Velocity $\dot{R}$ (Å/ps)",
    xlabel = L"Time $t$ (ps)",
)
ax = [ax_pos, ax_speed]
labs = ["(a)", "(b)"]

for ii = 1:2
    text!(
        ax[ii],
        0.95,
        0.95,
        text = labs[ii],
        align = (:right, :top),
        space = :relative,
        fontsize = 36,
        font = :latex,
        color = :black,
    )
end

lines!(
    ax_pos,
    δt .* (1:length(pos_full)),
    [x[3] for x in pos_full],
    color = CF_sky,
    linewidth = 4,
)
lines!(
    ax_speed,
    δt .* (1:length(speed_full)),
    [x[3] for x in speed_full],
    color = CF_sky,
    linewidth = 4,
)

lines!(
    ax_pos,
    δt .* (1:length(pos_local)),
    [x[3] for x in pos_local],
    color = CF_vermillion,
    linewidth = 4,
)
lines!(
    ax_speed,
    δt .* (1:length(speed_local)),
    [x[3] for x in speed_local],
    color = CF_vermillion,
    linewidth = 4,
)
hidexdecorations!(ax_pos)
polys =
    [PolyElement(color = c, strokecolor = :transparent) for c in [CF_sky, CF_vermillion]]

Legend(
    legend_grid[1, 1],
    [polys],
    [["Full solution", "Time-local solution"]],
    [""],
    halign = :center,
    valign = :center,
    tellheight = false,
    tellwidth = false,
    framevisible = false,
    orientation = :horizontal,
    titlevisible = false,
    titleposition = :left,
    titlefont = :latex,
)

fig
save("Trapped_Dissipation.pdf", fig)
