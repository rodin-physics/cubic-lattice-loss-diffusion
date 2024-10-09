include("../src/main.jl")
include("../src/plotting.jl")

## PARAMETERS
M = [Inf, Inf, 0.7]         # Mass of the particle in meV * (ps / Å)²
a = 3                       # Lattice constant in Å
k1 = 520                    # Nearest neighbor force constant in meV / Å²
k2 = 170                    # Next-nearest neighbor force constant in meV / Å²
m = 3.5                     # Lattice mass in meV * (ps / Å)²
ħ = 0.6582119569            # Planck constant in meV * ps
δt = 5e-3                   # Time step in ps

# System size
size_x = size_y = size_z = 50

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
sys = system(size_x, size_y, size_z, couplings)

# Simulation parameters
speeds = [2, 5]
λs = [1 / 2, 1]
Us = [2000]

params = [(s, λ, U0) for s in speeds, λ in λs, U0 in Us] |> vec

# for p in params
Threads.@threads for p in params
    init_speed = p[1]
    λ = p[2]
    U0 = p[3]
    if !isfile(
        "Data/Reflection/Full_Reflection_Yukawa_U$(U0)_λ$(λ)_init_speed$(init_speed).jld2",
    )
        @inline function U(r)
            res = U0 * exp(-norm(r) / λ) / norm(r)
            return res
        end

        # The particle interacts with a fixed collection of atoms
        atoms = [(1, 1, n) for n = 1:1]

        # Initial conditions
        init_pos = -12 * λ
        R = [0, 0, init_pos]
        R_dot = [0, 0, init_speed]

        # Simulation parameters
        t_max = 3 * abs(init_pos / init_speed) / 2  # Maximum time
        nPts = floor(t_max / δt) |> Int
        pos = Vector{Vector{Float64}}(undef, nPts)
        speed = Vector{Vector{Float64}}(undef, nPts)

        pos_atom = Vector{Vector{Float64}}(undef, nPts)
        speed_atom = Vector{Vector{Float64}}(undef, nPts)

        pos[1] = R
        speed[1] = R_dot

        r_init = zeros(length(sys[1]))
        r_dot_init = zeros(length(sys[1]))

        pos_atom[1] = r_init[1:3]
        speed_atom[1] = r_dot_init[1:3]

        current_state = (r_init, r_dot_init, R, R_dot)

        @showprogress for ii = 2:nPts
            current_state = RKstep(m, M, a, sys, atoms, U, current_state, δt)
            pos_atom[ii] = current_state[1][1:3]
            speed_atom[ii] = current_state[2][1:3]
            pos[ii] = current_state[3]
            speed[ii] = current_state[4]
        end

        save_object(
            "Data/Reflection/Full_Reflection_Yukawa_U$(U0)_λ$(λ)_init_speed$(init_speed).jld2",
            (pos_atom, speed_atom, pos, speed),
        )
    end

end

## FIGURES
W_mat = W(m, dyn_mat, 0, (0, 0, 0))[1]

T = diag(W_mat) |> mean
L = diag(loss_mat) |> mean

data = [
    "Data/Reflection/Full_Reflection_Yukawa_U2000_λ0.5_init_speed5.jld2",
    "Data/Reflection/Full_Reflection_Yukawa_U2000_λ1.0_init_speed5.jld2",
    "Data/Reflection/Full_Reflection_Yukawa_U2000_λ0.5_init_speed2.jld2",
    "Data/Reflection/Full_Reflection_Yukawa_U2000_λ1.0_init_speed2.jld2",
]

λs = [1 / 2, 1, 1 / 2, 1]

set_theme!(CF_theme)
colors = [CF_vermillion, CF_orange, CF_green, CF_sky]
fig = Figure(size = (2400, 1000))

supertitle = fig[1, 1]
Label(
    supertitle,
    "Yukawa Interaction",
    tellwidth = false,
    tellheight = false,
    font = :latex,
    fontsize = 40,
    valign = :center,
)

main_grid = fig[2:30, 1] = GridLayout()
legend_grid = fig[31, 1] = GridLayout()

deflection_grid = main_grid[1, 1] = GridLayout()
force_grid = main_grid[2, 1] = GridLayout()

titles = [
    L"\lambda = 1/2\mathrm{\AA},\,\dot{R} = 5\mathrm{\AA/ps}",
    L"\lambda = 1\mathrm{\AA},\,\dot{R} = 5\mathrm{\AA/ps}",
    L"\lambda = 1/2\mathrm{\AA},\,\dot{R} = 2\mathrm{\AA/ps}",
    L"\lambda = 1\mathrm{\AA},\,\dot{R} = 2\mathrm{\AA/ps}",
]

axs_deflection = [
    Axis(
        deflection_grid[1, n],
        title = titles[n],
        ylabel = "Deflection (Å)",
        xlabel = "Time (ps)",
    ) for n = 1:4
]

axs_force =
    [Axis(force_grid[1, n], ylabel = "Force (meV/Å)", xlabel = "Time (ps)") for n = 1:4]

for ii = 1:4

    pos_atom, speed_atom, pos, speed = load_object(data[ii])

    pos_atom = [x[3] for x in pos_atom]
    speed_atom = [x[3] for x in speed_atom]
    pos = [x[3] for x in pos]
    speed = [x[3] for x in speed]

    λ = λs[ii]
    U0 = 2000
    @inline function U(r)
        res = U0 * exp(-norm(r) / λ) / norm(r)
        return res
    end

    U_pr = [ForwardDiff.derivative(U, x) for x in pos .- pos_atom]
    U_d_pr = [
        ForwardDiff.derivative(x -> ForwardDiff.derivative(U, x), x) for
        x in pos .- pos_atom
    ]

    times = δt .* (0:length(pos_atom)-1)
    # Actual trajectory
    r_eff = pos_atom
    lines!(axs_deflection[ii], times, pos_atom, linewidth = 4, color = CF_sky)
    lines!(
        axs_force[ii],
        times,
        [-ForwardDiff.derivative(U, x) for x in (pos .- r_eff)],
        linewidth = 4,
        color = CF_sky,
    )

    # Time-local trajectory
    r_eff = T .* U_pr - L .* U_d_pr .* (speed .- 1 .* speed_atom)
    lines!(axs_deflection[ii], times, r_eff, linewidth = 4, color = CF_vermillion)
    lines!(
        axs_force[ii],
        times,
        [-ForwardDiff.derivative(U, x) for x in (pos .- r_eff)],
        linewidth = 4,
        color = CF_vermillion,
    )

    # Time-local trajectory with framework speed neglected
    r_eff = T .* U_pr - L .* U_d_pr .* (speed .- 0 .* speed_atom)
    lines!(
        axs_deflection[ii],
        times,
        r_eff,
        linewidth = 4,
        color = CF_green,
        linestyle = :dash,
    )
    lines!(
        axs_force[ii],
        times,
        [-ForwardDiff.derivative(U, x) for x in (pos .- r_eff)],
        linewidth = 4,
        color = CF_green,
        linestyle = :dash,
    )

    # Set r = rH = 0; ṙ = ṙH = 0
    U_pr = [ForwardDiff.derivative(U, x) for x in pos]
    U_d_pr = [ForwardDiff.derivative(x -> ForwardDiff.derivative(U, x), x) for x in pos]

    r_eff = T .* U_pr - L .* U_d_pr .* speed
    lines!(
        axs_deflection[ii],
        times,
        r_eff,
        linewidth = 6,
        color = CF_blue,
        linestyle = :dot,
    )

    lines!(
        axs_force[ii],
        times,
        [-ForwardDiff.derivative(U, x) for x in (pos .- r_eff)],
        linewidth = 6,
        color = CF_blue,
        linestyle = :dot,
    )

    hlines!(axs_deflection[ii], [0.001], linewidth = 2, color = CF_red)
    hlines!(axs_force[ii], [-1], linewidth = 2, color = CF_red)

    xlims!(axs_deflection[ii], (0, times[end]))
    xlims!(axs_force[ii], (0, times[end]))
end

labs_deflection = ["(a)", "(b)", "(c)", "(d)"]
labs_force = ["(e)", "(f)", "(g)", "(h)"]

for ii = 1:4
    text!(
        axs_deflection[ii],
        0.95,
        0.95,
        text = labs_deflection[ii],
        align = (:right, :top),
        space = :relative,
        fontsize = 36,
        font = :latex,
        color = :black,
    )

    text!(
        axs_force[ii],
        0.95,
        0.95,
        text = labs_force[ii],
        align = (:right, :top),
        space = :relative,
        fontsize = 36,
        font = :latex,
        color = :black,
    )
end

hidexdecorations!.(axs_deflection)
hideydecorations!.(axs_deflection[2:end])
hideydecorations!.(axs_force[2:end])
hideydecorations!.(axs_deflection[1], label = false)
hideydecorations!.(axs_force[1], label = false)

legs = [
    LineElement(color = CF_sky, linewidth = 4),
    LineElement(color = CF_vermillion, linewidth = 4),
    LineElement(color = CF_green, linestyle = :dash, linewidth = 4),
    LineElement(color = CF_blue, linestyle = :dot, linewidth = 6),
]

Legend(
    legend_grid[1, 1],
    [legs],
    [[
        "Full solution",
        "Time local solution",
        L"Time local, $\dot{r} = 0$",
        L"Time local,  $r = \dot{r} = 0$",
    ]],
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
save("Yukawa_Reflection.pdf", fig)
