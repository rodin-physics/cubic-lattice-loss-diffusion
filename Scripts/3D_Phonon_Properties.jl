using Random
include("../src/main.jl")
include("../src/plotting.jl")

Random.seed!(150)
## PARAMETERS
a = 3                       # Lattice constant in Å
k1 = 520                    # Nearest neighbor force constant in meV / Å²
k2 = 170                    # Next-nearest neighbor force constant in meV / Å²
m = 3.5                     # Lattice mass in meV * (ps / Å)²
ħ = 0.6582119569            # Planck constant in meV * ps

# LATTICE
couplings = cubic_lattice(k1, k2)
# Dynamical matrix
DynamicalMatrix = dynamical_matrix(m, couplings)

# Precalculate W

if !isfile("Data/W_precalc.jld2")
    println("Calculating W's")
    lim = 7
    displacements = Iterators.product([-lim:lim, -lim:lim, -lim:lim]...) |> collect |> vec
    Ws = Vector{Matrix{Float64}}(undef, length(displacements))
    pr = Progress(length(Ws))
    Threads.@threads for ii in eachindex(Ws)
        Ws[ii] = W(m, DynamicalMatrix, 0, displacements[ii])[1]
        next!(pr)
    end
    W_dict = Dict(zip(displacements, Ws))
    save_object("Data/W_precalc.jld2", W_dict)
end

# BAND STRUCTURE
Γ = [0, 0, 0]
X = [π, 0, 0]
M = [π, π, 0]
R = [π, π, π]

nPts = 200

ΓX = [Γ + (X - Γ) ./ nPts * n for n = 0:nPts]
XM = [X + (M - X) ./ nPts * n for n = 1:nPts]
MΓ = [M + (Γ - M) ./ nPts * n for n = 1:nPts]
ΓR = [Γ + (R - Γ) ./ nPts * n for n = 1:nPts]

path = vcat([ΓX, XM, MΓ, ΓR]...)
energies = [sqrt.(real.(eigen(DynamicalMatrix(v)).values)) for v in path]

plot_positions = zeros(length(path))

for ii in eachindex(path)[2:end]
    plot_positions[ii] = plot_positions[ii-1] + norm(path[ii] - path[ii-1])
end

## RECOIL
if !isfile("Data/System/Recoil.jld2")
    nPts = 200
    res = zeros(nPts)
    tmin = 0
    tmax = 2
    ts = range(tmin, tmax, length = nPts)
    pr = Progress(nPts)
    Threads.@threads for ii in eachindex(res)
        res[ii] = W(m, DynamicalMatrix, ts[ii], [0, 0, 0])[1][1, 1]
        next!(pr)
        GC.safepoint()
    end
    save_object("Data/System/Recoil.jld2", (ts, res))
end

## HOMOGENEOUS MOTION
ħΩTs = [1e-3, 5, 25, 50]
size_x = size_y = size_z = 60
Threads.@threads for ħΩT in ħΩTs
    if !isfile("Data/System/Homogeneous_ħT$(ħΩT).jld2")
        disp, speed =
            homogeneous_init(ħ, m, DynamicalMatrix, ħΩT / ħ, size_x, size_y, size_z)

        save_object("Data/System/Homogeneous_ħT$(ħΩT).jld2", (disp, speed))
    end
end

## FIGURES
set_theme!(CF_theme)
fig = Figure(size = (1200, 800))

ax_bands = Axis(fig[1, 1], title = "Phonon dispersion", ylabel = "Energy (meV)")

ax_recoil = Axis(
    fig[1, 2],
    title = "Recoil kernel",
    ylabel = L"W_{11}(t) / W_{11}(0)",
    xlabel = L"$t$ (ps)",
)

ax_pos = Axis(
    fig[2, 1],
    title = "Displacement distribution",
    xlabel = L"r \, (\AA)",
    ylabel = "arb. units.",
)

ax_speed =
    Axis(fig[2, 2], title = "Speed distribution", xlabel = L"\dot{r}\,(\AA/\mathrm{ps})")

## BAND PLOTTING
for ii = 1:3
    lines!(
        ax_bands,
        plot_positions,
        [e[ii] for e in energies] .* ħ,
        linewidth = 4,
        color = CF_vermillion,
    )
end

x_ticks = [
    plot_positions[1],
    plot_positions[nPts+1],
    plot_positions[2*nPts+1],
    plot_positions[3*nPts+1],
    plot_positions[4*nPts+1],
]
x_labels = [L"Γ", "X", "M", L"\Gamma", "R"]

ax_bands.xticks = (x_ticks, x_labels)
xlims!(ax_bands, (plot_positions[1], plot_positions[end]))
ylims!(ax_bands, (0, 22))
xlims!(ax_recoil, (0, 1.99))
## RECOIL PLOTTING
(ts, res) = load_object("Data/System/Recoil.jld2")
lines!(ax_recoil, ts, res ./ res[1], linewidth = 4, color = CF_sky)

## THERMAL STATISTICS PLOTTING

data = [load_object("Data/System/Homogeneous_ħT$(ħΩT).jld2") for ħΩT in ħΩTs]
σ_pos = [Corr_pos(ħ, m, DynamicalMatrix, ħΩT / ħ, 0, [0, 0, 0])[1][1, 1] for ħΩT in ħΩTs]
σ_speed =
    [Corr_speed(ħ, m, DynamicalMatrix, ħΩT / ħ, 0, [0, 0, 0])[1][1, 1] for ħΩT in ħΩTs]
ΩT_label = [L"10^{-3}", L"5", L"25", L"50"]
colors = [CF_red, CF_yellow, CF_green, CF_blue]

for jj in eachindex(data)
    pos_data = data[jj][1]
    speed_data = data[jj][2]

    pos_sq = pos_data .^ 2
    speed_sq = speed_data .^ 2

    N = size_x * size_y * size_z

    r = [sum(pos_sq[(3*ii+1):(3*ii+3)]) for ii = 0:(N-1)]
    v = [sum(speed_sq[(3*ii+1):(3*ii+3)]) for ii = 0:(N-1)]

    hist!(
        ax_pos,
        sqrt.(r),
        bins = 100,
        normalization = :pdf,
        color = alphacolor(colors[jj], 0.5),
    )
    hist!(
        ax_speed,
        sqrt.(v),
        bins = 100,
        normalization = :pdf,
        color = alphacolor(colors[jj], 0.5),
    )

    plotPts = 500
    pos = range(0, 1, length = plotPts)
    speed = range(0, 15, length = plotPts)

    lines!(
        ax_pos,
        pos,
        4 * π .* pos .^ 2 ./ √(2 * π * σ_pos[jj])^3 .* exp.(-pos .^ 2 ./ 2 ./ σ_pos[jj]),
        color = colors[jj],
        linewidth = 4,
    )

    lines!(
        ax_speed,
        speed,
        4 * π .* speed .^ 2 ./ √(2 * π * σ_speed[jj])^3 .*
        exp.(-speed .^ 2 ./ 2 ./ σ_speed[jj]),
        color = colors[jj],
        linewidth = 4,
        label = ΩT_label[jj],
    )
end
xlims!(ax_pos, (0, 0.749))
xlims!(ax_speed, (0, 14.9))
ylims!(ax_pos, (0, 9))
ylims!(ax_speed, (0, 1 / 2))
hideydecorations!(ax_pos, label = false)
hideydecorations!(ax_speed)
axislegend(L"\hbar\Omega_T (\mathrm{meV})", orientation = :horizontal, framevisible = false)

## PANEL LABELS
labs = ["(b)", "(c)", "(d)", "(e)"]
ax = [ax_bands, ax_recoil, ax_pos, ax_speed]
for jj in eachindex(ax)
    text!(
        ax[jj],
        0.11,
        0.975,
        text = labs[jj],
        align = (:right, :top),
        space = :relative,
        fontsize = 36,
        font = :latex,
        color = :black,
    )

end
fig

save("Framework_properties.pdf", fig)
