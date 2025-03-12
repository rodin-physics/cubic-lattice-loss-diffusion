include("3D_Gaussian_diffusion.jl")
include("../src/plotting.jl")

data = [
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT10_size20.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT15_size20.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT20_size20.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT25_size20.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT30_size20.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT40_size20.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT50_size20.jld2",
    "Data/Diffusion/Diffusion_U02000_λ0.75_ħΩT25_size20.jld2",
    "Data/Diffusion/Diffusion_U02000_λ0.75_ħΩT50_size20.jld2",
]

## DISPLACEMENT VARIANCE CALCULATION
nSteps = floor(Int, t_max / δt / 2)
skp = 1000

nPts = floor(Int, nSteps / skp)

times = skp .* (1:nPts) .* δt

saves = [
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT10_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT15_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT20_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT25_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT30_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT40_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT50_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U02000_λ0.75_ħΩT25_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U02000_λ0.75_ħΩT50_size20_Disp_Var.jld2",
]

for ii in eachindex(data)
    if !isfile(saves[ii])
        d = load_object(data[ii])
        loc = d[1]

        res = zeros(nPts)

        @showprogress for time = 1:nPts
            disps = loc[(1+skp*time):end] - loc[1:(end-skp*time)]
            res[time] = [dot(x, x) for x in disps] |> mean
        end
        save_object(saves[ii], (times, res))
    end

end


# DISPLACEMENT VARIANCE

data = [
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT25_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT50_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U02000_λ0.75_ħΩT25_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U02000_λ0.75_ħΩT50_size20_Disp_Var.jld2",
]

prefactor = [3.5, 14, 1 / 20, 2 / 3]

colors = [CF_sky, CF_vermillion, CF_sky, CF_vermillion]
markers = [:circle, :circle, :cross, :cross]
labs = [
    L"$\hbar\Omega_T = 25$meV, $U_0 = 1$eV",
    L"$\hbar\Omega_T = 50$meV, $U_0 = 1$eV",
    L"$\hbar\Omega_T = 25$meV, $U_0 = 2$eV",
    L"$\hbar\Omega_T = 50$meV, $U_0 = 2$eV",
]

set_theme!(CF_theme)
fig = Figure(size = (1200, 800))

ax = Axis(
    fig[1, 1],
    title = "Displacement Variance",
    titlesize = 44,
    xscale = log10,
    yscale = log10,
    xlabel = "Time (ps)",
    ylabel = L"$\langle |\Delta \mathbf{R}|^2\rangle$ (\AA)",
)

for ii in eachindex(data)
    d = load_object(data[ii])

    scatter!(
        ax,
        d[1],
        d[2],
        color = colors[ii],
        marker = markers[ii],
        markersize = 14,
        label = labs[ii],
    )
    lines!(ax, d[1], prefactor[ii] .* d[1], color = colors[ii], linewidth = 2)

end
axislegend(ax, position = :lt, framevisible = false)
fig
save("Disp_Var.pdf", fig)

## DIFFUSIVITY ESTIMATION

data = [
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT10_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT15_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT20_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT25_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT30_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT40_size20_Disp_Var.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT50_size20_Disp_Var.jld2",
]
data_pts = 5
ħΩs = [10, 15, 20, 25, 30, 40, 50]
ln_D = zeros(length(ħΩs))

for ii in eachindex(res)
    d = load_object(data[ii])
    ln_D[ii] = linear_fit(log.(d[1])[1:data_pts], log.(d[2])[1:data_pts])[1]
end

colors = [CF_red, CF_vermillion, CF_orange, CF_yellow, CF_green, CF_sky, CF_blue]
set_theme!(CF_theme)

fig = Figure(size = (1200, 1600))

supertitle = fig[1:2, 1]
Label(
    supertitle,
    "Diffusivity",
    tellwidth = false,
    tellheight = false,
    font = :latex,
    fontsize = 44,
    valign = :center,
)

main_grid = fig[3:30, 1] = GridLayout()
# legend_grid = fig[31, 1] = GridLayout()

ax_autocorr = Axis(
    main_grid[1, 1],
    xscale = log10,
    yscale = log10,
    xlabel = "Time (ps)",
    ylabel = L"$\langle |\Delta \mathbf{R}|^2\rangle$ (\AA)",
)

ax_diffusivity = Axis(
    main_grid[2, 1],
    # title = "Velocity",
    ylabel = L"Velocity $\dot{R}$ (Å/ps)",
    xlabel = L"$1 / \hbar\Omega_T\,\mathrm{(1/meV)}$",
)
ax = [ax_autocorr, ax_diffusivity]
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

labs = ["10", "15", "20", "25", "30", "40", "50"]
for ii in eachindex(data)
    d = load_object(data[ii])

    scatter!(
        ax_autocorr,
        d[1],
        d[2],
        color = colors[ii],
        # marker = markers[ii],
        markersize = 14,
        label = labs[ii],
    )
    # lines!(ax, d[1], prefactor[ii] .* d[1], color = colors[ii], linewidth = 2)

end
axislegend(
    ax_autocorr,
    L"\hbar \Omega_T\,\mathrm{(meV)}",
    position = :lt,
    framevisible = false,
    orientation = :horizontal,
)

scatter!(
    ax_diffusivity,
    inv_T,
    ln_D,
    color = CF_sky,
    # linewidth = 4,
)
# lines!(
#     ax_speed,
#     δt .* (1:length(speed_full)),
#     [x[3] for x in speed_full],
#     color = CF_sky,
#     linewidth = 4,
# )

# lines!(
#     ax_pos,
#     δt .* (1:length(pos_local)),
#     [x[3] for x in pos_local],
#     color = CF_vermillion,
#     linewidth = 4,
# )
# lines!(
#     ax_speed,
#     δt .* (1:length(speed_local)),
#     [x[3] for x in speed_local],
#     color = CF_vermillion,
#     linewidth = 4,
# )
# hidexdecorations!(ax_pos)
# polys =
#     [PolyElement(color = c, strokecolor = :transparent) for c in [CF_sky, CF_vermillion]]

# Legend(
#     legend_grid[1, 1],
#     [polys],
#     [["Full solution", "Time-local solution"]],
#     [""],
#     halign = :center,
#     valign = :center,
#     tellheight = false,
#     tellwidth = false,
#     framevisible = false,
#     orientation = :horizontal,
#     titlevisible = false,
#     titleposition = :left,
#     titlefont = :latex,
# )

fig


inv_T = 1 ./ ħΩs

scatter(ln_D, inv_T)
ln_D
# ln r^2 = ln d + ln t
# ln d is the intercept
# d ~ c exp^[-E / T]
# ln d ~ ln c - E / T
# plot ln d vs 1 / T
# ln_D = ln c - S * (inv_T)^b
linear_fit(1 ./ ħΩs[3:end], (res)[3:end])
# curve_fit(LinearFit, ln_D, inv_T)
1 ./ ħΩs
