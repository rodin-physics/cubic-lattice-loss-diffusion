include("3D_Gaussian_diffusion.jl")
include("../src/plotting.jl")

data = [
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT25_size20.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT50_size20.jld2",
    "Data/Diffusion/Diffusion_U02000_λ0.75_ħΩT25_size20.jld2",
    "Data/Diffusion/Diffusion_U02000_λ0.75_ħΩT50_size20.jld2",
]

data_reduced = [
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT25_size20_Reduced.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT50_size20_Reduced.jld2",
    "Data/Diffusion/Diffusion_U02000_λ0.75_ħΩT25_size20_Reduced.jld2",
    "Data/Diffusion/Diffusion_U02000_λ0.75_ħΩT50_size20_Reduced.jld2",
]

for jj = 1:4
    if !isfile(data_reduced[jj])
        d = load_object(data[jj])
        p = d[1]
        v = d[2]

        p = p[1:10:end]
        v = v[1:10:end]

        save_object(data_reduced[jj], (p, v))
    end
end

## FIGURES
colors = [CF_sky, CF_vermillion, CF_sky, CF_vermillion]
set_theme!(CF_theme)
fig = Figure(size = (1200, 1200))
ax_cold_weak = Axis3(
    fig[1, 1],
    aspect = :data,
    title = L"$U_0 = 1000$meV, $\hbar\Omega_T = 25$meV",
    # xlabel = L"$x$ \, (\AA)",
    # ylabel = L"$y$ \, (\AA)",
    # zlabel = L"$z$ \, (\AA)",
    xlabelvisible = false,
    ylabelvisible = false,
    zlabelvisible = false,
    xticks = WilkinsonTicks(3),
    yticks = WilkinsonTicks(3),
    zticks = WilkinsonTicks(3),
    # xlabeloffset = 75
)

ax_hot_weak = Axis3(
    fig[1, 2],
    aspect = :data,
    title = L"$U_0 = 1000$meV, $\hbar\Omega_T = 50$meV",
    # xlabel = L"$x$ \, (\AA)",
    # ylabel = L"$y$ \, (\AA)",
    # zlabel = L"$z$ \, (\AA)",
    xlabelvisible = false,
    ylabelvisible = false,
    zlabelvisible = false,
    xticks = WilkinsonTicks(3),
    yticks = WilkinsonTicks(3),
    zticks = WilkinsonTicks(3),
    # xlabeloffset = 75
)

ax_cold_strong = Axis3(
    fig[2, 1],
    aspect = :data,
    title = L"$U_0 = 2000$meV, $\hbar\Omega_T = 25$meV",
    # xlabel = L"$x$ \, (\AA)",
    # ylabel = L"$y$ \, (\AA)",
    # zlabel = L"$z$ \, (\AA)",
    xlabelvisible = false,
    ylabelvisible = false,
    zlabelvisible = false,
    xticks = WilkinsonTicks(3),
    yticks = WilkinsonTicks(3),
    zticks = WilkinsonTicks(3),
    # xlabeloffset = 75
)

ax_hot_strong = Axis3(
    fig[2, 2],
    aspect = :data,
    title = L"$U_0 = 2000$meV, $\hbar\Omega_T = 50$meV",
    # xlabel = L"$x$ \, (\AA)",
    # ylabel = L"$y$ \, (\AA)",
    # zlabel = L"$z$ \, (\AA)",
    xlabelvisible = false,
    ylabelvisible = false,
    zlabelvisible = false,
    xticks = WilkinsonTicks(3),
    yticks = WilkinsonTicks(3),
    zticks = WilkinsonTicks(3),
    # xlabeloffset = 75
)
colgap!(fig.layout, Relative(0.05))
rowgap!(fig.layout, Relative(0.1))

axs = [ax_cold_weak, ax_hot_weak, ax_cold_strong, ax_hot_strong]
for ii = 1:4
    d = load_object(data_reduced[ii])
    p = d[1]
    lines!(
        axs[ii],
        [r[1] for r in p],
        [r[2] for r in p],
        [r[3] for r in p],
        color = colors[ii],
    )

end

fig
save("Diffusion.pdf", fig)


## Histograms of speed
d = load_object("Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT25_size20.jld2")

## FIGURES
set_theme!(CF_theme)
fig = Figure(size = (1200, 800))

ax = Axis(
    fig[1, 1],
    # aspect = :data,
    # title = L"$U_0 = 1000$meV, $\hbar\Omega_T = 25$meV",
    # xlabel = L"$x$ \, (\AA)",
    # ylabel = L"$y$ \, (\AA)",
)

for ii = 1:4
    d = load_object(data[ii])

    pos_data = d[1]
    speed_data = d[2]

    # pos_sq = pos_data .^ 2
    # speed_sq = speed_data .^ 2
    hist!(
        ax,
        norm.(speed_data),
        bins = 100,
        normalization = :pdf,
        color = alphacolor(colors[ii], 0.5),
    )

end
# colors
fig

## AUTOCORRELATION
nSteps = floor(Int, t_max / δt / 2)
skip = 1000

nPts = floor(Int, nSteps / skip)

times = skip .* (1:nPts) .* δt

saves = [
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT25_size20_Autocorr.jld2",
    "Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT50_size20_Autocorr.jld2",
    "Data/Diffusion/Diffusion_U02000_λ0.75_ħΩT25_size20_Autocorr.jld2",
    "Data/Diffusion/Diffusion_U02000_λ0.75_ħΩT50_size20_Autocorr.jld2",
]

for ii = 1:4
    d = load_object(data[ii])
    loc = d[1]

    res = zeros(nPts)

    @showprogress for time = 1:nPts
        res[time] = (norm.(loc[(1+skip*time):end] - loc[1:(end-skip*time)])) .^ 2 |> mean
    end
    save_object(saves[ii], (times, res))

end

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
    # aspect = :data,
    title = "Autocorrelation",
    xscale = log10,
    yscale = log10,
    xlabel = L"$t$ \, (ps)",
    ylabel = L"$\langle R^2\rangle$ (\AA)",
)

for ii = 1:4
    d = load_object(saves[ii])

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
save("Autocorrelation.pdf", fig)


# d = load_object("Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT50_size20.jld2")
# loc = d[1]

# norm.(loc[300001:end] - loc[1:end-300000]) |> mean
# ST = 10






# speed = d[2]
# (norm.(speed) .^ 2 |> maximum) / 2 * 0.7
# hist(norm.(speed))
# U(1.3) - 8U(sqrt(3) * (1.5))



# data = [load_object("Data/System/Homogeneous_ħT$(ħΩT).jld2") for ħΩT in ħΩTs]
# σ_pos = [Corr_pos(ħ, m, DynamicalMatrix, ħΩT / ħ, 0, [0, 0, 0])[1][1, 1] for ħΩT in ħΩTs]
# σ_speed =
#     [Corr_speed(ħ, m, DynamicalMatrix, ħΩT / ħ, 0, [0, 0, 0])[1][1, 1] for ħΩT in ħΩTs]
# ΩT_label = [L"10^{-3}", L"5", L"25", L"50"]
# colors = [CF_red, CF_yellow, CF_green, CF_blue]

# for jj in eachindex(data)
#     pos_data = data[jj][1]
#     speed_data = data[jj][2]

#     pos_sq = pos_data .^ 2
#     speed_sq = speed_data .^ 2

#     N = size_x * size_y * size_z

#     r = [sum(pos_sq[(3*ii+1):(3*ii+3)]) for ii = 0:(N-1)]
#     v = [sum(speed_sq[(3*ii+1):(3*ii+3)]) for ii = 0:(N-1)]

#     hist!(
#         ax_pos,
#         sqrt.(r),
#         bins = 100,
#         normalization = :pdf,
#         color = alphacolor(colors[jj], 0.5),
#     )
#     hist!(
#         ax_speed,
#         sqrt.(v),
#         bins = 100,
#         normalization = :pdf,
#         color = alphacolor(colors[jj], 0.5),
#     )

#     plotPts = 500
#     pos = range(0, 1, length = plotPts)
#     speed = range(0, 15, length = plotPts)

#     lines!(
#         ax_pos,
#         pos,
#         4 * π .* pos .^ 2 ./ √(2 * π * σ_pos[jj])^3 .* exp.(-pos .^ 2 ./ 2 ./ σ_pos[jj]),
#         color = colors[jj],
#         linewidth = 4,
#     )

#     lines!(
#         ax_speed,
#         speed,
#         4 * π .* speed .^ 2 ./ √(2 * π * σ_speed[jj])^3 .*
#         exp.(-speed .^ 2 ./ 2 ./ σ_speed[jj]),
#         color = colors[jj],
#         linewidth = 4,
#         label = ΩT_label[jj],
#     )
# end
# xlims!(ax_pos, (0, 0.749))
# xlims!(ax_speed, (0, 14.9))
# ylims!(ax_pos, (0, 9))
# ylims!(ax_speed, (0, 1 / 2))
# hideydecorations!(ax_pos, label = false)
# hideydecorations!(ax_speed)
# axislegend(L"\hbar\Omega_T (\mathrm{meV})", orientation = :horizontal, framevisible = false)
d = load_object("Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT50_size20.jld2")
loc = d[1]

norm.(loc[300001:end] - loc[1:end-300000]) |> mean
ST = 10


tst = [(norm.(loc[(1+10000*ST):end] - loc[1:(end-10000*ST)])) .^ 2 |> mean for ST = 1:30]

scatter(10000 .* (1:30) .* 5e-3, sqrt.(tst))
