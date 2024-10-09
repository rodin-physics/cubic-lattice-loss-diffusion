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
