include("3D_Yukawa_diffusion.jl")
include("../src/plotting.jl")
using LsqFit
using Random
Random.seed!(250)
data = [
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT25_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT50_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U08000_λ0.5_ħΩT25_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U08000_λ0.5_ħΩT50_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT10_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT15_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT20_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT30_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT40_size20.jld2",
]

saves = [
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT25_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT50_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U08000_λ0.5_ħΩT25_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U08000_λ0.5_ħΩT50_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT10_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT15_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT20_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT30_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT40_size20_traj.jld2",
]
nBins = 10
nParts = 25
nTraj = nBins * nParts

for ii in eachindex(data)
    if !isfile(saves[ii])

        d = load_object(data[ii])
        pos = d[1]

        particle_traj = reshape(pos, :, nTraj)

        particle_traj =
            hcat([step .- particle_traj[1, :] for step in eachrow(particle_traj)]...)
        disp = norm.(particle_traj) .^ 2
        times = (1:size(disp)[2]) .* δt
        save_object(saves[ii], (times, disp))
    end

end

## MSD
colors = [CF_sky, CF_vermillion, CF_sky, CF_vermillion]
styles = [:solid, :solid, :dot, :dot]
labs = [
    L"$\hbar\Omega_T = 25$meV, $U_0 = 4$eV",
    L"$\hbar\Omega_T = 50$meV, $U_0 = 4$eV",
    L"$\hbar\Omega_T = 25$meV, $U_0 = 8$eV",
    L"$\hbar\Omega_T = 50$meV, $U_0 = 8$eV",
]
set_theme!(CF_theme)

fig = Figure(size = (1200, 800), figure_padding = 20)

ax = Axis(
    fig[1, 1],
    title = "Mean Squared Displacement",
    titlesize = 44,
    xlabel = "Time (ps)",
    ylabel = L"$\langle |\Delta \mathbf{R}|^2\rangle$ (\AA^2)",
)

for ii = 1:4
    d = load_object(saves[ii])
    times = d[1]
    nTraj = size(d[2], 1)
    MSD = [mean(d[2][:, step]) for step in eachindex(times)]
    STD = [std(d[2][:, step]) for step in eachindex(times)]
    ST_ERROR = STD ./ √(nTraj)

    lines!(
        ax,
        times,
        MSD,
        color = colors[ii],
        linestyle = styles[ii],
        linewidth = 4,
        label = labs[ii],
    )
    band!(ax, times, MSD .- ST_ERROR, MSD .+ ST_ERROR, color = colors[ii], alpha = 0.3)
    # LINEAR FIT
    function linear_fit(x, p)
        return (p[1] .* x .+ p[2])
    end
    function linear_fit(x, p)
        return (p[1] .* x)
    end
    p0 = [1.0]
    fit = LsqFit.curve_fit(linear_fit, times[200:end], MSD[200:end], p0)

    ps = fit.param
    lines!(ax, times, ps[1] .* times, linewidth = 2, color = CF_black, linestyle = :dash)

end
ylims!(ax, (0, 350))
xlims!(ax, (0, 12))
axislegend(ax, position = :lt, framevisible = false)
fig
save("MSD.pdf", fig)

## DIFFUSIVITY ESTIMATE
data = [
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT10_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT15_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT20_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT25_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT30_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT40_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT50_size20_traj.jld2",
]

ħΩs = [10, 15, 20, 25, 30, 40, 50]
inv_T = 1 ./ ħΩs
Ds = zeros(length(data))
Ds_Err = zeros(length(data))

colors = reverse([CF_red, CF_vermillion, CF_orange, CF_yellow, CF_green, CF_sky, CF_blue])
set_theme!(CF_theme)

fig = Figure(size = (1200, 800), figure_padding = 20)

ax_MSD = Axis(
    fig[1, 1],
    xlabel = "Time (ps)",
    ylabel = L"$\langle |\Delta \mathbf{R}|^2\rangle$ (\AA^2)",
    title = "Diffusion Coefficient",
    titlesize = 44,
)

ax_diffusivity = Axis(
    fig,
    bbox = BBox(240, 540, 400, 600),
    yscale = log10,
    xlabel = L"1 / \hbar\Omega_T\,\mathrm{(1/meV)}",
    ylabel = L"D\,\mathrm{(\AA^2/ps)}",
    xticks = WilkinsonTicks(3),
    # xaxisposition = :top,
    xlabelsize = 28,
    ylabelsize = 28,
    xticklabelsize = 28,
    yticklabelsize = 28,
)


labs = ["10", "15", "20", "25", "30", "40", "50"]
for ii in eachindex(data)
    d = load_object(data[ii])

    times = d[1]
    disp = d[2]
    nTraj = size(disp, 1)
    MSD = [mean(disp[:, step]) for step in eachindex(times)]
    STD = [std(disp[:, step]) for step in eachindex(times)]
    ST_ERROR = STD ./ √(nTraj)

    lines!(ax_MSD, times, MSD, color = colors[ii], linewidth = 4, label = labs[ii])
    band!(ax_MSD, times, MSD .- ST_ERROR, MSD .+ ST_ERROR, color = colors[ii], alpha = 0.3)

    # LINEAR FIT
    function linear_fit(x, p)
        return (p[1] .* x)
    end
    p0 = [1.0]
    fit = LsqFit.curve_fit(linear_fit, times[200:end], MSD[200:end], p0)

    ps = fit.param
    # println(standard_errors(fit))
    # println(confidence_interval(fit, 1))
    lines!(
        ax_MSD,
        times,
        ps[1] .* times,
        linewidth = 2,
        color = CF_black,
        linestyle = :dash,
    )
    Ds[ii] = ps[1] / 6

    # ERROR ESTIMATION

    disp = disp[shuffle(1:size(disp, 1)), :]

    bins = [disp[(1+(bin-1)*nParts):(bin*nParts), :] for bin = 1:nBins]
    MSDs = [sum(bin, dims = 1) ./ nParts for bin in bins]

    fits = [LsqFit.curve_fit(linear_fit, times[200:end], msd[200:end], p0) for msd in MSDs]
    ps = [f.param for f in fits]
    bin_Ds = [x[1] / 6 for x in ps]
    Ds_Err[ii] = std(bin_Ds) ./ sqrt(nBins)
    # Ds_Err[ii] = std(bin_Ds) ./ sqrt(nBins) * 1.96
end
axislegend(
    ax_MSD,
    L"\hbar \Omega_T\,\mathrm{(meV)}",
    position = :lt,
    framevisible = false,
    orientation = :horizontal,
)
ylims!(ax_MSD, (0, 350))
xlims!(ax_MSD, (0, 12))
scatter!(ax_diffusivity, inv_T, Ds, color = colors, markersize = 18)
errorbars!(
    ax_diffusivity,
    inv_T,
    Ds,
    [min(Ds_Err[ii], Ds[ii]) - 1e-4 for ii in eachindex(Ds)],
    Ds_Err,
    color = colors,
    whiskerwidth = 10,
)
ylims!(ax_diffusivity, (1e-2, 10))


function arrhenius(x, p)
    return p[1] * exp.(-(p[2] .* x))
end
p0 = [1.0, 1.0]
fit = LsqFit.curve_fit(arrhenius, inv_T[1:end], Ds[1:end], p0)

xs = range(minimum(inv_T) * 0.95, maximum(inv_T) * 1.05, length = 100)

ps = fit.param
lines!(
    ax_diffusivity,
    xs,
    [arrhenius(x, ps) for x in xs],
    color = CF_black,
    linewidth = 4,
    linestyle = :dash,
)
fig
save("Diffusivity.pdf", fig)

# ps
# fit

# confidence_interval(fit, 0.05)[1] .- ps[1]
# confidence_interval(fit, 0.05)[2] .- ps[2]
