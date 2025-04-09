include("3D_Yukawa_diffusion.jl")
include("3D_Yukawa_diffusion_full_periodic.jl")
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

data = [
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT10_size20.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT15_size20.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT20_size20.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT25_size20.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT30_size20.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT40_size20.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT50_size20.jld2",
]

saves = [
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT10_size20_traj.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT15_size20_traj.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT20_size20_traj.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT25_size20_traj.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT30_size20_traj.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT40_size20_traj.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT50_size20_traj.jld2",
]

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

data = [
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT25_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT50_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U08000_λ0.5_ħΩT25_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U08000_λ0.5_ħΩT50_size20_traj.jld2",
]
colors = [CF_sky, CF_vermillion, CF_sky, CF_vermillion]
styles = [:solid, :solid, :dash, :dash]
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
    d = load_object(data[ii])
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
    band!(ax, times, MSD .- ST_ERROR, MSD .+ ST_ERROR, color = colors[ii], alpha = 0.6)

    # LINEAR FIT

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
data_time_local = [
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT10_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT15_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT20_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT25_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT30_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT40_size20_traj.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT50_size20_traj.jld2",
]

data_full = [
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT10_size20_traj.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT15_size20_traj.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT20_size20_traj.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT25_size20_traj.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT30_size20_traj.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT40_size20_traj.jld2",
    "Data/Diffusion_Full_Periodic/Yukawa_Diffusion_Full_Periodic_U04000_λ0.5_ħΩT50_size20_traj.jld2",
]

ħΩs = [10, 15, 20, 25, 30, 40, 50]
inv_T = 1 ./ ħΩs
Ds = zeros(2, length(ħΩs))
Ds_Err = zeros(2, length(ħΩs))

colors = reverse([CF_red, CF_vermillion, CF_orange, CF_yellow, CF_green, CF_sky, CF_blue])
set_theme!(CF_theme)

fig = Figure(size = (1200, 1600), figure_padding = 20)

ax_MSD_time_local = Axis(
    fig[1, 1],
    xlabel = "Time (ps)",
    ylabel = L"Mean squared displacement $\langle |\Delta \mathbf{R}|^2\rangle$ (\AA^2)",
    title = "Diffusion Coefficient Temperature Dependence",
    titlesize = 44,
    yticks = 0:50:350,
)

ax_MSD_full = Axis(
    fig[2, 1],
    xlabel = "Time (ps)",
    ylabel = L"Mean squared displacement $\langle |\Delta \mathbf{R}|^2\rangle$ (\AA^2)",
    # title = "Diffusion Coefficient",
    titlesize = 44,
    yticks = 0:150:900,
)


ax_diffusivity_time_local = Axis(
    fig[1, 1],
    width = Relative(0.3),
    height = Relative(0.3),
    halign = 0.2,
    valign = 0.75,
    yscale = log10,
    xlabel = L"1 / \hbar\Omega_T\,\mathrm{(meV^{-1})}",
    ylabel = L"D\,\mathrm{(\AA^2/ps)}",
    xticks = WilkinsonTicks(3),
    yticks = ([0.1,1,10], ["0.1", "1", "10"]),
    xlabelsize = 28,
    ylabelsize = 28,
    xticklabelsize = 28,
    yticklabelsize = 28,
)

ax_diffusivity_full = Axis(
    fig[2, 1],
    width = Relative(0.3),
    height = Relative(0.3),
    halign = 0.2,
    valign = 0.75,
    yscale = log10,
    xlabel = L"1 / \hbar\Omega_T\,\mathrm{(meV^{-1})}",
    ylabel = L"D\,\mathrm{(\AA^2/ps)}",
    xticks = WilkinsonTicks(3),
    yticks = ([0.1,1,10], ["0.1", "1", "10"]),
    xlabelsize = 28,
    ylabelsize = 28,
    xticklabelsize = 28,
    yticklabelsize = 28,
)

labs = ["10", "15", "20", "25", "30", "40", "50"]
ax = [ax_MSD_time_local, ax_MSD_full]
insets = [ax_diffusivity_time_local, ax_diffusivity_full]
for ii in eachindex(data_full)
    ds = [load_object(data_time_local[ii]), load_object(data_full[ii])]

    for jj = 1:2
        d = ds[jj]
        times = d[1]
        disp = d[2]
        nTraj = size(disp, 1)
        MSD = [mean(disp[:, step]) for step in eachindex(times)]
        STD = [std(disp[:, step]) for step in eachindex(times)]
        ST_ERROR = STD ./ √(nTraj)

        lines!(ax[jj], times, MSD, color = colors[ii], linewidth = 4, label = labs[ii])
        band!(
            ax[jj],
            times,
            MSD .- ST_ERROR,
            MSD .+ ST_ERROR,
            color = colors[ii],
            alpha = 0.6,
        )

        # LINEAR FIT
        function linear_fit(x, p)
            return (p[1] .* x)
        end
        p0 = [1.0]
        fit = LsqFit.curve_fit(linear_fit, times[200:end], MSD[200:end], p0)

        ps = fit.param
        lines!(
            ax[jj],
            times,
            ps[1] .* times,
            linewidth = 2,
            color = CF_black,
            linestyle = :dash,
        )
        Ds[jj, ii] = ps[1] / 6

        # ERROR ESTIMATION

        disp = disp[shuffle(1:size(disp, 1)), :]

        bins = [disp[(1+(bin-1)*nParts):(bin*nParts), :] for bin = 1:nBins]
        MSDs = [sum(bin, dims = 1) ./ nParts for bin in bins]

        fits =
            [LsqFit.curve_fit(linear_fit, times[200:end], msd[200:end], p0) for msd in MSDs]
        ps = [f.param for f in fits]
        bin_Ds = [x[1] / 6 for x in ps]
        Ds_Err[jj, ii] = std(bin_Ds) ./ sqrt(nBins)
    end

end

axislegend(
    ax_MSD_time_local,
    L"\hbar \Omega_T\,\mathrm{(meV)}",
    position = :lt,
    framevisible = false,
    orientation = :horizontal,
)
for jj = 1:2
    scatter!(insets[jj], inv_T, Ds[jj, :], color = colors, markersize = 18)
    errorbars!(
        insets[jj],
        inv_T,
        Ds[jj, :],
        [min(Ds_Err[jj, ii], Ds[jj, ii]) - 1e-4 for ii in eachindex(Ds[jj, :])],
        Ds_Err[jj, :],
        color = colors,
        whiskerwidth = 10,
    )
    p0 = [1.0, 1.0]
    function arrhenius(x, p)
        return (-p[2] .* x .+ p[1])
    end
    fit = LsqFit.curve_fit(
        arrhenius,
        inv_T,
        log.(Ds[jj, :]),
        (Ds[jj, :] ./ Ds_Err[jj, :]) .^ 2,
        p0,
    )
    xs = range(minimum(inv_T) * 0.95, maximum(inv_T) * 1.05, length = 100)
    ps = fit.param
    lines!(
        insets[jj],
        xs,
        [exp.(arrhenius(x, ps)) for x in xs],
        color = CF_black,
        linewidth = 4,
        linestyle = :dash,
    )
    xlims!(ax[jj], (0, 12))
    print(ps)
    print(confidence_interval(fit, 0.05)[1] .- ps[1])
    print(confidence_interval(fit, 0.05)[2] .- ps[2])
end
# hidexdecorations(ax[1])
ylims!(ax[1], (0, 350))
ylims!(ax[2], (0, 1000))
ylims!(insets[1], (0.1, 30))
ylims!(insets[2], (0.1, 30))
labs = ["(a)", "(b)"]
for ii = 1:2
    text!(
        ax[ii],
        0.06,
        0.98,
        text = labs[ii],
        align = (:right, :top),
        fontsize = 36,
        font = :latex,
        color = :black,
        space = :relative,
    )

end

fig
save("Diffusivity.pdf", fig)
# exp(2.279622128837239+0.21740460475056222)
# exp(2.279622128837239-0.21740460475056222)
# exp(3.1353353267937805 + 0.17528142514903644)
# exp(3.1353353267937805 - 0.17528142514903644)

# exp(2.279622128837239+0)
# exp(3.1353353267937805)



# exp(2.279622128837239+0.21740460475056222) - exp(2.279622128837239)
# exp(2.279622128837239-0.21740460475056222)- exp(2.279622128837239)


# exp(3.1353353267937805 + 0.17528142514903644)-exp(3.1353353267937805 )
# exp(3.1353353267937805 - 0.17528142514903644)-exp(3.1353353267937805 )