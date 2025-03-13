include("3D_Yukawa_diffusion_full.jl")
include("../src/plotting.jl")
using LsqFit
using Random
Random.seed!(250)

data = [
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT10_size30.jld2",
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT15_size30.jld2",
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT20_size30.jld2",
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT25_size30.jld2",
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT30_size30.jld2",
    # "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT40_size30.jld2",
    # "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT50_size30.jld2",
]

saves = [
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT10_size30_traj.jld2",
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT15_size30_traj.jld2",
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT20_size30_traj.jld2",
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT25_size30_traj.jld2",
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT30_size30_traj.jld2",
    # "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT40_size30_traj.jld2",
    # "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT50_size30_traj.jld2",
]


nBins = 10
nParts = 25
nTraj = nBins * nParts
## Prepare the trajectory data
for ii in eachindex(data)
    # Get the trajectory
    if !isfile(saves[ii])
        d = load_object(data[ii])[1]
        # Get the displacement arrays 
        disp = hcat([[norm(pt .- traj[1]) for pt in traj] for traj in d]...)
        times = (1:size(disp)[1]) .* δt
        save_object(saves[ii], (times, disp))
    end
end

## DIFFUSIVITY ESTIMATE
data = [
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT10_size30_traj.jld2",
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT15_size30_traj.jld2",
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT20_size30_traj.jld2",
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT25_size30_traj.jld2",
    "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT30_size30_traj.jld2",
    # "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT40_size30_traj.jld2",
    # "Data/Diffusion_Full/Yukawa_Diffusion_Full_U04000_λ0.5_ħΩT50_size30_traj.jld2",
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
    MSD = [mean(disp[step, :]) for step in eachindex(times)]
    STD = [std(disp[step, :]) for step in eachindex(times)]
    ST_ERROR = STD ./ √(nTraj)

    lines!(ax_MSD, times, MSD, color = colors[ii], linewidth = 4, label = labs[ii])
    band!(ax_MSD, times, MSD .- ST_ERROR, MSD .+ ST_ERROR, color = colors[ii], alpha = 0.3)

    # LINEAR FIT
    function linear_fit(x, p)
        return (p[1] .* x .+ p[2])
    end
    p0 = [1.0, 0.0]
    fit = LsqFit.curve_fit(linear_fit, times[800:end], MSD[800:end], p0)

    ps = fit.param
    # println(standard_errors(fit))
    # println(confidence_interval(fit, 1))
    lines!(
        ax_MSD,
        (times),
        ps[1] .* times .+ ps[2],
        linewidth = 2,
        color = CF_black,
        linestyle = :dash,
    )
    Ds[ii] = ps[1] / 6

    # ERROR ESTIMATION

    # disp = disp[:,shuffle(1:size(disp, 2))]

    # bins = [disp[(1+(bin-1)*nParts):(bin*nParts), :] for bin = 1:nBins]
    # MSDs = [sum(bin, dims = 1) ./ nParts for bin in bins]

    # fits = [LsqFit.curve_fit(linear_fit, times[200:end], msd[200:end], p0) for msd in MSDs]
    # ps = [f.param for f in fits]
    # bin_Ds = [x[1] / 6 for x in ps]
    # Ds_Err[ii] = std(bin_Ds) ./ sqrt(nBins)
    # Ds_Err[ii] = std(bin_Ds) ./ sqrt(nBins) * 1.96
end
axislegend(
    ax_MSD,
    L"\hbar \Omega_T\,\mathrm{(meV)}",
    position = :lt,
    framevisible = false,
    orientation = :horizontal,
)
ylims!(ax_MSD, (0, 10))
xlims!(ax_MSD, (0, 6))
scatter!(ax_diffusivity, inv_T[1:5], Ds, color = colors[1:5], markersize = 18)
# errorbars!(
#     ax_diffusivity,
#     inv_T,
#     Ds,
#     [min(Ds_Err[ii], Ds[ii]) - 1e-4 for ii in eachindex(Ds)],
#     Ds_Err,
#     color = colors,
#     whiskerwidth = 10,
# )
# ylims!(ax_diffusivity, (1e-2, 10))


function arrhenius(x, p)
    return p[1] * exp.(-(p[2] .* x))
end
p0 = [1.0, 1.0]
fit = LsqFit.curve_fit(arrhenius, inv_T[1:5], Ds[1:5], p0)

# xs = range(minimum(inv_T) * 0.95, maximum(inv_T) * 1.05, length = 100)

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
save("Diffusivity_Full.pdf", fig)

# ps
# fit

# confidence_interval(fit, 0.05)[1] .- ps[1]
# confidence_interval(fit, 0.05)[2] .- ps[2]
