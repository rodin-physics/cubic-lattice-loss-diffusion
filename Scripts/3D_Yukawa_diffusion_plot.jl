include("3D_Yukawa_diffusion.jl")
include("../src/plotting.jl")

data = [
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT10_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT15_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT20_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT25_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT30_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT40_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT50_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U08000_λ0.5_ħΩT25_size20.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U08000_λ0.5_ħΩT50_size20.jld2",
]

data_reduced = [
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT10_size20_Reduced.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT15_size20_Reduced.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT20_size20_Reduced.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT25_size20_Reduced.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT30_size20_Reduced.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT40_size20_Reduced.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT50_size20_Reduced.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U08000_λ0.5_ħΩT25_size20_Reduced.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U08000_λ0.5_ħΩT50_size20_Reduced.jld2",
]

for jj in eachindex(data)
    if !isfile(data_reduced[jj])
        d = load_object(data[jj])
        p = d[1]
        v = d[2]

        p = p[1:10:end]
        v = v[1:10:end]

        save_object(data_reduced[jj], (p, v))
    end
end


data_reduced = [
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT25_size20_Reduced.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U04000_λ0.5_ħΩT50_size20_Reduced.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U08000_λ0.5_ħΩT25_size20_Reduced.jld2",
    "Data/Diffusion/Yukawa_Diffusion_U08000_λ0.5_ħΩT50_size20_Reduced.jld2",
]

## FIGURES
colors = [CF_sky, CF_vermillion, CF_sky, CF_vermillion]
colors_markers = reverse(colors)
set_theme!(CF_theme)
fig = Figure(size = (1200, 1200))

supertitle = fig[1:2, 1]
Label(
    supertitle,
    "Particle Diffusion",
    tellwidth = false,
    tellheight = false,
    font = :latex,
    fontsize = 44,
    valign = :center,
)

main_grid = fig[3:30, 1] = GridLayout()

ax_cold_weak = Axis3(
    main_grid[1, 1],
    aspect = :data,
    azimuth = 3pi / 4,
    title = L"$U_0 = 4$eV, $\hbar\Omega_T = 25$meV",
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
    main_grid[1, 2],
    aspect = :data,
    azimuth = 3pi / 4,
    title = L"$U_0 = 4$eV, $\hbar\Omega_T = 50$meV",
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
    main_grid[2, 1],
    aspect = :data,
    azimuth = 3pi / 4,
    title = L"$U_0 = 8$eV, $\hbar\Omega_T = 25$meV",
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
    main_grid[2, 2],
    aspect = :data,
    azimuth = 3pi / 4,
    title = L"$U_0 = 8$eV, $\hbar\Omega_T = 50$meV",
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

colsize!(main_grid, 1, Relative(0.45))
colsize!(main_grid, 2, Relative(0.45))
# colgap!(main_grid, Relative(0.05))
# rowgap!(main_grid, Relative(0.05))

axs = [ax_cold_weak, ax_hot_weak, ax_cold_strong, ax_hot_strong]
for ii = 1:4
    d = load_object(data_reduced[ii])
    p = d[1]
    xs = [r[1] for r in p]
    ys = [r[2] for r in p]
    zs = [r[3] for r in p]


    # SPLIT THE PATH INTO SMALLER SEGMENTS FOR ILLUSTRATOR
    nPaths = 10
    per_path = Int(length(p) / nPaths)
    for jj = 1:10
        lines!(
            axs[ii],
            xs[(jj-1)*per_path+1:jj*per_path],
            ys[(jj-1)*per_path+1:jj*per_path],
            zs[(jj-1)*per_path+1:jj*per_path],
            color = colors[ii],
        )
        # lines!(axs[ii], xs, ys, zs, color = colors[ii])

    end

    lines!(axs[ii], [0, 0, 10, 10, 0], [0, 10, 10, 0, 0], [0, 0, 0, 0, 0], color = CF_black)
    lines!(
        axs[ii],
        [0, 0, 10, 10, 0],
        [0, 10, 10, 0, 0],
        [10, 10, 10, 10, 10],
        color = CF_black,
    )

    lines!(axs[ii], [0, 0], [0, 0], [10, 0], color = CF_black)
    lines!(axs[ii], [10, 10], [0, 0], [10, 0], color = CF_black)
    lines!(axs[ii], [0, 0], [10, 10], [10, 0], color = CF_black)
    lines!(axs[ii], [10, 10], [10, 10], [10, 0], color = CF_black)


    # xlims!(axs[ii], (minimum(xs), maximum(xs)))
    # xlabel!(axs[ii], floor(Int, maximum(xs)-minimum(xs)))
    # ylims!(axs[ii], (minimum(ys), maximum(ys)))
    # zlims!(axs[ii], (minimum(zs), maximum(zs)))
    scatter!(
        axs[ii],
        [xs[end]],
        [ys[end]],
        [zs[end]],
        markersize = 16,
        color = colors_markers[ii],
        marker = :cross,
    )

    scatter!(
        axs[ii],
        [xs[1]],
        [ys[1]],
        [zs[1]],
        markersize = 16,
        color = colors_markers[ii],
    )
    hidedecorations!(axs[ii])
end

fig
save("Yukawa_Diffusion.pdf", fig)
