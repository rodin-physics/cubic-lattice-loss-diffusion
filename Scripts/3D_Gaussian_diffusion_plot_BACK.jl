include("3D_Gaussian_diffusion.jl")

ħΩTs = [25]

## FIGURES
set_theme!(CF_theme)
fig = Figure(size = (1200, 400))


ax_pos = Axis(
    fig[1, 1],
    title = "Displacement distribution",
    xlabel = L"r \, (\AA)",
    ylabel = "arb. units.",
)

ax_speed =
    Axis(fig[1, 2], title = "Speed distribution", xlabel = L"\dot{r}\,(\AA/\mathrm{ps})")

## THERMAL STATISTICS PLOTTING
"Data/Diffusion/Diffusion_U01000_λ0.75_ħT25.0_size20.jld2"
data = [
    load_object("Data/Diffusion/Diffusion_U0$(U0)_λ$(λ)_ħT$(ħΩT)_size$(size_x).jld2")
    for ħΩT in ħΩTs
]
# data = [load_object("Data/Diffusion/Diffusion_U0$(U0)_λ$(λ)_ħT$(ħΩT)_size$(size_x).jld2") for ħΩT in ħΩTs]
σ_pos = [Corr_pos(ħ, m, dyn_mat, ħΩT / ħ, 0, [0, 0, 0])[1][1, 1] for ħΩT in ħΩTs]
σ_speed = [Corr_speed(ħ, m, dyn_mat, ħΩT / ħ, 0, [0, 0, 0])[1][1, 1] for ħΩT in ħΩTs]
ΩT_label = [L"10^{-3}", L"5", L"25", L"50"]
colors = [CF_red, CF_yellow, CF_green, CF_blue]

for jj in eachindex(data)
    pos_data = data[jj][1]
    speed_data = data[jj][2]

    hist!(
        ax_pos,
        (pos_data),
        bins = 100,
        normalization = :pdf,
        color = alphacolor(colors[jj], 0.5),
    )
    hist!(
        ax_speed,
        (speed_data),
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

# ## PANEL LABELS
# labs = ["(b)", "(c)", "(d)", "(e)"]
# ax = [ax_bands, ax_recoil, ax_pos, ax_speed]
# for jj in eachindex(ax)
#     text!(
#         ax[jj],
#         0.11,
#         0.975,
#         text = labs[jj],
#         align = (:right, :top),
#         space = :relative,
#         fontsize = 36,
#         font = :latex,
#         color = :black,
#     )

# end
fig

# save("Framework_properties.pdf", fig)


res = load_object("Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT25_size20.jld2")
lines(res[1])



# # # # # # lines([r[1] for r in v],[r[2] for r in v],[r[3] for r in v])
# # acc = ((v[2:end].-v[1:end-1])./δt)[1:end]
# # vv = norm.(v)
# # hist(vv)
# # # # # # # # scatter([r[1] for r in acc],[r[2] for r in acc],[r[3] for r in acc])
# # lines(norm.(acc))
# # lines(norm.(p))
# # # # # lines(norm.(v))
# # lines(p)
# # # # # # # (norm.(v).^2/2*0.7)[10:20]
# # # # # # U([a /2,a/2,0a/2])
# # # # vcat(acc...)|>minimum



## FIGURES
set_theme!(CF_theme)
fig = Figure(size = (1200, 400))


ax = Axis3(
    fig[1, 1],
    title = "Diffusion",
    titlefont = :latex,

    # xlabel = L"",
    # ylabel = L"y \, (\AA)",
    # zlabel = L"z \, (\AA)",
    xticklabelfont = :latex,
    yticklabelfont = :latex,
    zticklabelfont = :latex,
    xlabelfont = :latex,
    ylabelfont = :latex,
    zlabelfont = :latex,
    aspect = :data,
)
#     xgridvisible = false,
#     ygridvisible = false,
#     xlabelfont = :latex,
#     ylabelfont = :latex,
#     xticklabelfont = :latex,
#     yticklabelfont = :latex,
#     xticklabelsize = 36,
#     yticklabelsize = 36,
#     titlesize = 36,
#     xlabelsize = 36,
#     ylabelsize = 36,
# )

hidedecorations!(ax, ticks = false, ticklabels = false, grid = false)
rr = load_object("Data/Diffusion/Diffusion_U01000_λ0.75_ħΩT25_size20.jld2")
p = rr[1]
lines!(ax, [r[1] for r in p], [r[2] for r in p], [r[3] for r in p], color = CF_vermillion)
fig

v = norm.(rr[2])

hist(v, bins = 200)
# Base.summarysize(rr)
