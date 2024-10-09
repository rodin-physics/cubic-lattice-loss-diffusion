include("lattice.jl")
# Function to calculate the effective positions of the framework and the mobile particle
function pos_eff(atoms, r, r_dot, R, R_dot, U, loss_mat)

    # Get ∇ᵣ assuming a pairwise interaction for a 
    # set of framework positions r and the particle position R
    gradient_r = [ForwardDiff.gradient(U, x .- R) for x in r]
    # The correction to nth position is given by ∑ₖ Wₙ₋ₖ grad_r[k]
    r_corr = [
        sum([
            get(W_dict, atoms[n] .- atoms[k], nothing) * gradient_r[k] for
            k in eachindex(atoms)
        ]) for n in eachindex(atoms)
    ]

    # The correction to R is given by loss_mat * ∇_R dU/dt = loss_mat * (∇_R (∇_r U ⋅ ṙ) + Hess(U)⋅Ṙ)
    R_corr =
        loss_mat * (
            ForwardDiff.gradient(
                R_ -> sum([
                    dot(ForwardDiff.gradient(U, r[n] .- R_), r_dot[n]) for n in eachindex(r)
                ]),
                R,
            ) + ForwardDiff.hessian(R -> sum([U(R - x) for x in r]), R) * R_dot
        ) |> real

    r_eff = r .- r_corr
    R_eff = R + R_corr

    return (r_eff, R_eff)
end

function derivatives_local(atoms, r, r_dot, R, R_dot, U, loss_mat)

    r_eff, R_eff = pos_eff(atoms, r, r_dot, R, R_dot, U, loss_mat)
    # Force is computed using effective positions
    res = (
        R_dot,
        -ForwardDiff.gradient(R_eff -> sum([U(R_eff - x) for x in r_eff]), R_eff) ./ M,
    )
    return res
end

# r and r_dot are functions of time that give the homogeneous framework trajectory for atoms
function RKstep_local(current_state, ζsηs, ζsηs_dot, qs, Ωs, atoms, U, t, δt, loss_mat)
    # Times used in RK routine for the time-dependent portion
    ts = [t, t + δt / 4, t + δt / 4, t + δt / 2, t + 3 * δt / 4, t + δt]

    # Compute rH for the times
    rHs = Vector{Tuple{Vector{Vector{Float64}},Vector{Vector{Float64}}}}(undef, 6)

    if isnothing(ζsηs) | isnothing(ζsηs_dot)
        for ii in eachindex(ts)
            res_pos = [atom .- [1, 1, 1] for atom in atoms] .* a |> real
            res_vel = [atom .- [1, 1, 1] for atom in atoms] .* 0 |> real

            rHs[ii] = (res_pos, res_vel)
        end
    else
        phases = [exp(1im * dot(q, atom)) for atom in atoms, q in qs]
        for ii in eachindex(ts)
            # Threads.@threads for ii in eachindex(ts)
            time_phase = exp.(-1im * ts[ii] .* Ωs)
            harmonics_pos = ζsηs .* time_phase
            harmonics_vel = ζsηs_dot .* time_phase
            res_pos =
                eachrow(phases * transpose(hcat(harmonics_pos...))) .+
                [atom .- [1, 1, 1] for atom in atoms] .* a |> real
            res_vel =
                eachrow(phases * transpose(hcat(harmonics_vel...))) .+
                [atom .- [1, 1, 1] for atom in atoms] .* 0 |> real

            rHs[ii] = (res_pos, res_vel)
        end
    end

    s1 = current_state
    k1 = derivatives_local(atoms, rHs[1]..., s1..., U, loss_mat)

    s2 = current_state .+ k1 .* (δt / 4)
    k2 = derivatives_local(atoms, rHs[2]..., s2..., U, loss_mat)

    s3 = current_state .+ (k1 .+ k2) .* (δt / 8)
    k3 = derivatives_local(atoms, rHs[3]..., s3..., U, loss_mat)

    s4 = current_state .+ k3 .* δt .- k2 .* (δt / 2)
    k4 = derivatives_local(atoms, rHs[4]..., s4..., U, loss_mat)

    s5 = current_state .+ k1 .* (δt * 3 / 16) .+ k4 .* (δt * 9 / 16)
    k5 = derivatives_local(atoms, rHs[5]..., s5..., U, loss_mat)

    s6 =
        current_state .- k1 .* (3 / 7 * δt) .+ k2 .* (2 / 7 * δt) .+ k3 .* (12 / 7 * δt) .-
        k4 .* (12 / 7 * δt) .+ k5 .* (8 / 7 * δt)
    k6 = derivatives_local(atoms, rHs[6]..., s6..., U, loss_mat)

    res =
        current_state .+
        (δt / 90) .* (7 .* k1 .+ 32 .* k3 .+ 12 .* k4 .+ 32 .* k5 .+ 7 .* k6)
    return res
end
