include("lattice.jl")
# Function to calculate the effective positions of the framework and the mobile particle
function pos_eff(atoms, r, r_dot, R, R_dot, U, loss_mat)

    # Get ∇ᵣ assuming a pairwise interaction for a 
    # set of framework positions r and the particle position R
    # The correction to nth position is given by ∑ₖ Wₙ₋ₖ grad_r[k]
    for ii = 1:10
        δ =
            -[
                sum([
                    get(W_dict, atoms[n] .- atoms[k], nothing) *
                    ForwardDiff.gradient(U, r[k] .- R) for k in eachindex(atoms)
                ]) for n in eachindex(atoms)
            ]
        r = r .+ δ

    end

    Hess_r = [ForwardDiff.hessian(U, x .- R) for x in r]
    Hess_R = sum(Hess_r)
    Grad_r_R = -hcat(Hess_r...)
    Grad_R_r = Grad_r_R'

    Hess_r = Matrix(blockdiag(sparse.(Hess_r)...))

    V_inv = [
        get(W_dict, atoms[n] .- atoms[k], nothing) for k in eachindex(atoms),
        n in eachindex(atoms)
    ]

    V_inv = hcat([vcat(row...) for row in eachrow(V_inv)]...)
    V_inv = (V_inv + V_inv') / 2

    δ_dot =
        -inv(I(length(atoms) * 3) + V_inv * Hess_r) *
        V_inv *
        (Hess_r * vcat(r_dot...) + Grad_R_r * R_dot)


    Δ = loss_mat * (Hess_R * R_dot + Grad_r_R * (vcat(r_dot...) + δ_dot))

    # return (Δ, R_corr)
    return (r, R + Δ)
end
# R = [0,0,10.5]
# pos_eff(atoms, rH, rH_dot, R, R_dot, U, loss_mat)
# pos_eff_TST(atoms, rH, rH_dot, R, R_dot, U, loss_mat)
# pos_eff(atoms, rH, rH_dot, R, R_dot, U, loss_mat)[2] ==
# pos_eff_TST(atoms, rH, rH_dot, R, R_dot, U, loss_mat)[2]
# pos_eff_TST(atoms, rH, rH_dot, R, R_dot, U, loss_mat)[1] ==
# pos_eff_TST(atoms, rH, rH_dot, R, R_dot, U, loss_mat)[2]
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
