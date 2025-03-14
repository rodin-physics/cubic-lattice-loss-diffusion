include("lattice.jl")

function derivatives_periodic(m, M, a, index, V, U, r, rdot, R, Rdot)
    # Intrinsic framework acceleration
    F_r = -V * r

    pos = ceil.(Int, (R ./ a))
    atoms =
        Iterators.product([pos[1]:pos[1]+1, pos[2]:pos[2]+1, pos[3]:pos[3]+1]...) |>
        collect |>
        vec

    atoms = [mod1.(atom, size_x) for atom in atoms]

    # Get the indices of the coordinates of each atom in 'atoms' interacting with the mobile particle
    atom_indices = [
        [
            get(index, c, ErrorException("Coordinate not found")) for
            c in [(a..., d) for d in [1, 2, 3]]
        ] for a in atoms
    ]

    # Get the positions of each atom in 'atoms'. Atom (1,1,1) is at the origin
    atom_positions = [r[atom_indices[n]] .+ (atoms[n] .- 1) .* a for n in eachindex(atoms)]

    # Fold R back into the system ONLY FOR THE FORCES
    R = mod.(R, size_x * a) # R components can be between 0 and size_x * a
    F_r_particle = Vector{Vector{Float64}}(undef, length(atoms))
    for (i, atom_pos) in enumerate(atom_positions)
        disp = mod.(atom_pos .- R .+ size_x * a / 2, size_x * a) .- size_x * a / 2
        F_r_particle[i] = -ForwardDiff.gradient(U, disp)
    end
    # Force on mobile particles due to interaction using Newton's 3d law
    F_R = -sum(F_r_particle)
    view(F_r, vcat(atom_indices...)) .+= vcat(F_r_particle...)

    return (rdot, F_r ./ m, Rdot, F_R ./ M)
end

# 5th order Runge-Kutta step. 
# current_state is (r, rdot, R, Rdot)
# sys is (index, V)
function RKstep_periodic(m, M, a, sys, U, current_state, δt)
    k1 = derivatives_periodic(m, M, a, sys..., U, current_state...)
    k2 = derivatives_periodic(m, M, a, sys..., U, (current_state .+ k1 .* (δt / 4))...)
    k3 = derivatives_periodic(
        m,
        M,
        a,
        sys...,
        U,
        (current_state .+ (k1 .+ k2) .* (δt / 8))...,
    )
    k4 = derivatives_periodic(
        m,
        M,
        a,
        sys...,
        U,
        (current_state .+ k3 .* δt .- k2 .* (δt / 2))...,
    )
    k5 = derivatives_periodic(
        m,
        M,
        a,
        sys...,
        U,
        (current_state .+ k1 .* (δt * 3 / 16) .+ k4 .* (δt * 9 / 16))...,
    )
    k6 = derivatives_periodic(
        m,
        M,
        a,
        sys...,
        U,
        (
            current_state .- k1 .* (3 / 7 * δt) .+ k2 .* (2 / 7 * δt) .+
            k3 .* (12 / 7 * δt) .- k4 .* (12 / 7 * δt) .+ k5 .* (8 / 7 * δt)
        )...,
    )
    res =
        current_state .+
        (δt / 90) .* (7 .* k1 .+ 32 .* k3 .+ 12 .* k4 .+ 32 .* k5 .+ 7 .* k6)
    return res
end
