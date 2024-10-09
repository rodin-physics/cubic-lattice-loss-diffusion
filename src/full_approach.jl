include("lattice.jl")
# A function to produce a coupling matrix for a cubic lattice
function system(size_x::Int, size_y::Int, size_z::Int, coupling::Vector{Coupling})
    # Enumerate all the atoms in the system
    atoms = Iterators.product([1:size_x, 1:size_y, 1:size_z]...) |> collect |> vec
    # Enumerate all the coordinates in the system
    coords = [(a..., c) for c in [1, 2, 3], a in atoms] |> vec

    # Assign each coordinate an ordinal index
    index = Dict(zip(coords, 1:length(coords)))

    coupling_elements = vcat(
        [
            [
                (
                    # For each atom, take the cartesian coordinate from the couping term and get the corresponding index
                    get(index, (a..., c.d1), nothing),
                    # Using the displacement from the coupling, 
                    # get the atom to which "a" couples, as well as its cartesian coordinate (c.d2 field)
                    # To make the system periodic, let b = a + c.disp, the lattice coordinate of the target atom
                    # Atoms are enumerated from 1 to size_{x/y/z}. For periodicity, we use modulo, requiring
                    # enumeration from 0 to size_{x/y/z} - 1. Hence, we subtract (1,1,1) from b before applying
                    # the modulo and then add (1,1,1) after.
                    get(
                        index,
                        (
                            (1, 1, 1) .+
                            mod.(a .+ c.disp .- (1, 1, 1), (size_x, size_y, size_z))...,
                            c.d2,
                        ),
                        nothing,
                    ),
                    # Set the coupling between the two indices to the magnitude from c
                    c.k,
                ) for c in coupling
            ] for a in atoms
        ]...,
    )
    row = [c[1] for c in coupling_elements]
    col = [c[2] for c in coupling_elements]
    val = [c[3] for c in coupling_elements]
    # Assemble to coupling tuples into a sparse matrix
    V = sparse(row, col, val)
    return (index, V)
end

function derivatives(m, M, a, index, V, atoms, U, r, rdot, R, Rdot)
    # Intrinsic framework acceleration
    F_r = -V * r

    # Get the indices of the coordinates of each atom in 'atoms' interacting with the mobile particle
    atom_indices = [
        [
            get(index, c, ErrorException("Coordinate not found")) for
            c in [(a..., d) for d in [1, 2, 3]]
        ] for a in atoms
    ]

    # Get the equilibrium positions of each atom in 'atoms'. Atom (1,1,1) is at the origin
    atom_positions = [r[atom_indices[n]] .+ (atoms[n] .- 1) .* a for n in eachindex(atoms)]

    # Calculate the negative gradient of the interaction with each atom in 'atoms'
    F_r_particle = [-1 .* ForwardDiff.gradient(U, pos - R) for pos in atom_positions]
    # Force on mobile particles due to interaction using Newton's 3d law
    F_R = -sum(F_r_particle)
    view(F_r, vcat(atom_indices...)) .+= vcat(F_r_particle...)

    return (rdot, F_r ./ m, Rdot, F_R ./ M)
end

# 5th order Runge-Kutta step. 
# current_state is (r, rdot, R, Rdot)
# sys is (index, V)
function RKstep(m, M, a, sys, atoms, U, current_state, δt)
    k1 = derivatives(m, M, a, sys..., atoms, U, current_state...)
    k2 = derivatives(m, M, a, sys..., atoms, U, (current_state .+ k1 .* (δt / 4))...)
    k3 =
        derivatives(m, M, a, sys..., atoms, U, (current_state .+ (k1 .+ k2) .* (δt / 8))...)
    k4 = derivatives(
        m,
        M,
        a,
        sys...,
        atoms,
        U,
        (current_state .+ k3 .* δt .- k2 .* (δt / 2))...,
    )
    k5 = derivatives(
        m,
        M,
        a,
        sys...,
        atoms,
        U,
        (current_state .+ k1 .* (δt * 3 / 16) .+ k4 .* (δt * 9 / 16))...,
    )
    k6 = derivatives(
        m,
        M,
        a,
        sys...,
        atoms,
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
