include("../src/main.jl")
include("../src/plotting.jl")

## PARAMETERS
a = 3                       # Lattice constant in Å
k1 = 520                    # Nearest neighbor force constant in meV / Å²
k2 = 170                    # Next-nearest neighbor force constant in meV / Å²
m = 3.5                     # Lattice mass in meV * (ps / Å)²
ħ = 0.6582119569            # Planck constant in meV * ps

size_x = size_y = size_z = 40

# LATTICE
couplings = cubic_lattice(k1, k2)
sys = system(size_x, size_y, size_z, couplings)

R_mid = [(size_x + 1) / 2 - 1, (size_y + 1) / 2 - 1, (size_z + 1) / 2 - 1] .* a
R_edge = [(size_x + 0) / 2 - 1, (size_y + 0) / 2 - 1, (size_z + 1) / 2 - 1] .* a

λ = 3 / 4
# λ = 1 / sqrt(3)
U0 = 2000
# U(1.5)*2
m_unrelaxed = relaxation(sys, R_mid, λ, U0, 0)
e_unrelaxed = relaxation(sys, R_edge, λ, U0, 0)

diff_unrelaxed = sum(e_unrelaxed[1:2]) - sum(m_unrelaxed[1:2])
m_unrelaxed[4] |> maximum

m_relaxed = relaxation(sys, R_mid, λ, U0, 30)
e_relaxed = relaxation(sys, R_edge, λ, U0, 30)

diff_relaxed = sum(e_relaxed[1:2]) - sum(m_relaxed[1:2])


relaxed = relaxation(sys, R, λ, U0, 30)


atoms = Iterators.product([1:size_x, 1:size_y, 1:size_z]...) |> collect |> vec
atom_indices = [
    [
        get(sys[1], c, ErrorException("Coordinate not found")) for
        c in [(a..., d) for d in [1, 2, 3]]
    ] for a in atoms
]
@inline function U(r)
    res = U0 * exp(-dot(r, r) / 2 / λ^2)
    return res
end

r = zeros(length(sys[1]))
r_dot = zeros(length(sys[1]))

atom_positions = [r[atom_indices[n]] .+ (atoms[n] .- 1) .* a for n in eachindex(atoms)]
interaction = [U(pos - R) for pos in atom_positions] |> sum
interaction = [U(pos - R_mid) for pos in atom_positions] |> sum
U([60, 60, 60] - R_edge)
U([1.5, 1.5, 0])

exp(-dot([1.5, 1.5, 0], [1.5, 1.5, 0]) / 2 / λ^2)
dot([1.5, 1.5, 0], [1.5, 1.5, 0])
[60, 60, 60] - R_edge
R

@inline function U(r)
    res = U0 * exp(-dot(r, r) / 2 / λ^2)
    return res
end
