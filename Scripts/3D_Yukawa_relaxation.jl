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
R_edge = [(size_x + 0) / 2 - 1, (size_y + 1) / 2 - 1, (size_z + 1) / 2 - 1] .* a

λ = 1 / 2
U0 = 4000

m_unrelaxed = relaxation(sys, R_mid, λ, U0, 0)
e_unrelaxed = relaxation(sys, R_edge, λ, U0, 0)

diff_unrelaxed = sum(e_unrelaxed[1:2]) - sum(m_unrelaxed[1:2])
m_unrelaxed[4] |> maximum

m_relaxed = relaxation(sys, R_mid, λ, U0, 30)
e_relaxed = relaxation(sys, R_edge, λ, U0, 30)

diff_relaxed = sum(e_relaxed[1:2]) - sum(m_relaxed[1:2])
diff_relaxed = e_relaxed[1:2] .- m_relaxed[1:2]
