struct Coupling
    disp::Tuple{Int64,Int64,Int64}  # Displacement
    d1::Int64                       # Coordinate of the first atom
    d2::Int64                       # Coordinate of the second atom
    k::Float64                      # Coupling strength
end

function cubic_lattice(k1, k2)
    # Nearest-neighbor and next-nearest neighbor displacements
    N_disp = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    NN_disp = [
        (1, 1, 0),
        (-1, 1, 0),
        (1, -1, 0),
        (-1, -1, 0),
        (0, 1, 1),
        (0, -1, 1),
        (0, 1, -1),
        (0, -1, -1),
        (1, 0, 1),
        (-1, 0, 1),
        (1, 0, -1),
        (-1, 0, -1),
    ]

    # Couplings for nearest and next-nearest atoms
    N_couplings = vcat(
        [
            [Coupling(v, n, m, -k1 * v[n] * v[m] / dot(v, v)) for n = 1:3, m = 1:3] |> vec
            for v in N_disp
        ]...,
    )
    N_couplings = filter(x -> x.k != 0, N_couplings)

    NN_couplings = vcat(
        [
            [Coupling(v, n, m, -k2 * v[n] * v[m] / dot(v, v)) for n = 1:3, m = 1:3] |> vec
            for v in NN_disp
        ]...,
    )

    NN_couplings = filter(x -> x.k != 0, NN_couplings)
    # The zero coupling terms are filtered out

    # Normalized vectors to nearest and next-nearest neighbors
    N_vec = [[x for x in v] ./ norm(v) for v in N_disp]
    NN_vec = [[x for x in v] ./ norm(v) for v in NN_disp]

    self_coupling_N = [v * v' for v in N_vec] |> sum
    self_coupling_NN = [v * v' for v in NN_vec] |> sum
    self_coupling_matrix = k1 .* self_coupling_N + k2 .* self_coupling_NN
    self_coupling =
        [Coupling((0, 0, 0), n, m, self_coupling_matrix[n, m]) for n = 1:3, m = 1:3] |> vec
    self_coupling = filter(x -> x.k != 0, self_coupling)
    couplings = vcat(N_couplings, NN_couplings, self_coupling)

    return couplings
end

# q here is q * a
function dynamical_matrix(m, coupling::Vector{Coupling})
    res = function f(q)
        D = zeros(ComplexF64, 3, 3)
        for c in coupling
            D[c.d1, c.d2] += c.k * exp(-1im * dot(q, c.disp)) ./ m
        end
        return D
    end

    return res

end

# Expanded dynamical matrix to q². The matrix has been divided by q².
function dynamical_matrix_small(a, m, coupling::Vector{Coupling})
    res = function f(θ, ϕ)
        D = zeros(ComplexF64, 3, 3)
        for c in coupling
            D[c.d1, c.d2] +=
                a^2 .* c.k *
                (-1im * dot([cos(ϕ) * sin(θ), sin(ϕ) * sin(θ), cos(θ)], c.disp))^2 /
                2 ./ m
        end
        return D
    end

    return res

end

# Friction matrix used in local-time formalism
function Loss_Matrix(a, m, DynamicalMatrix_Small)
    ρ = m / a^3
    function fun_int(q)
        # Solve the eigenproblem
        eig = eigen(DynamicalMatrix_Small(q...))
        vs = sqrt.(eig.values)
        ηs = eig.vectors
        res = sum([(ηs[:, jj] * ηs[:, jj]') ./ vs[jj]^3 for jj = 1:3])
        return res
    end
    res = hcubature(
        q -> sin(q[1]) .* fun_int(q) / 16 / π^2 / ρ,
        [0, 0],
        [π, 2 * π],
        rtol = 1e-4,
        # initdiv = 20,
    )
    return res
end

# Displacement correlation matrix for a given dynamical matrix, temperature, time, and separation between the atoms
function Corr_pos(ħ, m, dyn_mat, ΩT, t, disp)
    function fun_int(q)
        # Solve the eigenproblem
        eig = eigen(dyn_mat(q))
        Ωs = sqrt.(eig.values)
        ηs = eig.vectors

        res =
            [
                (ηs[:, ii] * ηs[:, ii]') / 2 / Ωs[ii] * coth(Ωs[ii] / 2 / ΩT) .*
                exp(-1im * Ωs[ii] * t) for ii = 1:3
            ] |> sum
        res = res .* exp(1im * dot(q, disp))
        return real(res)
    end
    res =
        hcubature(
            q -> fun_int(q) / (2 * π)^3,
            0.0 .* ones(3),
            2.0 .* π .* ones(3),
            rtol = 1e-3,
            # initdiv = 20,
            atol = 1e-6,
        ) ./ m .* ħ
    return res
end

# Velocity correlation matrix for a given dynamical matrix, temperature, time, and separation between the atoms
function Corr_speed(ħ, m, dyn_mat, ΩT, t, disp)
    function fun_int(q)
        # Solve the eigenproblem
        eig = eigen(dyn_mat(q))
        Ωs = sqrt.(eig.values)
        ηs = eig.vectors

        res =
            [
                (ηs[:, ii] * ηs[:, ii]') / 2 * Ωs[ii] * coth(Ωs[ii] / 2 / ΩT) .*
                exp(-1im * Ωs[ii] * t) for ii = 1:3
            ] |> sum
        res = res .* exp(1im * dot(q, disp))
        return real(res)
    end
    res =
        hcubature(
            q -> fun_int(q) / (2 * π)^3,
            0.0 .* ones(3),
            2.0 .* π .* ones(3),
            rtol = 1e-3,
            initdiv = 20,
            atol = 1e-6,
        ) ./ m .* ħ
    return res
end

# Recoil kernel 
function W(m, dyn_mat, t, disp)
    function fun_int(q)
        # Solve the eigenproblem
        eig = eigen(dyn_mat(q))
        Ωs = sqrt.(eig.values)
        ηs = eig.vectors

        res = [(ηs[:, ii] * ηs[:, ii]') / Ωs[ii]^2 .* cos(Ωs[ii] * t) for ii = 1:3] |> sum
        res = res .* exp(1im * dot(q, disp))
        return real(res)
    end
    res =
        hcubature(
            q -> fun_int(q) / (2 * π)^3,
            0.0 .* ones(3),
            2.0 .* π .* ones(3),
            rtol = 1e-2,
            # initdiv = 30,
            # atol = 1e-6,
        ) ./ m
    return res
end

# Homogeneous properties
# Mode amplitude
function Amplitude(ħ, Ωq, ΩT)
    η = 1e-12
    # Subtract a small number from p. The reason is that for low ΩT, p ≈ 1,
    # causing issues with the rand() generator
    n = rand(Geometric(1 - exp(-Ωq / ΩT) - η))
    res = √(n + 1 / 2) * √(2 * ħ / Ωq)
    return res
end

function homogeneous_init(ħ, m, dyn_mat, ΩT, size_x, size_y, size_z)

    atoms = Iterators.product([1:size_x, 1:size_y, 1:size_z]...) |> collect |> vec
    N = length(atoms)

    r_init = zeros(Float64, 3 * N)
    r_dot_init = zeros(Float64, 3 * N)

    momenta =
        Iterators.product(
            [
                2 * π .* (0:(size_x-1)) ./ size_x,
                2 * π .* (0:(size_y-1)) ./ size_y,
                2 * π .* (0:(size_z-1)) ./ size_z,
            ]...,
        ) |>
        collect |>
        vec
    # Drop the (0,0,0) momentum
    momenta = momenta[2:end]

    for q in momenta
        # @showprogress for q in momenta
        # Get the eigenvalues and 3-component eigenvectors for each q
        eigs = dyn_mat(q) |> eigen
        ηs = eigs.vectors
        # Get the corresponding frequencies
        Ωs = sqrt.(eigs.values |> real)

        # Generate a thermal amplitude for each frequency and multiply by a random phase factor
        ζs = [Amplitude(ħ, Ωq, ΩT) * exp(-2im * π * rand()) for Ωq in Ωs]

        # Multiply each ζ by the corresponding eigenvector and sum
        disp = [ζs[ii] .* ηs[:, ii] for ii = 1:3] |> sum
        # Multiply each ζ by the corresponding eigenvector and -i * Ω to get the speed, then sum
        speed = [-1im * Ωs[ii] * ζs[ii] .* ηs[:, ii] for ii = 1:3] |> sum

        # Lattice phase vectors
        ε = exp.([1im * dot(q, a) for a in atoms]) ./ sqrt(N)
        # Multiply each element of ε by disp and aggregate
        view(r_init, :) .+= real(kron(ε, disp))
        # Multiply each element of ε by speed and aggregate
        view(r_dot_init, :) .+= real(kron(ε, speed))

    end

    return r_init ./ √(m), r_dot_init ./ √(m)
end
