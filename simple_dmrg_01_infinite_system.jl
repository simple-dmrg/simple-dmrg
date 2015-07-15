#!/usr/bin/env julia
#
# Simple DMRG tutorial.  This code contains a basic implementation of the
# infinite system algorithm
#
# Copyright 2013-2015 James R. Garrison and Ryan V. Mishmash.
# Open source under the MIT license.  Source code at
# <https://github.com/simple-dmrg/simple-dmrg/>

using Compat # allows compatibility with Julia versions 0.3 and 0.4

# Data structures to represent the block and enlarged block objects.
immutable Block
    length::Int
    basis_size::Int
    operator_dict::Dict{Symbol,AbstractMatrix{Float64}}
end

immutable EnlargedBlock
    length::Int
    basis_size::Int
    operator_dict::Dict{Symbol,AbstractMatrix{Float64}}
end

# For these objects to be valid, the basis size must match the dimension of
# each operator matrix.
isvalid(block::Union(Block,EnlargedBlock)) =
    all(op -> size(op) == (block.basis_size, block.basis_size), values(block.operator_dict))

# Model-specific code for the Heisenberg XXZ chain
model_d = 2  # single-site basis size

Sz1 = [0.5 0.0; 0.0 -0.5]  # single-site S^z
Sp1 = [0.0 1.0; 0.0 0.0]  # single-site S^+

H1 = [0.0 0.0; 0.0 0.0]  # single-site portion of H is zero

function H2(Sz1, Sp1, Sz2, Sp2)  # two-site part of H
    # Given the operators S^z and S^+ on two sites in different Hilbert spaces
    # (e.g. two blocks), returns a Kronecker product representing the
    # corresponding two-site term in the Hamiltonian that joins the two sites.
    const J = 1.0
    const Jz = 1.0
    return (J / 2) * (kron(Sp1, Sp2') + kron(Sp1', Sp2)) + Jz * kron(Sz1, Sz2)
end

# conn refers to the connection operator, that is, the operator on the edge of
# the block, on the interior of the chain.  We need to be able to represent S^z
# and S^+ on that site in the current basis in order to grow the chain.
initial_block = Block(1, model_d, @compat Dict{Symbol,AbstractMatrix{Float64}}(
    :H => H1,
    :conn_Sz => Sz1,
    :conn_Sp => Sp1,
))

function enlarge_block(block::Block)
    # This function enlarges the provided Block by a single site, returning an
    # EnlargedBlock.
    mblock = block.basis_size
    o = block.operator_dict

    # Create the new operators for the enlarged block.  Our basis becomes a
    # Kronecker product of the Block basis and the single-site basis.  NOTE:
    # `kron` uses the tensor product convention making blocks of the second
    # array scaled by the first.  As such, we adopt this convention for
    # Kronecker products throughout the code.
    enlarged_operator_dict = @compat Dict{Symbol,AbstractMatrix{Float64}}(
        :H => kron(o[:H], speye(model_d)) + kron(speye(mblock), H1) + H2(o[:conn_Sz], o[:conn_Sp], Sz1, Sp1),
        :conn_Sz => kron(speye(mblock), Sz1),
        :conn_Sp => kron(speye(mblock), Sp1),
    )

    return EnlargedBlock(block.length + 1,
                         block.basis_size * model_d,
                         enlarged_operator_dict)
end

function rotate_and_truncate(operator, transformation_matrix)
    # Transforms the operator to the new (possibly truncated) basis given by
    # `transformation_matrix`.
    return transformation_matrix' * (operator * transformation_matrix)
end

function single_dmrg_step(sys::Block, env::Block, m::Int)
    # Performs a single DMRG step using `sys` as the system and `env` as the
    # environment, keeping a maximum of `m` states in the new basis.

    @assert isvalid(sys)
    @assert isvalid(env)

    # Enlarge each block by a single site.
    sys_enl = enlarge_block(sys)
    if sys === env  # no need to recalculate a second time
        env_enl = sys_enl
    else
        env_enl = enlarge_block(env)
    end

    @assert isvalid(sys_enl)
    @assert isvalid(env_enl)

    # Construct the full superblock Hamiltonian.
    m_sys_enl = sys_enl.basis_size
    m_env_enl = env_enl.basis_size
    sys_enl_op = sys_enl.operator_dict
    env_enl_op = env_enl.operator_dict
    superblock_hamiltonian = kron(sys_enl_op[:H], speye(m_env_enl)) + kron(speye(m_sys_enl), env_enl_op[:H]) +
                             H2(sys_enl_op[:conn_Sz], sys_enl_op[:conn_Sp], env_enl_op[:conn_Sz], env_enl_op[:conn_Sp])

    # Call ARPACK to find the superblock ground state.  (:SR means find the
    # eigenvalue with the "smallest real" value.)
    #
    # But first, we explicitly modify the matrix so that it will be detected as
    # Hermitian by `eigs`.  (Without this step, the matrix is effectively
    # Hermitian but won't be detected as such due to small roundoff error.)
    superblock_hamiltonian = (superblock_hamiltonian + superblock_hamiltonian') / 2
    (energy,), psi0 = eigs(superblock_hamiltonian, nev=1, which=:SR)

    # Construct the reduced density matrix of the system by tracing out the
    # environment
    #
    # We want to make the (sys, env) indices correspond to (row, column) of a
    # matrix, respectively.  Since the environment (column) index updates most
    # quickly in our Kronecker product structure, psi0 is thus row-major.
    # However, Julia stores matrices in column-major format, so we first
    # construct our matrix in (env, sys) form and then take the transpose.
    psi0 = transpose(reshape(psi0, (env_enl.basis_size, sys_enl.basis_size)))
    rho = Hermitian(psi0 * psi0')

    # Diagonalize the reduced density matrix and sort the eigenvectors by
    # eigenvalue.
    fact = eigfact(rho)
    evals, evecs = fact[:values], fact[:vectors]
    permutation = sortperm(evals, rev=true)

    # Build the transformation matrix from the `m` overall most significant
    # eigenvectors.
    my_m = min(length(evals), m)
    indices = permutation[1:my_m]
    transformation_matrix = evecs[:, indices]

    truncation_error = 1 - sum(evals[indices])
    println("truncation error: ", truncation_error)

    # Rotate and truncate each operator.
    new_operator_dict = Dict{Symbol,AbstractMatrix{Float64}}()
    for (name, op) in sys_enl.operator_dict
        new_operator_dict[name] = rotate_and_truncate(op, transformation_matrix)
    end

    newblock = Block(sys_enl.length, my_m, new_operator_dict)

    return newblock, energy
end

function infinite_system_algorithm(L::Int, m::Int)
    block = initial_block
    # Repeatedly enlarge the system by performing a single DMRG step, using a
    # reflection of the current block as the environment.
    while 2 * block.length < L
        println("L = ", block.length * 2 + 2)
        block, energy = single_dmrg_step(block, block, m)
        println("E/L = ", energy / (block.length * 2))
    end
end

infinite_system_algorithm(100, 20)
