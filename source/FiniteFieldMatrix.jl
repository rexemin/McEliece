#=
Module with functions for matrices in finite fields.

Author: Ivan A. Moreno Soto
Last updated: 2019/July/15
=#

module FiniteFieldMatrix

using Random
using Nemo
using LinearAlgebra

include("./ExtensionField.jl")
using .ExtensionField

export
    multiplyPolyMatrices!,
    gaussianEliminationColumnPivoting!,
    GeneratorInverse!,
    getPermutation,
    getScrambler

"""
    multiplyPolyMatrices!(A, B, g)

Multiplies two matrices with coefficients in GF(2^m).
"""
function multiplyPolyMatrices!(A, B, g)
    C = zeros(UInt32, (size(A)[1], size(B)[2]))

    for i in 1:size(A)[1]
        for j in 1:size(A)[2]
            for k in 1:size(B)[2]
                temp = multiply(A[i, j], B[j, k])
                C[i, k] = add(C[i, k], temp)
            end
        end
    end

    # Applying modulo to every entry.
    C = modulo.(C, g)

    return C
end

"""
    searchPivot!(H, i, n)

Searches for the ith pivot in H in the n-i columns left.
"""
function searchPivot!(H, i, n)
    pivot = i
    for j in (i+1):n
        if H[i, j] == 1
            pivot = j
            break
        end
    end

    return pivot
end

"""
    permuteColumns!(H, i, pivot, k)

Permutes columns i and pivot with k entries each one in H.
"""
function permuteColumns!(H, i, pivot, k)
    for j in 1:k
        temp = H[j, i]
        H[j, i] = H[j, pivot]
        H[j, pivot] = temp
    end
end

"""
    pivotColumns!(H, i, n)

Pivots column i of matrix H. The function can search up to n-i columns
right.
"""
function pivotColumns!(H, i, n)
    k = size(H)[1]
    pivot = searchPivot!(H, i, n)

    if i == pivot
        error("Pivoting failed.")
    end

    permuteColumns!(H, i, pivot, k)

    return pivot
end

"""
    swapCodeSupport!(codeSupport, i, j)

Swaps two elements of a code support.
"""
function swapCodeSupport!(codeSupport, i, j)
    t = codeSupport[i]
    codeSupport[i] = codeSupport[j]
    codeSupport[j] = t
end

"""
    gaussianEliminationColumnPivoting!(H, n, codeSupport)

Applies Gauss-Jordan elimination to matrix H of dimension nxk. It also permutes
a given code support to keep a coherent goppa code.
"""
function gaussianEliminationColumnPivoting!(H, n, codeSupport)
    k = size(H)[2]
    # In GF(2)
    for i in 1:n
        # Pivoting the matrix to ensure that it can be reduced.
        if H[i, i] != 1
            pivot = pivotColumns!(H, i, k)
            swapCodeSupport!(codeSupport, i, pivot)
        end

        if H[i, i] != 1
            error("The matrix is not full rank.")
        end

        for j in (i+1):n
            if H[j, i] != 0
                for l in i:k
                    H[j, l] = (H[j, l] + H[i, l]) % 2
                end
            end
        end
    end

    for i in 2:n
        for j in 1:(i-1)
            if H[j, i] != 0
                for l in i:k
                    H[j, l] = (H[j, l] + H[i, l]) % 2
                end
            end
        end
    end
end

"""
    GeneratorInverse!(G)

Finds the right inverse of a generator matrix.
"""
function GeneratorInverse!(G)
    A = (G * transpose(G)) .% 2
    AInverse = BinaryInverse!(A)
    return (transpose(G) * AInverse) .% 2
end

"""
    BinaryInverse!(A)

Finds the inverse of a binary matrix A. Can throw an exception
if the inverse doesn't exist.
"""
function BinaryInverse!(A)
    B = Array{Int}(A)
    n = size(B)[1]
    R = ResidueRing(ZZ, 2)
    S = MatrixSpace(R, n, n)
    NemoInverse = inv(S(B))
    Inverse = zeros(UInt32, (n, n))
    for i in 1:n
        for j in 1:n
            if NemoInverse[i, j] == 1
                Inverse[i, j] = 1
            end
        end
    end
    return Inverse
end


"""
    getPermutation(n)

Returns a random permutation matrix of dimension nxn.
"""
function getPermutation(n)
    P = Matrix{UInt32}(I, n, n)
    return P[randperm(n), :]
end

"""
    fillBinaryMatrixRandomly!(M)

Fills randomly the entries of a matrix M with 0s and 1s.
"""
function fillBinaryMatrixRandomly!(M, k)
    for i in 1:k
        M[i, i] = 1
    end

    for i in 1:k
        for j in i:k
            if rand() <= 0.5
                M[i, j] = 1
            end
        end
    end

    for i in 2:k
        for j in 1:k
            M[i, j] = (M[i, j] + M[1, j]) % 2
        end
    end
end

"""
    getScrambler(k)

Returns a random binary non-singular matrix of dimension kxk.
"""
function getScrambler(k)
    S = SInverse = zeros(UInt32, (k, k))
    fillBinaryMatrixRandomly!(S, k)
    SInverse = BinaryInverse!(S)

    return S, SInverse
end

end # FiniteFieldMatrix
