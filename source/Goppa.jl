#=
Functions that automatically create Goppa codes.

Author: Ivan A. Moreno Soto
Last updated: 2019/July/15
=#

module GoppaCode

export
    generateGoppaCode,
    decode!

include("./Polynomials.jl")
using .Polynomials

include("./ExtensionField.jl")
using .ExtensionField

include("./FiniteFieldMatrix.jl")
using .FiniteFieldMatrix

include("./Misc.jl")
using .Misc

using StatsBase
using LinearAlgebra

using DelimitedFiles

"""
    selectCodeSupport(n::Int, m::Int)

Randomly selects and returns n distinct elements from the
extension field 2^m.
"""
function selectCodeSupport(n::Int, m::Int)
    return sample(Array{UInt32, 1}(0:(2^m - 1)), n, replace=false)
end

"""
    getY!(codeSupport, t, g)

Creates the Y matrix needed to compute the parity check matrix H. Basically
computes a Vandermonde matrix with the code support.
"""
function getY!(codeSupport, t, g)
    n = length(codeSupport)
    Y = zeros(UInt32, (t, n))

    for i in 1:n
        Y[1, i] = 1
    end

    for i in 1:n
        for j in 2:t
            Y[j, i] = power(codeSupport[i], (j-1), g)
        end
    end

    return Y
end

"""
    getZ!(goppaPolynomial, moduloPolynomial, codeSupport)

Creates the Z matrix needed to compute the parity check matrix H.
"""
function getZ!(goppaPolynomial, moduloPolynomial, codeSupport)
    n = length(codeSupport)
    Z = zeros(UInt32, (n, n))

    for i in 1:n
        Z[i, i] = inverse(evaluate!(codeSupport[i], goppaPolynomial, moduloPolynomial), moduloPolynomial)
    end

    return Z
end

"""
    toBits!(matrix, t, n, m)

Takes a matrix of dimension txn and returns a binary matrix
with m-bits representations of every element spread vertically.
"""
function toBits!(matrix, t, n, m)
    bits = falses(m*t, n)

    for i in 1:t
        for j in 1:n
            bits[((i - 1) * m + 1):(i * m), j] = digits(matrix[i, j], base=2, pad=m)
        end
    end

    return bits
end

"""
    generateGoppaCode(n::Int, t::Int, m::Int)

Generates a (n, k) irreducible linear Goppa code with minimum
distance d >= 2t + 1. Returns the generator matrix G, the goppa
polynomial, the code support, and the number of rows k of G.
"""
function generateGoppaCode(n::Int, t::Int, m::Int)
    # Get polynomials
    moduloPolynomial = primitivePolynomials[m]
    goppaPolynomial = randomPoly(t, m)
    # Get code support
    codeSupport = selectCodeSupport(n, m)

    # Compute H = YZ
    Y = getY!(codeSupport, t, moduloPolynomial)
    Z = getZ!(goppaPolynomial, moduloPolynomial, codeSupport)
    T = multiplyPolyMatrices!(Y, Z, moduloPolynomial)

    H = toBits!(T, t, n, m)
    gaussianEliminationColumnPivoting!(H, m*t, codeSupport)

    # Compute G
    k = n - m*t
    ATransposed = transpose(H[:, (m*t + 1):end])
    G = hcat(ATransposed, Matrix{UInt32}(I, k, k))

    return G, goppaPolynomial, codeSupport, k
end

"""
    syndrome!(msg, codeSupport, goppaPolynomial, m)

Computes the syndrome of a vector using a binary irreducible goppa code.
"""
function syndrome!(msg, codeSupport, goppaPolynomial, m)
    syndrome = Array{UInt32, 1}([0])
    for (c, γ) in zip(msg, codeSupport)
        if c == 1
            currentInverse = inverseRing!(Array{UInt32, 1}([γ, 1]), goppaPolynomial, m)
            syndrome = addRing!(syndrome, currentInverse)
        end
    end
    return syndrome
end

"""
    syndromeSquareRoot!(syndrome, goppaPolynomial, m)

Computes the square root of syndrome + x mod goppaPolynomial.
"""
function syndromeSquareRoot!(syndrome, goppaPolynomial, m)
    # Computing T(x) + x.
    syndromeInverse = inverseRing!(syndrome, goppaPolynomial, m)
    h = addRing!(syndromeInverse, PolynomialConstants['x'])

    # Computing the square root of x mod goppaPolynomial.
    goppaEven, goppaOdd = splitPolynomial!(goppaPolynomial, m)
    sx = multiplyRing!(goppaEven, inverseRing!(goppaOdd, goppaPolynomial, m), m)

    # Computing the odd and even parts of the square root of T(x) + x.
    evenh = [multiplyRing!(power(coeff, 2^(m-1), primitivePolynomials[m]), monomial(i-1), m) for (i, coeff) in enumerate(h[1:2:end])]
    oddh = [multiplyRing!(multiplyRing!([power(coeff, 2^(m-1), primitivePolynomials[m])], monomial(i-1), m), sx, m) for (i, coeff) in enumerate(h[2:2:end])]
    if length(oddh) == 0
        oddh = [copy(PolynomialConstants['0'])]
    end

    # Computing the square root by adding the even and odd parts.
    τ = copy(PolynomialConstants['0'])
    for term in vcat(evenh, oddh)
        τ = addRing!(τ, term)
    end

    return τ
end

"""
    halfwayXGCD!(goppaPolynomial, tau, t, m)

Runs the extended eucliedean algorithm halfway to compute the least
degree polynomials that can define the error locator polynomial for
a codeword with errors.
"""
function halfwayXGCD!(goppaPolynomial, tau, t, m)
    r = copy(goppaPolynomial)
    R = copy(tau)
    s = Array{UInt32, 1}([1])
    S = Array{UInt32, 1}([0])
    u = Array{UInt32, 1}([0])
    U = Array{UInt32, 1}([1])

    while length(R) - 1 > div(t, 2)
        q = divRing!(r, R, m)

        x = R; R = addRing!(r, multiplyRing!(q, R, m)); r = x
        x = S; S = addRing!(s, multiplyRing!(q, S, m)); s = x
        x = U; U = addRing!(u, multiplyRing!(q, U, m)); u = x
    end

    return R, U
end

"""
    findRootsErrorLocatorPolynomial!(sigma, codeSupport, m)

Finds the roots of sigma by evaluating every element of the code support.
"""
function findRootsErrorLocatorPolynomial!(sigma, codeSupport, m)
    totalErrors = 0
    errors = []
    for (i, gamma) in enumerate(codeSupport)
        if evaluate!(gamma, sigma, primitivePolynomials[m]) == 0
            push!(errors, i)
            totalErrors += 1
        end
    end

    return errors, totalErrors
end

"""
    decode!(msg, goppaPolynomial, codeSupport, m, t)

Given a msg with noise, finds a codeword in the goppa code defined
by the goppaPolynomial and codeSupport in GF(2^m).
"""
function decode!(msg, goppaPolynomial, codeSupport, m, t)
    # Compute syndrome
    msgSyndrome = syndrome!(msg, codeSupport, goppaPolynomial, m)
    if msgSyndrome == PolynomialConstants['0']
        throw(InvalidCodewordError())
    end

    tau = syndromeSquareRoot!(msgSyndrome, goppaPolynomial, m)

    # Using the extended euclidean algorithm to find the error locator polynomial.
    a, b = halfwayXGCD!(goppaPolynomial, tau, t, m)
    sigma = addRing!(multiplyRing!(a, a, m), multiplyRing!(PolynomialConstants['x'], multiplyRing!(b, b, m), m))
    sigma = makePolynomialMonic!(sigma, m)

    # Finding the roots of the error locator polynomial.
    errors, totalErrors = findRootsErrorLocatorPolynomial!(sigma, codeSupport, m)

    # Checking that it was a valid cyphertext.
    if totalErrors < t
        throw(InvalidCodewordError())
    end

    # Correcting the t errors.
    for position in errors
        msg[position] = (msg[position] + 1) % 2
    end

    return msg
end

end # GoppaCode
