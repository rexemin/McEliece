#=
Polynomial operations made over unsigned integers in extension fields GF(2^m).

Author: Ivan A. Moreno Soto
Last updated: 2019/July/08
=#

module Polynomials

export
    PolynomialConstants,
    monomial,
    evaluate!,
    addRing!,
    multiplyRing!,
    divRing!,
    modRing!,
    makePolynomialMonic!,
    inverseRing!,
    isPolynomialIrreducible!,
    randomPoly,
    splitPolynomial!

include("./ExtensionField.jl")
using .ExtensionField

PolynomialConstants = Dict(
    '0' => Array{UInt32, 1}([0]),
    '1' => Array{UInt32, 1}([1]),
    'x' => Array{UInt32, 1}([0, 1])
)

"""
    monomial(degree)

Returns a monomial of the given degree.
"""
function monomial(degree)
    monomial = Array{UInt32, 1}([0 for _ in 0:degree])
    monomial[end] = 1
    return monomial
end

"""
    evaluate!(x, p, g)

Evaluates p(x) modulo g(x).
"""
function evaluate!(x, p, g)
    y = UInt32(0)

    for i in length(p):-1:1
        y = add(modulo(multiply(x, y), g), p[i])
    end

    return y
end

"""
    reducePolynomialArray!(p)

Takes a polynomial p and returns an array without unnecessary zero
coefficients.
"""
function reducePolynomialArray!(p)
    if length(p) == 1
        return p
    elseif !any(x -> x > 0, p)
        return [p[1]]
    else
        return p[1:findlast(x -> x > 0, p)]
    end
end

"""
    addRing!(a, b)

Adds two polynomials with coefficients in GF(2^m).
"""
function addRing!(a, b)
    if length(a) >= length(b)
        longPolynomial = a
        shortPolynomial = b
    else
        longPolynomial = b
        shortPolynomial = a
    end

    total = Array{UInt32, 1}([coeff for coeff in longPolynomial])

    for i in 1:length(shortPolynomial)
        total[i] = add(shortPolynomial[i], longPolynomial[i])
    end

    return reducePolynomialArray!(total)
end

"""
    multiplyRing!(a, b, m)

Multiplies two polynomials with coefficients in GF(2^m).
"""
function multiplyRing!(a, b, m)
    product = Array{UInt32, 1}([0 for _ in 1:(length(a) + length(b) - 1)])

    for i in 1:length(a)
        for j in 1:length(b)
            term = modulo(multiply(a[i], b[j]), primitivePolynomials[m])
            product[i+j - 1] = add(product[i+j - 1], term)
        end
    end

    return reducePolynomialArray!(product)
end

"""
    polynomialLongDivision!(a, b, m)

Returns the quotient and the modulo of a/b with coefficients
in GF(2^m).
"""
function polynomialLongDivision!(a, b, m)
    if a == b
        return Array{UInt32, 1}([1]), Array{UInt32, 1}([0])
    elseif length(a) < length(b) #a < b
        return Array{UInt32, 1}([0]), a
    end

    quotient = Array{UInt32, 1}([])
    remainder = a
    for i in (length(a) - length(b)):-1:0
        if length(remainder) >= length(b) + i
            step = Array{UInt32, 1}([0 for _ in 0:i])
            step[end] = modulo(multiply(inverse(b[end], primitivePolynomials[m]), remainder[end]), primitivePolynomials[m])

            partial = multiplyRing!(b, step, m)

            remainder = addRing!(remainder, partial)
            push!(quotient, step[end])
        else
            push!(quotient, 0)
        end
    end

    reverse!(quotient) # push! adds the coefficients in reverse order.

    return reducePolynomialArray!(quotient), reducePolynomialArray!(remainder)
end

"""
    divRing!(a, b, m)

Returns the quotient of a/b with coefficients in GF(2^m).
"""
function divRing!(a, b, m)
    quotient, _ = polynomialLongDivision!(a, b, m)
    return quotient
end

"""
    modRing!(a, b, m)

Returns the remainder of a/b with coefficients in GF(2^m).
"""
function modRing!(a, b, m)
    _, remainder = polynomialLongDivision!(a, b, m)
    return remainder
end

"""
    makePolynomialMonic!(p, m)

Returns p as a monic polynomial by dividing p by the inverse of its leading
coefficient.
"""
function makePolynomialMonic!(p, m)
    return multiplyRing!([inverse(p[end], primitivePolynomials[m])], p, m)
end

"""
    ringgcs!(a, b, m)

Computes the GCD of a and b with coefficients in GF(2^m).
"""
function ringgcd!(a, b, m)
    while b != PolynomialConstants['0']
        t = copy(b)
        b = modRing!(a, b, m)
        a = t
    end

    return makePolynomialMonic!(a, m)
end

"""
    inverseRing!(a, g, m)

Computes the inverse of a mod g with coefficients in GF(2^m). Throws
an exception if the inverse doesn't exist.
"""
function inverseRing!(a, g, m)
    # Extended euclidean algorithm.
    r = copy(a)
    R = copy(g)
    s = Array{UInt32, 1}([1])
    S = Array{UInt32, 1}([0])
    t = Array{UInt32, 1}([0])
    T = Array{UInt32, 1}([1])

    while R != Array{UInt32, 1}([0])
        q = divRing!(r, R, m)

        x = R; R = addRing!(r, multiplyRing!(q, R, m)); r = x
        x = S; S = addRing!(s, multiplyRing!(q, S, m)); s = x
        x = T; T = addRing!(t, multiplyRing!(q, T, m)); t = x
    end

    c = inverse(r[end], primitivePolynomials[m])

    if makePolynomialMonic!(r, m) != PolynomialConstants['1']
        error("The element doesn't have an inverse.")
    end

    return multiplyRing!([c], s, m)
end

"""
    isPolynomialIrreducible!(polynomial, n, m)

Checks if the given polynomial has at least one root in GF(2^m). Returns
true if it doesn't have roots, false otherwise.
"""
function isPolynomialIrreducible!(polynomial, n, m)
    for i in 0:(2^m - 1)
        if evaluate!(UInt32(i), polynomial, primitivePolynomials[m]) == 0
            return false
        end
    end

    return true
end

"""
    getRandomCoefficients(t, m)

Returns a random sequence of elements from GF(2^m)
to serve as the coefficients of polynomial of degree t.
"""
function getRandomCoefficients(t, m)
    polynomial = [UInt32(0) for _ in 1:(t + 1)]
    polynomial[end] = UInt32(1) # Making the polynomial monic.
    # Adding random coefficients to the polynomial.
    for i in 1:t
        if rand() >= 0.5
            polynomial[i] = UInt32(rand(1:(2^m - 1)))
        end
    end

    return polynomial
end

"""
    randomPoly(t::Int, m::Int)

Returns a random monic irreducible polynomial of degree t
over the extension field GF(2^m).
"""
function randomPoly(t::Int, m::Int)
    polynomial = getRandomCoefficients(t, m)

    while !isPolynomialIrreducible!(polynomial, t, m)
        polynomial = getRandomCoefficients(t, m)
    end

    return polynomial
end

"""
    splitPolynomial!(p, m)

Splits the polynomial p into even and odd terms in GF(2^m).
"""
function splitPolynomial!(p, m)
    even = [squareroot(coeff, m) for coeff in p[1:2:end]]
    odd = [squareroot(coeff, m) for coeff in p[2:2:end]]
    return even, odd
end

end # Module Polynomials
