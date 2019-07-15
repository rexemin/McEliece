#=
This module contains functions that allow operations in extension fields
GF(2^m). It also contains a dictionary with primitive polynomials for
small extension fields.

Author: Ivan A. Moreno Soto
Last updated: 2019/July/02
=#

module ExtensionField

export primitivePolynomials,
    add,
    multiply,
    divide,
    modulo,
    power,
    inverse,
    squareroot,
    polygcd

primitivePolynomials = Dict(
    2 => UInt32(7),
    3 => UInt32(11),
    4 => UInt32(19),
    5 => UInt32(37),
    6 => UInt32(103),
    7 => UInt32(131),
    8 => UInt32(487),
    9 => UInt32(901),
    10 => UInt32(1033),
    11 => UInt32(2053),
    12 => UInt32(4179),
    13 => UInt32(8219),
    14 => UInt32(16707),
    15 => UInt32(32785),
)

"""
    add(a::UInt32, b::UInt32)

XORs two polynomials in GF(2).
"""
function add(a::UInt32, b::UInt32)
    return a ⊻ b
end

"""
    multiply(a::UInt32, b::UInt32)

Multiplies two polynomials in GF(2). It doesn't reduce modulo
an irreducible polynomial.
"""
function multiply(a::UInt32, b::UInt32)
    partial = Array{UInt32, 1}()
    for (i, bit) in enumerate(digits(b, base=2))
        if bit == 1
            push!(partial, a << (i - 1))
        end
    end

    c = UInt32(0)
    for i in partial
        c ⊻= i
    end

    return c
end

"""
    division(a::UInt32, b::UInt32)

Computes a/b. Returns the quotient and the remainder.
"""
function division(a::UInt32, b::UInt32)
    if a == b
        return UInt32(1), UInt32(0)
    elseif length(digits(a, base=2)) < length(digits(b, base=2)) #a < b
        return UInt32(0), a
    end

    quotient = UInt32(0)
    remainder = a
    for i in (length(digits(a, base=2)) - length(digits(b, base=2))):-1:0
        if length(digits(remainder, base=2)) >= length(digits(b << i, base=2))
            remainder ⊻= b << i
            quotient ⊻= UInt32(2^i)
        end
    end

    return quotient, remainder
end

"""
    divide(a::UInt32, b::UInt32)

Computes a/b but only returns the quotient.
"""
function divide(a::UInt32, b::UInt32)
    quotient, _ = division(a, b)
    return quotient
end

"""
    modulo(a::UInt32, b::UInt32)

Computes a/b but only returns the remainder.
"""
function modulo(a::UInt32, b::UInt32)
    _, remainder = division(a, b)
    return remainder
end

"""
    power(x::UInt32, n::Int, g::UInt32)

Computes and returns x^n modulo g in GF(2^m).
"""
function power(x::UInt32, n::Int, g::UInt32)
    if n == 0
        return UInt32(1)
    end

    y = UInt32(1)

    while n > 1
        if n % 2 == 0
            n = n / 2
        else
            y = modulo(multiply(x, y), g)
            n = (n - 1)/2
        end

        x = modulo(multiply(x, x), g)
    end

    return modulo(multiply(x, y), g)
end

"""
    inverse(a::UInt32, g::UInt32)

Computes and returns the inverse of a mod g(x).
"""
function inverse(a::UInt32, g::UInt32)
    # Extended euclidean algorithm.
    t = UInt32(0)
    T = UInt32(1)
    r = g
    R = a

    while R != 0
        q = divide(r, R)
        x = T; T = add(t, multiply(q, T)); t = x
        x = R; R = add(r, multiply(q, R)); r = x
    end

    if r > 1
        error("The element $(a) doesn't have an inverse mod $(g).")
    end

    t = modulo(t, g)

    return t
end

"""
    squareroot(p::UInt32, m)

Computes the square root of p in GF(2^m).
"""
function squareroot(p::UInt32, m)
    if p == 0
        return UInt32(0)
    elseif p == 1
        return UInt32(1)
    end

    return power(p, 2^(m - 1), primitivePolynomials[m])
end

"""
    polygcd(a::UInt32, b::UInt32)

Computes the greatest common divisor of polynomials a and b.
"""
function polygcd(a::UInt32, b::UInt32)
    while b != 0
        t = b
        b = modulo(a, b)
        a = t
    end

    return a
end

end # ExtensionField
