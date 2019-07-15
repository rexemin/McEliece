#=
File with utility functions that don't really fit anywhere else. May change
locations/name in the future.

Author: Ivan Alejandro Moreno Soto
Last updated: 2019/July/15
=#

module Misc

export InvalidCodewordError,
    string2bits,
    bits2string

struct InvalidCodewordError <: Exception end

"""
    string2bits(plaintext, k)

Takes a plaintext and returns its binary representation with as much padding as
needed to make its length divisible by k.
"""
function string2bits(plaintext, k)
    plaintextBytes = Vector{UInt8}(plaintext)
    bits = Array{UInt32, 1}([digit for number in digits.(plaintextBytes, base=2, pad=8) for digit in number])

    while length(bits) % k != 0
        push!(bits, false)
    end

    return bits
end

"""
    bits2string(msg)

Returns the string that msg represents.
"""
function bits2string(msg)
    decomposed_msg = []
    # Taking 8 bits segments.
    for i in 1:div(length(msg), 8)
        current_character = msg[(8*(i-1) + 1):(8*i)]

        # Stopping if padding is found.
        if !(1 in current_character)
            break
        end

        current_string = ""
        for c in current_character
            current_string = current_string * string(c)
        end
        current_string = reverse(current_string)
        push!(decomposed_msg, Char(parse(Int, current_string, base = 2)))
    end

    # Concatenating every character in a single string.
    text = ""
    for c in decomposed_msg
        text = text * c
    end
    return text
end

end
