#=
Companion file to the McEliece implementation to test the main functions of it.
It's not really a benchmark. Know that if you run the whole set of parameters,
it's going to take a long time, it isn't the most optimized implementation.

Author: Ivan Alejandro Moreno Soto
Last updated: 2019/July/15
=#

using DelimitedFiles

include("./McEliece.jl")
using .McEliece

"""
    testKeyGeneration(n, t, m, teststrings, publicDirectory, privateDirectory, saveKeysToFile=false)

Times key generation, encryption and decryption of the implementation. If you want, you can save the
keys to text files during the process.

n, t, and m are the parameters for the goppa code. teststrings is the list of strings to encrypt and
decrypt.
"""
function testKeyGeneration(n, t, m, teststrings, publicDirectory, privateDirectory, saveKeysToFile=false)
    publicKey, privateKey = generateKeys(n, t, m)

    if saveKeysToFile
        saveKeys(publicKey, privateKey, publicDirectory, privateDirectory)
    end

    for plaintext in teststrings
        println("Timing encryption of $(plaintext)")
        @time cyphertext = encrypt(plaintext, publicKey[1], t)

        println("Timing decryption of $(plaintext)")
        @time decryptedtext = decrypt(cyphertext,
            publicKey[1],
            privateKey[1],
            privateKey[2],
            privateKey[3],
            privateKey[4],
            m, t)

        println("decrypted == plaintext? $(decryptedtext == plaintext)")
    end
end

GPubLoc = "./keys/public_g.txt"
GPrivLoc = "./keys/private_g.txt"
tLoc = "./keys/public_t.txt"
mLoc = "./keys/public_m.txt"
SLoc = "./keys/private_s.txt"
PLoc = "./keys/private_p.txt"
supLoc = "./keys/private_support.txt"
goppaLoc = "./keys/private_polynomial.txt"

parameters = [
    [2048, 31, 11],
    # [2048, 51, 11],
    # [2048, 81, 11],
    # [4096, 41, 12],
    # [4096, 62, 12],
    # [4096, 101, 12],
]

directories = [
    [".\\keys\\2048-1\\public", ".\\keys\\2048-1\\private"],
    # [".\\keys\\2048-2\\public", ".\\keys\\2048-2\\private"],
    # [".\\keys\\2048-3\\public", ".\\keys\\2048-3\\private"],
    # [".\\keys\\4096-1\\public", ".\\keys\\4096-1\\private"],
    # [".\\keys\\4096-2\\public", ".\\keys\\4096-2\\private"],
    # [".\\keys\\4096-3\\public", ".\\keys\\4096-3\\private"],
]

teststrings = [
    "Hello!",
    "Where are we meeting?",
    "I was wondering when you would show up."
]

for ((n, t, m), d) in zip(parameters, directories)
    println("Timing generation of keys with n=$(n), t=$(t), m=$(m)")
    testKeyGeneration(n, t, m, teststrings, d[1], d[2])
end
