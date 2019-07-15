# McEliece
A proof-of-concept implementation of a McEliece cryptosystem in Julia.

Made under supervision of Dr. Kirill Morozov and mentorship of his PhD student Franz Aguirre Farro, while doing a
summer internship at the University of North Texas, in Denton, Texas.

Almost everything is coded from scratch. The whole system is just 6 modules that are (mostly) self-contained. `Test.jl` is an auxiliar file to time the performance of key generation, encryption, and decryption. The system depends on the following non-default packages:
- DelimitedFiles
- StatsBase
- Primes
- Nemo

At this moment in time, this implementation is not optimized at all, so expect slow running times for everything except encryption.
