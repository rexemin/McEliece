# McEliece
A proof-of-concept implementation of the McEliece cryptosystem in Julia.

Made under supervision of Dr. Kirill Morozov and mentorship of his PhD student Franz Aguirre Farro, while doing a
summer internship at the University of North Texas, in Denton, Texas.

Almost everything is coded from scratch. The whole system is just 6 modules that are (mostly) self-contained. `Test.jl` is an auxiliar file to time key generation, encryption, and decryption. The system depends on the following non-default packages:
- DelimitedFiles
- StatsBase
- Nemo

At this moment in time, the implementation is not optimized at all, so expect slow running times for everything except encryption.

You can read about the McEliece cryptosystem in the PDF accompanying this repository. You can see the current mean times
for various operations in the next table:

| Parameters (n, k, m, t) | Key generation | Encryption | Decryption |
| :---------------------: | :------------: | :--------: | :--------: |
| (2048, 1707, 11, 31)    | 62s            | 0.57ms     | 1.4s       |
| (2048, 1487, 11, 51)    | 109s           | 0.45ms     | 2.9s       |
| (2048, 1157, 11, 81)    | 171s           | 0.32ms     | 5.9s       |
| (4096, 3604, 12, 41)    | 312s           | 3.21ms     | 4.8s       |
| (4096, 3352, 12, 62)    | 601s           | 3.26ms     | 8.2s       |
| (4096, 2884, 12, 101)   | 552s           | 2.84ms     | 15.7s      |
| (6960, 5413, 13, 119)   | 1550s          | 19.60ms    | 36.5s      |
