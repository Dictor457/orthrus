# klpt-engine

Didactic C++20 implementation of the Kohel-Lauter-Petit-Tignol (KLPT) algorithm for navigating quaternion ideal lattices in supersingular isogeny cryptography.

## What is this?

This repository contains a standalone implementation of the quaternion navigation core used in the SQISign signature scheme, based on Christophe Petit's 2014 decomposition method.

Given a prime p = 3 (mod 4) and a left ideal I of norm N in the maximal quaternion order B_{p, inf}, the engine finds an element gamma in I such that:
    Nrd(gamma) = N * 2^e

This yields an equivalent ideal J = I * (gamma_bar / N) with smooth degree 2^e, which corresponds to a navigable chain of 2-isogenies via the Deuring correspondence.

## Current Scope & Limitations

What is implemented:
- Arithmetic of quaternion algebra B_{p, inf} over Q.
- Tonelli-Shanks algorithm for modular square roots.
- Modified Cornacchia algorithm for quadratic forms x^2 + d*y^2 = m.
- Petit's decomposition gamma = C + D*j for ideal representation.

What is NOT implemented:
- Supersingular elliptic curve point arithmetic.
- Finite field arithmetic over F_{p^2}.
- Velu / sqrt(Velu) isogeny evaluation.
- Full SQISign signature protocol.

This is strictly a mathematical and algorithmic prototype demonstrating the quaternion lattice reduction step.

## Build and Usage

Requirements:
- Linux (tested on Arch Linux)
- GCC 13+ / Clang 16+
- CMake 3.20+
- GMP library (libgmp, gmpxx)

Build:
    cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
    cmake --build build

Run:
    ./build/klpt_test <prime_p> <ideal_norm_N> <degree_e>

Example:
    ./build/klpt_test 431 97 32

## Benchmark

Hardware: AMD Ryzen 5 5600 @ 4.47 GHz (Arch Linux)
- Target: p = 431, N = 97, e = 32 (Target norm = 416,611,827,712)
- Resolution time: < 1 ms

## References
- Kohel, Lauter, Petit, Tignol. "Quaternions, covers, and cryptography" (2014).
- De Feo, Kohel, Leriche, Petit, Wesolowski. "SQISign: compact post-quantum signatures from quaternions and isogenies" (2020).

## License
MIT
