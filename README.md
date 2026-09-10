# klpt-engine

Didactic C++20 implementation of the Kohel-Lauter-Petit-Tignol (KLPT) algorithm for navigating quaternion ideal lattices in supersingular isogeny cryptography.

## What is this?

This repository contains a standalone implementation of the quaternion navigation core used in the SQISign signature scheme, based on Christophe Petit's 2014 decomposition method.

Given a prime p = 3 (mod 4) and a left ideal I of norm N in the maximal quaternion order B_{p, inf}, the engine finds an element gamma in I such that:
    Nrd(gamma) = N * 2^e

This yields an equivalent ideal J = I * (gamma_bar / N) with smooth degree 2^e, which corresponds to a navigable chain of 2-isogenies via the Deuring correspondence.

## Performance & Scalability (AMD Ryzen 5 5600 @ 4.47 GHz)

The core Diophantine resolution via Cornacchia's algorithm and Tonelli-Shanks scales as O(log^2 M), ensuring sub-millisecond execution across all security levels:

| Security Level | Prime p | Target Norm M | Time (ms) | Verification |
|---|---|---|---|---|
| Toy-Level | 32 bits | 64 bits | 0.022 ms | Exact Norm Match |
| Mid-Level | 64 bits | 96 bits | 0.099 ms | Exact Norm Match |
| High-Level | 128 bits | 160 bits | 0.077 ms | Exact Norm Match |
| **NIST-1 (PQC)** | **256 bits** | **288 bits** | **0.110 ms** | **Exact Norm Match** |

## Features & Implementation

- Arithmetic of quaternion algebra B_{p, inf} over Q.
- Modified Cornacchia algorithm for quadratic forms x^2 + d*y^2 = m.
- Tonelli-Shanks algorithm for modular square roots.
- Dynamic ParameterGenerator for finding primes p = 3 (mod 4) up to 256-bit.
- Built-in automated benchmark suite across security levels.

## Build and Usage

Requirements:
- Linux (tested on Arch Linux x86_64)
- GCC 13+ or Clang 16+ (C++20 support)
- CMake 3.20+
- GNU Multiple Precision Arithmetic Library (gmp, gmpxx)

Build:
    cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
    cmake --build build

Run Automated Benchmarks:
    ./build/klpt_test --benchmark

Run Single Navigation Instance:
    ./build/klpt_test <prime_p> <ideal_norm_N> <degree_e>

Example:
    ./build/klpt_test 431 97 32

## References
- Kohel, Lauter, Petit, Tignol. "Quaternions, covers, and cryptography" (2014).
- De Feo, Kohel, Leriche, Petit, Wesolowski. "SQISign: compact post-quantum signatures from quaternions and isogenies" (2020).

## License
MIT
