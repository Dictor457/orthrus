# Orthrus

Experimental C++20 implementation of Petit's heuristic quaternion lattice sampling algorithm (the core Diophantine step of the KLPT algorithm in SQISign).

## Algorithmic Scope

This repository implements the heuristic lattice sampling algorithm proposed by Christophe Petit (2014) for resolving the norm equation in maximal orders of quaternion algebras:

Given a prime p = 3 (mod 4) and a cyclic left ideal I = O_0<X0 + j, N> of norm N in the maximal order O_0 of B_{p, inf}, the engine finds an element gamma in I such that:
    Nrd(gamma) = N * 2^e

This yields an equivalent ideal J = I * (gamma_bar / N) of smooth power-of-two norm 2^e.

## Limitations & Missing Components

To be strictly clear, this repository is NOT a complete SQISign implementation and NOT a full academic KLPT solver. Specifically:

1. No Elliptic Curves: Does not evaluate supersingular elliptic curve points or isogeny chains over F_{p^2} (no Velu / sqrt(Velu) formulas).
2. No Signature Scheme: Does not implement the SQISign identification or digital signature protocols.
3. Cyclic Ideals Only: The current 4x4 Hermite Normal Form (HNF) construction assumes cyclic ideals of the form O_0<X0 + j, N>.
4. No 4D LLL Reduction: General arbitrary ideals require lattice basis reduction (Lenstra-Lenstra-Lovasz) in dimension 4, which is not yet implemented.
5. No Eichler Orders: Does not compute arbitrary right orders O_R(I) = {q in B_{p, inf} | I*q subset I}.
6. No Ideal Filtration: Does not factor ideals of composite norm into chains of prime-degree ideals.
7. Heuristic Search: The resolution of the modular quadratic form relies on 2D lattice coordinate sampling and retry loops; it is not deterministic polynomial time for worst-case inputs.

## Implemented Components

- Quaternion arithmetic in B_{p, inf} (conjugation, reduced trace, reduced norm, multiplication).
- Maximal order O_0 = Z*1 + Z*i + Z*(1+j)/2 + Z*(i+k)/2 for p = 3 (mod 4).
- 4x4 Hermite Normal Form (HNF) lattice representation for cyclic left ideals.
- Tonelli-Shanks algorithm for modular square roots (prime moduli only).
- Modified Cornacchia algorithm for integer quadratic forms x^2 + d*y^2 = m.
- Stack-allocated 256-bit unsigned integer arithmetic (uint256) using x86-64 carry intrinsics (_addcarry_u64, _subborrow_u64).
- C-compatible ABI (orthrus.h) with string-based GMP serialization.
- Arch Linux PKGBUILD manifest.

## Measured Benchmarks

Hardware: AMD Ryzen 5 5600 @ 4.47 GHz, Arch Linux x86_64, GCC -O3:
- Single resolution (p = 431, N = 97, e = 32): ~0.45 ms
- 1,000 batch resolutions: 608 ms (~0.608 ms per operation, 100% convergence rate)
- Fuzzing suite: 100,000 mutation cycles over boundary values: 0 crashes / 0 memory leaks
- uint256 stack addition benchmark: 0.45 ms per 1,000,000 additions (vs 5.08 ms in mpz_class)

## Building

Dependencies:
- Linux (tested on Arch Linux)
- GCC 13+ or Clang 16+
- CMake 3.20+
- GNU Multiple Precision Arithmetic Library (libgmp, gmpxx)

Commands:
    cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
    cmake --build build

Binaries generated:
- build/liborthrus.so (Shared library)
- build/liborthrus.a (Static library)
- build/orthrus-cli (Command-line tool)
- build/c_demo (C-ABI test)
- build/orthrus_fuzz (Fuzz testing harness)
- build/orthrus_bench (Benchmark runner for Linux perf)
- build/test_uint256 (uint256 vs GMP benchmark)

Run CLI:
    ./build/orthrus-cli <prime_p> <ideal_norm_N> <degree_e>

Example:
    ./build/orthrus-cli 431 97 32

## References

- Christophe Petit. "Faster algorithms for isogeny problems using quaternions." (2014).
- David Kohel, Kristin Lauter, Christophe Petit, Jean-Pierre Tignol. "Quaternions, covers, and cryptography." (2014).
- Luca De Feo, David Kohel, Antonin Leriche, Christophe Petit, Benjamin Wesolowski. "SQISign: compact post-quantum signatures from quaternions and isogenies." (2020).

## License

MIT
