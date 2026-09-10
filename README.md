# SQISign KLPT-Engine

A high-performance C++20 implementation of the **KLPT (Kohel-Lauter-Petit-Tignol)** algorithm, providing the quaternion navigation core for the **SQISign** post-quantum digital signature scheme.

## Overview

In isogeny-based cryptography, navigating the supersingular l-isogeny graph directly over finite fields takes exponential time O(sqrt(p)).

Via the **Deuring Correspondence**, supersingular elliptic curves correspond to maximal orders in the quaternion algebra B_{p, inf}, and isogenies correspond to integral ideals. The KLPT algorithm solves the navigation problem in polynomial time on the quaternion side:

Given an ideal I of arbitrary degree N, KLPT computes an equivalent ideal:
    J = I * (gamma_bar / N)
such that Nrd(J) = 2^e, reducing the graph navigation to an efficient evaluation of a degree-2^e isogeny chain on the curve.

## Mathematical Architecture

1. Quaternion Algebra B_{p, inf}:
   For p = 3 mod 4, generators <1, i, j, k> with:
   i^2 = -1, j^2 = -p, ij = -ji = k.
2. Cornacchia Algorithm & Tonelli-Shanks:
   Solving diophantine equations x^2 + d*y^2 = M over Z.
3. Petit Decomposition:
   Decomposing elements gamma = C + D*j in I with C = D*X0 (mod N) such that gamma * alpha_bar = 0 (mod N) holds identically.

## Build and Run

Requirements:
- Arch Linux / Modern Linux
- GCC 13+ or Clang 16+ (C++20 support)
- CMake 3.20+
- GNU Multiple Precision Arithmetic Library (gmp, gmpxx)

Build commands:
    cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
    cmake --build build

Run with custom parameters:
    ./build/klpt_test <prime_p> <ideal_norm_N> <degree_e>
Example:
    ./build/klpt_test 431 97 32

## Performance

Tested on AMD Ryzen 5 5600 running Arch Linux:
- RepresentInteger (128-bit norm): < 1 ms
- KLPT Path Resolution (2^32-isogeny): < 2 ms

## License
MIT
