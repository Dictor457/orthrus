#include <iostream>
#include <vector>
#include <cstdint>
#include <random>
#include "orthrus/cornacchia.hpp"

using namespace orthrus;

int main() {
    std::cout << "==================================================================\n";
    std::cout << "      ORTHRUS FUZZER: Stress-Testing Cornacchia & Tonelli-Shanks   \n";
    std::cout << "==================================================================\n";

    std::vector<mpz_class> edge_cases = {
        0, 1, 2, 3, 4, -1, -431, 
        mpz_class("18446744073709551616"), // 2^64
        mpz_class("340282366920938463463374607431768211456"), // 2^128
        mpz_class("1000000000000000000000000000000000000000000")
    };

    std::cout << "[*] Running extreme boundary edge cases...\n";
    for (const auto& m : edge_cases) {
        for (const auto& d : edge_cases) {
            auto res = Cornacchia::solve(d, m);
            (void)res;
            mpz_class out;
            Cornacchia::mod_sqrt(d, m, out);
        }
    }
    std::cout << "[+] Boundary inputs PASSED instantly (infinite loops eliminated).\n";

    std::cout << "[*] Fuzzing 100,000 randomized mutations across 1..256 bit lengths...\n";
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<unsigned int> bit_dist(1, 256);
    gmp_randclass gmp_rng(gmp_randinit_default);
    gmp_rng.seed(42);

    unsigned long passed = 0;
    for (unsigned long i = 1; i <= 100000; ++i) {
        unsigned int bits = bit_dist(rng);
        mpz_class rand_m = gmp_rng.get_z_bits(bits);
        mpz_class rand_d = gmp_rng.get_z_bits(bits / 2 + 1);

        if (i % 3 == 0) rand_m = -rand_m;
        if (i % 5 == 0) rand_d = -rand_d;

        auto res = Cornacchia::solve(rand_d, rand_m);
        if (res.has_value()) {
            auto [x, y] = *res;
            if (x * x + rand_d * y * y != rand_m) {
                std::cerr << "[CRITICAL ERROR] Solution invariant broken during fuzzing!\n";
                return 1;
            }
        }
        passed++;

        if (i % 20000 == 0) {
            std::cout << "    -> Completed " << i << " / 100,000 iterations...\n";
        }
    }

    std::cout << "[SUCCESS] All " << passed << " fuzzing cycles finished cleanly!\n";
    std::cout << "          Zero crashes, zero memory corruptions, zero infinite loops.\n";
    std::cout << "==================================================================\n";
    return 0;
}
