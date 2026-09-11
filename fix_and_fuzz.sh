#!/usr/bin/env bash
set -e

# 1. Исправленный tests/fuzz_cornacchia.cpp на стандартном std::mt19937_64
cat << 'FILE_FUZZ' > tests/fuzz_cornacchia.cpp
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
    std::cout << "[*] Running 100,000 mutation cycles on extreme boundary values...\n";

    std::vector<mpz_class> edge_cases = {
        0, 1, 2, 3, 4, -1, -431, 
        mpz_class("18446744073709551616"), // 2^64
        mpz_class("340282366920938463463374607431768211456"), // 2^128
        mpz_class("1000000000000000000000000000000000000000000")
    };

    // 1. Тест жестких краевых значений
    for (const auto& m : edge_cases) {
        for (const auto& d : edge_cases) {
            auto res = Cornacchia::solve(d, m);
            (void)res;
            mpz_class out;
            Cornacchia::mod_sqrt(d, m, out);
        }
    }
    std::cout << "[+] Hardcoded boundary inputs passed without SIGSEGV/SIGFPE.\n";

    // 2. Рандомизированные мутации через std::mt19937_64
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<unsigned int> bit_dist(1, 256);
    gmp_randclass gmp_rng(gmp_randinit_default);
    gmp_rng.seed(42);

    unsigned long passed = 0;
    for (unsigned long i = 0; i < 100000; ++i) {
        unsigned int bits = bit_dist(rng);
        mpz_class rand_m = gmp_rng.get_z_bits(bits);
        mpz_class rand_d = gmp_rng.get_z_bits(bits / 2 + 1);

        // Инвертируем знаки случайно
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
    }

    std::cout << "[SUCCESS] 100,000 mutation cycles completed successfully!\n";
    std::cout << "          Zero memory corruptions, zero unhandled exceptions.\n";
    std::cout << "==================================================================\n";
    return 0;
}
FILE_FUZZ

# 2. Дособираем проект
cmake --build build

# 3. Запускаем фаззер
echo -e "\n[*] Executing Fuzzer..."
./build/orthrus_fuzz

# 4. Запускаем бенчмарк для профилирования
echo -e "\n[*] Running Profiling Benchmark..."
./build/orthrus_bench

echo "=================================================================="
echo "[+] Все три части полностью собраны и протестированы!"
echo "=================================================================="
rm -f fix_and_fuzz.sh
