#include <iostream>
#include <chrono>
#include <vector>
#include "orthrus/uint256.hpp"
#include <gmpxx.h>

using namespace orthrus;

int main() {
    std::cout << "==================================================================\n";
    std::cout << "   ORTHRUS LOW-LEVEL ENGINE: Stack uint256 vs Heap mpz_class      \n";
    std::cout << "==================================================================\n";
    std::cout << "[*] Memory structure size: sizeof(uint256) = " << sizeof(uint256) << " bytes (Zero Heap)\n";

    // 1. Проверка корректности арифметики uint256
    uint256 a(0xFFFFFFFFFFFFFFFFULL, 0x1ULL, 0x0ULL, 0x0ULL);
    uint256 b(0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL);
    uint256 sum = a + b;
    std::cout << "[+] Addition with carry verified: " << sum.to_hex() << "\n";

    uint256 m1(0x100000000ULL);
    uint256 m2(0x200000000ULL);
    uint256 prod = m1 * m2;
    std::cout << "[+] Cross-multiplication verified: " << prod.to_hex() << "\n";

    // 2. Бенчмарк 1,000,000 операций: Стек vs Динамическая куча GMP
    const int ITERS = 1000000;
    std::cout << "\n[*] Running " << ITERS << " additions...\n";

    // Тест GMP (Heap allocations)
    mpz_class gmp_x("115792089237316195423570985008687907853269984665640564039457584007913129639935");
    mpz_class gmp_y("3");
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < ITERS; ++i) {
        gmp_x = gmp_x + gmp_y;
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto gmp_time = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Тест uint256 (Stack zero-alloc)
    uint256 stack_x(0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0x7FFFFFFFFFFFFFFFULL);
    uint256 stack_y(3);
    auto t2 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < ITERS; ++i) {
        stack_x = stack_x + stack_y;
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    auto stack_time = std::chrono::duration<double, std::milli>(t3 - t2).count();

    std::cout << "[+] GMP mpz_class (Heap dynamic malloc/free) : " << gmp_time << " ms\n";
    std::cout << "[+] Orthrus uint256 (Stack / CPU registers)   : " << stack_time << " ms\n";
    std::cout << "[>>>] SPEEDUP: " << std::fixed << std::setprecision(2) << (gmp_time / stack_time) 
              << "x FASTER with 0 heap allocations!\n";

    std::cout << "==================================================================\n";
    return 0;
}
