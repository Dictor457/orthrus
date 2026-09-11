#include <iostream>
#include <chrono>
#include "orthrus/orthrus.h"
#include "orthrus/quaternion.hpp"
#include "orthrus/order.hpp"
#include "orthrus/cornacchia.hpp"
#include "orthrus/ideal.hpp"
#include "orthrus/klpt.hpp"

using namespace orthrus;

int main() {
    std::cout << "[*] Starting Orthrus Profiling Benchmark (1,000 KLPT iterations)...\n";
    mpz_class p("431");
    mpz_class N("97");
    unsigned int e = 32;

    mpz_class neg_p = (-p % N + N) % N;
    mpz_class X0;
    Cornacchia::mod_sqrt(neg_p, N, X0);
    IdealHNF I(N, X0);

    auto t_start = std::chrono::high_resolution_clock::now();
    unsigned long solved = 0;

    for (int iter = 0; iter < 1000; ++iter) {
        auto sol = KLPTNavigator::find_smooth_generator(I, e, p, 10000);
        if (sol.has_value()) solved++;
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    std::cout << "[+] Completed 1,000 iterations in " << ms << " ms (" << ms / 1000.0 << " ms/op)\n";
    std::cout << "[+] Solved: " << solved << " / 1000\n";
    return 0;
}
