#include <iostream>
#include <chrono>
#include <string>
#include "orthrus/orthrus.h"
#include "orthrus/quaternion.hpp"
#include "orthrus/order.hpp"
#include "orthrus/cornacchia.hpp"
#include "orthrus/ideal.hpp"
#include "orthrus/klpt.hpp"

using namespace orthrus;

int main(int argc, char* argv[]) {
    std::cout << "==================================================================\n";
    std::cout << "         ORTHRUS-CLI: Quaternion Navigation Engine v" << orthrus_version() << "\n";
    std::cout << "==================================================================\n";

    mpz_class p("431");
    mpz_class N("97");
    unsigned int e = 32;

    if (argc >= 2) p = mpz_class(argv[1]);
    if (argc >= 3) N = mpz_class(argv[2]);
    // Фикс варнинга -Wconversion через static_cast
    if (argc >= 4) e = static_cast<unsigned int>(std::stoul(argv[3]));

    std::cout << "[*] System Parameters:\n";
    std::cout << "    Base Characteristic (p) : " << p << " (p = 3 mod 4)\n";
    std::cout << "    Input Ideal Norm (N)    : " << N << "\n";
    std::cout << "    Smooth Target Degree    : 2^" << e << "\n";

    mpz_class neg_p = (-p % N + N) % N;
    mpz_class X0;
    if (!Cornacchia::mod_sqrt(neg_p, N, X0)) {
        std::cerr << "[ERROR] -p is not a QR mod N!\n";
        return 1;
    }

    IdealHNF I(N, X0);
    std::cout << "[+] 4x4 Hermite Normal Form Lattice constructed. det(H) = " << I.determinant() << "\n";

    std::cout << "[*] Executing KLPT 2D Lattice Sampling...\n";
    auto t_start = std::chrono::high_resolution_clock::now();
    auto sol = KLPTNavigator::find_smooth_generator(I, e, p);
    auto t_end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    if (!sol.has_value()) {
        std::cerr << "[FAIL] No solution found.\n";
        return 1;
    }

    Quaternion sc = sol->to_scaled_quaternion();
    std::cout << "[SUCCESS] Converged in " << ms << " ms!\n";
    std::cout << "    gamma = (" << sc.a0 << " + " << sc.a1 << "*i + " << sc.a2 << "*j + " << sc.a3 << "*k) / 2\n";
    std::cout << "    Norm: " << sol->norm(p) << " (Expected: " << N * (mpz_class(1) << e) << ")\n";
    std::cout << "==================================================================\n";
    return 0;
}
