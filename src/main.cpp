#include <iostream>
#include <chrono>
#include "orthrus/quaternion.hpp"
#include "orthrus/order.hpp"
#include "orthrus/cornacchia.hpp"
#include "orthrus/ideal.hpp"
#include "orthrus/klpt.hpp"

using namespace orthrus;

int main(int argc, char* argv[]) {
    std::cout << "==================================================================\n";
    std::cout << "     ORTHRUS: High-Assurance KLPT Quaternion Navigation Core      \n";
    std::cout << "==================================================================\n";

    mpz_class p("431");
    mpz_class N("97");
    unsigned int e = 32;

    if (argc >= 2) p = mpz_class(argv[1]);
    if (argc >= 3) N = mpz_class(argv[2]);
    if (argc >= 4) e = std::stoul(argv[3]);

    std::cout << "[*] System Parameters:\n";
    std::cout << "    Base Characteristic (p) : " << p << " (p = 3 mod 4)\n";
    std::cout << "    Input Ideal Norm (N)    : " << N << "\n";
    std::cout << "    Smooth Target Degree    : 2^" << e << "\n";

    // 1. Построение идеала через HNF 4x4
    mpz_class neg_p = (-p % N + N) % N;
    mpz_class X0;
    if (!Cornacchia::mod_sqrt(neg_p, N, X0)) {
        std::cerr << "[ERROR] -p is not a QR mod N!\n";
        return 1;
    }

    IdealHNF I(N, X0);
    mpz_class det = I.determinant();
    std::cout << "\n[+] 4x4 Hermite Normal Form (HNF) Lattice constructed:\n";
    std::cout << "    det(H) = " << det << " (Expected N^2 = " << N * N << ")\n";

    if (det != N * N) {
        std::cerr << "[CRITICAL] HNF Determinant mismatch!\n";
        return 1;
    }

    // 2. Декомпозиция Пети
    std::cout << "\n[*] Running 2D Lattice Sampling & Retry Loop for Nrd(gamma) = N * 2^" << e << "...\n";

    auto t_start = std::chrono::high_resolution_clock::now();
    auto sol = KLPTNavigator::find_smooth_generator(I, e, p);
    auto t_end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    if (!sol.has_value()) {
        std::cerr << "[FAIL] Solution not found within sampling budget.\n";
        return 1;
    }

    OrderElement gamma = *sol;
    mpz_class target_norm = N * (mpz_class(1) << e);
    mpz_class actual_norm = gamma.norm(p);
    Quaternion scaled = gamma.to_scaled_quaternion();

    std::cout << "[SUCCESS] Lattice Sampling converged in " << ms << " ms!\n";
    std::cout << "    gamma = (" << scaled.a0 << " + " 
                               << scaled.a1 << "*i + " 
                               << scaled.a2 << "*j + " 
                               << scaled.a3 << "*k) / 2\n";
    std::cout << "    Computed Norm : " << actual_norm << "\n";
    std::cout << "    Target Norm   : " << target_norm << "\n";

    if (actual_norm == target_norm) {
        std::cout << "[VERIFIED] Exact norm match! Path of length " << e << " ready for curve evaluation.\n";
    } else {
        std::cerr << "[CRITICAL ERROR] Norm equation violated!\n";
        return 1;
    }

    std::cout << "==================================================================\n";
    return 0;
}
