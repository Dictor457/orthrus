#!/usr/bin/env bash
set -e

mkdir -p include/orthrus src

# 1. CMakeLists.txt с именем проекта ORTHRUS и жесткими флагами контроля
cat << 'FILE_CMAKE' > CMakeLists.txt
cmake_minimum_required(VERSION 3.20)
project(orthrus CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# "Control your environment, or the defaults will control you."
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra -Wpedantic -Wshadow -Wconversion")

find_path(GMP_INCLUDE_DIR NAMES gmp.h gmpxx.h)
find_library(GMP_LIBRARY NAMES gmp)
find_library(GMPXX_LIBRARY NAMES gmpxx)

if (NOT GMP_LIBRARY OR NOT GMPXX_LIBRARY)
    message(FATAL_ERROR "GMP/GMPXX not found! Run: sudo pacman -S gmp")
endif()

include_directories(include ${GMP_INCLUDE_DIR})

add_executable(orthrus_test src/main.cpp)
target_link_libraries(orthrus_test PRIVATE ${GMPXX_LIBRARY} ${GMP_LIBRARY})
FILE_CMAKE

# 2. include/orthrus/quaternion.hpp - Арифметика кватернионов
cat << 'FILE_QUAT' > include/orthrus/quaternion.hpp
#pragma once

#include <gmpxx.h>

namespace orthrus {

struct Quaternion {
    mpz_class a0{0}; // 1
    mpz_class a1{0}; // i
    mpz_class a2{0}; // j
    mpz_class a3{0}; // k

    [[nodiscard]] Quaternion conjugate() const {
        return {a0, -a1, -a2, -a3};
    }

    [[nodiscard]] mpz_class norm(const mpz_class& p) const {
        return a0 * a0 + a1 * a1 + p * (a2 * a2 + a3 * a3);
    }

    [[nodiscard]] Quaternion multiply(const Quaternion& r, const mpz_class& p) const {
        return {
            a0 * r.a0 - a1 * r.a1 - p * (a2 * r.a2 + a3 * r.a3),
            a0 * r.a1 + a1 * r.a0 + p * (a2 * r.a3 - a3 * r.a2),
            a0 * r.a2 - a1 * r.a3 + a2 * r.a0 + a3 * r.a1,
            a0 * r.a3 + a1 * r.a2 - a2 * r.a1 + a3 * r.a0
        };
    }

    bool operator==(const Quaternion& o) const = default;
};

} // namespace orthrus
FILE_QUAT

# 3. include/orthrus/order.hpp - Максимальный порядок O_0 с базисом Гурвица/Эйхлера
cat << 'FILE_ORDER' > include/orthrus/order.hpp
#pragma once

#include "quaternion.hpp"

namespace orthrus {

// Базис O_0 при p = 3 mod 4:
// w1 = 1, w2 = i, w3 = (1 + j)/2, w4 = (i + k)/2
// q = c1*w1 + c2*w2 + c3*w3 + c4*w4
// 2*q = (2*c1 + c3) + (2*c2 + c4)*i + c3*j + c4*k
struct OrderElement {
    mpz_class c1{0};
    mpz_class c2{0};
    mpz_class c3{0};
    mpz_class c4{0};

    [[nodiscard]] Quaternion to_scaled_quaternion() const {
        return {
            2 * c1 + c3,
            2 * c2 + c4,
            c3,
            c4
        };
    }

    [[nodiscard]] mpz_class norm(const mpz_class& p) const {
        Quaternion q2 = to_scaled_quaternion();
        return q2.norm(p) / 4;
    }

    OrderElement operator+(const OrderElement& o) const {
        return {c1 + o.c1, c2 + o.c2, c3 + o.c3, c4 + o.c4};
    }

    OrderElement operator*(const mpz_class& s) const {
        return {c1 * s, c2 * s, c3 * s, c4 * s};
    }
};

} // namespace orthrus
FILE_ORDER

# 4. include/orthrus/cornacchia.hpp - Диофантов решатель
cat << 'FILE_CORN' > include/orthrus/cornacchia.hpp
#pragma once

#include <gmpxx.h>
#include <optional>
#include <utility>

namespace orthrus {

class Cornacchia {
public:
    static bool mod_sqrt(const mpz_class& a, const mpz_class& p, mpz_class& out) {
        if (a == 0) { out = 0; return true; }
        if (p % 4 == 3) {
            mpz_class exp = (p + 1) / 4;
            mpz_powm(out.get_mpz_t(), a.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
            return true;
        }
        mpz_class q = p - 1;
        unsigned long s = 0;
        while (q % 2 == 0) { q /= 2; s++; }
        mpz_class z = 2;
        while (mpz_legendre(z.get_mpz_t(), p.get_mpz_t()) != -1) z++;
        mpz_class c, r, t, m_exp = (q + 1) / 2;
        mpz_powm(c.get_mpz_t(), z.get_mpz_t(), q.get_mpz_t(), p.get_mpz_t());
        mpz_powm(r.get_mpz_t(), a.get_mpz_t(), m_exp.get_mpz_t(), p.get_mpz_t());
        mpz_powm(t.get_mpz_t(), a.get_mpz_t(), q.get_mpz_t(), p.get_mpz_t());
        unsigned long m_val = s;
        while (t != 1) {
            mpz_class t2 = t;
            unsigned long i = 0;
            for (; i < m_val; ++i) {
                if (t2 == 1) break;
                t2 = (t2 * t2) % p;
            }
            if (i == m_val) return false;
            mpz_class b = c;
            for (unsigned long j = 0; j < m_val - i - 1; ++j) b = (b * b) % p;
            r = (r * b) % p;
            c = (b * b) % p;
            t = (t * c) % p;
            m_val = i;
        }
        out = r;
        return true;
    }

    static std::optional<std::pair<mpz_class, mpz_class>> solve(const mpz_class& d, const mpz_class& m) {
        if (m <= 0 || d <= 0 || d >= m) return std::nullopt;
        mpz_class neg_d = (-d % m + m) % m;
        if (mpz_legendre(neg_d.get_mpz_t(), m.get_mpz_t()) != 1) return std::nullopt;

        mpz_class r0;
        if (!mod_sqrt(neg_d, m, r0)) return std::nullopt;
        if (r0 < m - r0) r0 = m - r0;

        mpz_class r_prev = m, r_curr = r0, limit;
        mpz_sqrt(limit.get_mpz_t(), m.get_mpz_t());
        while (r_curr > limit) {
            mpz_class rem = r_prev % r_curr;
            r_prev = r_curr;
            r_curr = rem;
        }

        mpz_class x = r_curr, rem_sq = m - x * x;
        if (rem_sq < 0 || rem_sq % d != 0) return std::nullopt;
        mpz_class y2 = rem_sq / d;

        if (mpz_perfect_square_p(y2.get_mpz_t()) != 0) {
            mpz_class y;
            mpz_sqrt(y.get_mpz_t(), y2.get_mpz_t());
            return std::make_pair(x, y);
        }
        return std::nullopt;
    }
};

} // namespace orthrus
FILE_CORN

# 5. include/orthrus/ideal.hpp - Базис решётки HNF 4x4
cat << 'FILE_IDEAL' > include/orthrus/ideal.hpp
#pragma once

#include "order.hpp"
#include <array>

namespace orthrus {

// Идеал I максимального порядка O_0 в виде Z-решетки ранга 4
// Задается верхнетреугольной матрицей в канонической форме Эрмита (HNF) 4x4:
// H = [ [ N, 0, h13,   0 ],
//       [ 0, N,   0, h24 ],
//       [ 0, 0,   1,   0 ],
//       [ 0, 0,   0,   1 ] ]
// det(H) = N * N * 1 * 1 = N^2
class IdealHNF {
public:
    mpz_class norm;
    mpz_class h13;
    mpz_class h24;
    mpz_class H13; // 2*h13 + 1
    mpz_class H24; // 2*h24 + 1
    std::array<OrderElement, 4> basis;

    IdealHNF(const mpz_class& N, const mpz_class& X0) : norm(N) {
        // Вычисляем 2^-1 mod N
        mpz_class inv_2;
        mpz_invert(inv_2.get_mpz_t(), mpz_class(2).get_mpz_t(), N.get_mpz_t());

        // h13 = h24 = (X0 - 1) * 2^-1 mod N
        h13 = ((X0 - 1) * inv_2) % N;
        if (h13 < 0) h13 += N;
        h24 = h13;

        H13 = 2 * h13 + 1;
        H24 = 2 * h24 + 1;

        basis[0] = OrderElement{N, 0, 0, 0};
        basis[1] = OrderElement{0, N, 0, 0};
        basis[2] = OrderElement{h13, 0, 1, 0};
        basis[3] = OrderElement{0, h24, 0, 1};
    }

    [[nodiscard]] OrderElement element(const mpz_class& c1, const mpz_class& c2, 
                                       const mpz_class& c3, const mpz_class& c4) const {
        return basis[0] * c1 + basis[1] * c2 + basis[2] * c3 + basis[3] * c4;
    }

    [[nodiscard]] mpz_class determinant() const {
        return basis[0].c1 * basis[1].c2 * basis[2].c3 * basis[3].c4;
    }
};

} // namespace orthrus
FILE_IDEAL

# 6. include/orthrus/klpt.hpp - Цикл декомпозиции Пети (Lattice Sampling & Retry Loop)
cat << 'FILE_KLPT' > include/orthrus/klpt.hpp
#pragma once

#include "ideal.hpp"
#include "cornacchia.hpp"
#include <random>
#include <optional>

namespace orthrus {

class KLPTNavigator {
public:
    // Полноценная декомпозиция Пети: 2D Lattice Sampling + Retry Loop
    static std::optional<OrderElement> find_smooth_generator(
        const IdealHNF& I, 
        unsigned int e, 
        const mpz_class& p, 
        unsigned long max_retries = 500000) 
    {
        mpz_class two_pow_e = mpz_class(1) << e;
        mpz_class target_norm = I.norm * two_pow_e;
        mpz_class two_N = 2 * I.norm;

        gmp_randclass rng(gmp_randinit_default);
        rng.seed(static_cast<unsigned long>(time(nullptr)));

        mpz_class bound;
        mpz_class ratio = (4 * target_norm) / p;
        mpz_sqrt(bound.get_mpz_t(), ratio.get_mpz_t());
        bound /= 2;
        if (bound < 2) bound = 2;

        // --- RETRY LOOP ---
        for (unsigned long attempt = 0; attempt < max_retries; ++attempt) {
            // 1. Фиксируем свободные коэффициенты (c3, c4) - сэмплирование решётки
            mpz_class c3 = rng.get_z_range(bound);
            mpz_class c4 = rng.get_z_range(bound);

            mpz_class p_part = p * (c3 * c3 + c4 * c4);
            if (p_part >= 4 * target_norm) continue;

            mpz_class M_prime = 4 * target_norm - p_part;
            if (M_prime % 4 != 1) continue;
            if (mpz_probab_prime_p(M_prime.get_mpz_t(), 15) == 0) continue;

            // 2. Решаем бинарную квадратичную форму через Корнаккью
            auto xy = Cornacchia::solve(mpz_class(1), M_prime);
            if (!xy.has_value()) continue;

            auto [x, y] = *xy;

            // 3. Подъем в целые числа и сопоставление вычетов по модулю 2N
            std::pair<mpz_class, mpz_class> candidates[8] = {
                { x,  y}, {-x,  y}, { x, -y}, {-x, -y},
                { y,  x}, {-y,  x}, { y, -x}, {-y, -x}
            };

            for (const auto& [cand_X, cand_Y] : candidates) {
                mpz_class diff_X = cand_X - c3 * I.H13;
                mpz_class diff_Y = cand_Y - c4 * I.H24;

                if (diff_X % two_N == 0 && diff_Y % two_N == 0) {
                    mpz_class c1 = diff_X / two_N;
                    mpz_class c2 = diff_Y / two_N;

                    OrderElement gamma = I.element(c1, c2, c3, c4);
                    if (gamma.norm(p) == target_norm) {
                        return gamma; // Точное решение найдено!
                    }
                }
            }
        }

        return std::nullopt;
    }
};

} // namespace orthrus
FILE_KLPT

# 7. src/main.cpp - Точка входа
cat << 'FILE_MAIN' > src/main.cpp
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
FILE_MAIN

# 8. Сборка и тестовый прогон
cmake -B build -S .
cmake --build build
./build/orthrus_test 431 97 32

echo "=================================================================="
echo "[+] ЧАСТЬ 1 успешно внедрена и протестирована!"
echo "=================================================================="
rm -f setup_orthrus_part1.sh
