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
