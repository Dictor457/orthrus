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
