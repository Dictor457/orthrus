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
