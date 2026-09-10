#pragma once

#include "quaternion.hpp"

struct OrderElement {
    mpz_class c0{0};
    mpz_class c1{0};
    mpz_class c2{0};
    mpz_class c3{0};

    [[nodiscard]] Quaternion to_scaled_quaternion() const {
        return {
            2 * c0 + c3, // X
            2 * c1 + c2, // Y
            c2,          // Z
            c3           // W
        };
    }

    [[nodiscard]] mpz_class norm(const mpz_class& p) const {
        Quaternion q2 = to_scaled_quaternion();
        return q2.norm(p) / 4;
    }

    static bool is_in_order(const mpz_class& X, const mpz_class& Y, const mpz_class& Z, const mpz_class& W) {
        return ((X - W) % 2 == 0) && ((Y - Z) % 2 == 0);
    }

    static OrderElement from_scaled(const mpz_class& X, const mpz_class& Y, const mpz_class& Z, const mpz_class& W) {
        return {
            (X - W) / 2,
            (Y - Z) / 2,
            Z,
            W
        };
    }

    // Умножение двух элементов порядка O_0: (q1 * q2) in O_0
    [[nodiscard]] OrderElement multiply(const OrderElement& rhs, const mpz_class& p) const {
        Quaternion q1_scaled = to_scaled_quaternion();
        Quaternion q2_scaled = rhs.to_scaled_quaternion();
        Quaternion prod_scaled = q1_scaled.multiply(q2_scaled, p); // 4 * (q1 * q2)
        // 2 * (q1 * q2) = prod_scaled / 2
        return OrderElement::from_scaled(prod_scaled.a0 / 2, prod_scaled.a1 / 2, prod_scaled.a2 / 2, prod_scaled.a3 / 2);
    }

    OrderElement operator+(const OrderElement& o) const {
        return {c0 + o.c0, c1 + o.c1, c2 + o.c2, c3 + o.c3};
    }

    OrderElement operator*(const mpz_class& scalar) const {
        return {c0 * scalar, c1 * scalar, c2 * scalar, c3 * scalar};
    }
};
