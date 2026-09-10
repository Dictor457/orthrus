#pragma once

#include <vector>
#include <array>
#include "order.hpp"

class LeftIdeal {
public:
    OrderElement alpha;
    mpz_class norm;
    // 4 базисных элемента O_0, образующих Z-решетку идеала I
    std::array<OrderElement, 4> basis;

    LeftIdeal(OrderElement gen, mpz_class n, const mpz_class& p) 
        : alpha(std::move(gen)), norm(std::move(n)) 
    {
        compute_hnf_basis(p);
    }

    // Получение элемента идеала по координатам решетки: gamma = sum(x_i * b_i) in I
    [[nodiscard]] OrderElement element_from_lattice(const mpz_class& x0, const mpz_class& x1, 
                                                    const mpz_class& x2, const mpz_class& x3) const {
        return basis[0] * x0 + basis[1] * x1 + basis[2] * x2 + basis[3] * x3;
    }

    [[nodiscard]] bool contains(const OrderElement& gamma, const mpz_class& p) const {
        Quaternion q_gamma = gamma.to_scaled_quaternion();
        Quaternion q_alpha = alpha.to_scaled_quaternion();
        Quaternion q_alpha_conj = q_alpha.conjugate();

        Quaternion prod = q_gamma.multiply(q_alpha_conj, p);
        mpz_class mod = 4 * norm;

        return (prod.a0 % mod == 0) && 
               (prod.a1 % mod == 0) && 
               (prod.a2 % mod == 0) && 
               (prod.a3 % mod == 0);
    }

private:
    // Построение канонического 4x4 HNF базиса из 8 образующих
    void compute_hnf_basis(const mpz_class& p) {
        std::array<OrderElement, 4> w = {
            OrderElement{1, 0, 0, 0}, // w0 = 1
            OrderElement{0, 1, 0, 0}, // w1 = i
            OrderElement{0, 0, 1, 0}, // w2 = (i+j)/2
            OrderElement{0, 0, 0, 1}  // w3 = (1+k)/2
        };

        // 8 образующих идеала: w_i * alpha и N * w_i
        std::vector<std::array<mpz_class, 4>> M(8);
        for (size_t i = 0; i < 4; ++i) {
            OrderElement prod = w[i].multiply(alpha, p);
            M[i] = {prod.c0, prod.c1, prod.c2, prod.c3};
            M[i + 4] = {0, 0, 0, 0};
            M[i + 4][i] = norm;
        }

        // Целочисленное приведение матрицы 8x4 к верхнетреугольной HNF
        size_t lead = 0;
        for (size_t col = 0; col < 4 && lead < 8; ++col) {
            while (true) {
                size_t pivot = 8;
                for (size_t r = lead; r < 8; ++r) {
                    if (M[r][col] != 0) {
                        if (pivot == 8 || abs(M[r][col]) < abs(M[pivot][col])) {
                            pivot = r;
                        }
                    }
                }
                if (pivot == 8) break; // Все нули в столбце

                std::swap(M[lead], M[pivot]);
                if (M[lead][col] < 0) {
                    for (size_t c = 0; c < 4; ++c) M[lead][c] = -M[lead][c];
                }

                bool all_zero = true;
                for (size_t r = lead + 1; r < 8; ++r) {
                    if (M[r][col] != 0) {
                        mpz_class q = M[r][col] / M[lead][col];
                        for (size_t c = 0; c < 4; ++c) {
                            M[r][c] -= q * M[lead][c];
                        }
                        if (M[r][col] != 0) all_zero = false;
                    }
                }
                if (all_zero) break;
            }

            // Редукция верхних строк
            for (size_t r = 0; r < lead; ++r) {
                mpz_class q = M[r][col] / M[lead][col];
                for (size_t c = 0; c < 4; ++c) {
                    M[r][c] -= q * M[lead][c];
                }
                if (M[r][col] < 0) {
                    for (size_t c = 0; c < 4; ++c) M[r][c] += M[lead][c];
                }
            }
            lead++;
        }

        for (size_t i = 0; i < 4; ++i) {
            basis[i] = OrderElement{M[i][0], M[i][1], M[i][2], M[i][3]};
        }
    }
};
