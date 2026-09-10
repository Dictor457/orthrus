#pragma once

#include "order.hpp"
#include "ideal.hpp"
#include "cornacchia.hpp"

class KLPTNavigator {
public:
    // Поиск элемента gamma in I такого, что N(gamma) = N(I) * 2^e
    static std::optional<OrderElement> find_smooth_generator(const LeftIdeal& I, 
                                                             unsigned int e, 
                                                             const mpz_class& p, 
                                                             unsigned long max_trials = 200000) {
        mpz_class two_pow_e;
        mpz_ui_pow_ui(two_pow_e.get_mpz_t(), 2, e);
        mpz_class target_norm = I.norm * two_pow_e;

        gmp_randclass rng(gmp_randinit_default);
        rng.seed(static_cast<unsigned long>(time(nullptr)));

        // Ограничение диапазона координат решетки
        mpz_class search_bound;
        mpz_class ratio = target_norm / I.norm;
        mpz_sqrt(search_bound.get_mpz_t(), ratio.get_mpz_t());
        if (search_bound < 2) search_bound = 2;

        for (unsigned long trial = 0; trial < max_trials; ++trial) {
            // Сэмплируем вектор решетки (x2, x3)
            mpz_class x2 = rng.get_z_range(search_bound) - (search_bound / 2);
            mpz_class x3 = rng.get_z_range(search_bound) - (search_bound / 2);

            OrderElement partial = I.basis[2] * x2 + I.basis[3] * x3;
            mpz_class partial_norm = partial.norm(p);
            if (partial_norm >= target_norm) continue;

            mpz_class rem_norm = target_norm - partial_norm;
            if (rem_norm % 4 != 1) continue;
            if (mpz_probab_prime_p(rem_norm.get_mpz_t(), 15) == 0) continue;

            // Разрешаем остаточную 2D-форму через Корнаккью
            auto xy = CornacchiaSolver::solve(mpz_class(1), rem_norm);
            if (!xy.has_value()) continue;

            auto [x0, x1] = *xy;
            OrderElement candidate = I.element_from_lattice(x0, x1, x2, x3);

            if (candidate.norm(p) == target_norm && I.contains(candidate, p)) {
                return candidate;
            }
        }

        return std::nullopt;
    }
};
