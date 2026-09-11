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
