#pragma once

#include <gmpxx.h>
#include <optional>
#include <utility>

namespace orthrus {

class Cornacchia {
public:
    static bool mod_sqrt(const mpz_class& a, const mpz_class& p, mpz_class& out) {
        if (p <= 1) return false;
        // Защита от зависания: Тоннелли-Шенкс математически определен только над простыми модулями
        if (mpz_probab_prime_p(p.get_mpz_t(), 5) == 0) return false;

        mpz_class a_norm = (a % p + p) % p;
        if (a_norm == 0) { out = 0; return true; }
        if (p == 2) { out = a_norm; return true; }

        if (p % 4 == 3) {
            mpz_class exp = (p + 1) / 4;
            mpz_powm(out.get_mpz_t(), a_norm.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
            return true;
        }

        mpz_class q = p - 1;
        unsigned long s = 0;
        while (q % 2 == 0) { q /= 2; s++; }
        
        // Ограничение поиска невычета (защита от нечетных составных чисел)
        mpz_class z = 2;
        unsigned long search_limit = 128;
        while (z < search_limit && z < p && mpz_legendre(z.get_mpz_t(), p.get_mpz_t()) != -1) {
            z++;
        }
        if (z >= search_limit || z >= p) return false;

        mpz_class c, r, t, m_exp = (q + 1) / 2;
        mpz_powm(c.get_mpz_t(), z.get_mpz_t(), q.get_mpz_t(), p.get_mpz_t());
        mpz_powm(r.get_mpz_t(), a_norm.get_mpz_t(), m_exp.get_mpz_t(), p.get_mpz_t());
        mpz_powm(t.get_mpz_t(), a_norm.get_mpz_t(), q.get_mpz_t(), p.get_mpz_t());
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
        if (m <= 1 || d <= 0 || d >= m) return std::nullopt;

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
