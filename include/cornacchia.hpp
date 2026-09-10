#pragma once

#include <gmpxx.h>
#include <optional>
#include <utility>

class CornacchiaSolver {
public:
    // Решает диофантово уравнение: x^2 + d * y^2 = m
    // Для m - простого, 0 < d < m
    static std::optional<std::pair<mpz_class, mpz_class>> solve(const mpz_class& d, const mpz_class& m) {
        if (m <= 0 || d <= 0 || d >= m) {
            return std::nullopt;
        }

        // Проверяем квадратичный вычет: (-d / m) == 1
        mpz_class neg_d = (-d % m + m) % m;
        if (mpz_legendre(neg_d.get_mpz_t(), m.get_mpz_t()) != 1) {
            return std::nullopt;
        }

        // Вычисляем корень r0^2 = -d (mod m)
        mpz_class r0;
        if (!modular_sqrt(neg_d, m, r0)) {
            return std::nullopt;
        }

        if (r0 < m - r0) {
            r0 = m - r0;
        }

        // Евклидово усечение до границы sqrt(m)
        mpz_class r_prev = m;
        mpz_class r_curr = r0;
        mpz_class limit;
        mpz_sqrt(limit.get_mpz_t(), m.get_mpz_t());

        while (r_curr > limit) {
            mpz_class rem = r_prev % r_curr;
            r_prev = r_curr;
            r_curr = rem;
        }

        mpz_class x = r_curr;
        mpz_class rem_sq = m - x * x;

        if (rem_sq < 0 || rem_sq % d != 0) {
            return std::nullopt;
        }

        mpz_class y2 = rem_sq / d;

        // Фикс GMP: mpz_perfect_square_p возвращает не ноль, если y2 - точный квадрат
        if (mpz_perfect_square_p(y2.get_mpz_t()) != 0) {
            mpz_class y;
            mpz_sqrt(y.get_mpz_t(), y2.get_mpz_t());
            return std::make_pair(x, y);
        }

        return std::nullopt;
    }

private:
    // Алгоритм Тоннелли-Шенкса для r^2 = a mod p
    static bool modular_sqrt(const mpz_class& a, const mpz_class& p, mpz_class& out) {
        if (a == 0) {
            out = 0;
            return true;
        }
        if (p == 2) {
            out = a;
            return true;
        }
        if (p % 4 == 3) {
            mpz_class exp = (p + 1) / 4;
            mpz_powm(out.get_mpz_t(), a.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
            return true;
        }

        // p - 1 = Q * 2^S
        mpz_class q = p - 1;
        unsigned long s = 0;
        while (q % 2 == 0) {
            q /= 2;
            s++;
        }

        // Поиск невычета z
        mpz_class z = 2;
        while (mpz_legendre(z.get_mpz_t(), p.get_mpz_t()) != -1) {
            z++;
        }

        mpz_class c, r, t, m_exp;
        mpz_powm(c.get_mpz_t(), z.get_mpz_t(), q.get_mpz_t(), p.get_mpz_t());
        m_exp = (q + 1) / 2;
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
            for (unsigned long j = 0; j < m_val - i - 1; ++j) {
                b = (b * b) % p;
            }

            r = (r * b) % p;
            c = (b * b) % p;
            t = (t * c) % p;
            m_val = i;
        }

        out = r;
        return true;
    }
};
