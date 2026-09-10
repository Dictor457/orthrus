#include <iostream>
#include <chrono>
#include <optional>
#include <vector>
#include <string>
#include <iomanip>
#include <gmpxx.h>

struct Quat {
    mpz_class a0{0}, a1{0}, a2{0}, a3{0};

    [[nodiscard]] Quat conj() const { return {a0, -a1, -a2, -a3}; }

    [[nodiscard]] mpz_class norm(const mpz_class& p) const {
        return a0 * a0 + a1 * a1 + p * (a2 * a2 + a3 * a3);
    }

    [[nodiscard]] Quat mul(const Quat& r, const mpz_class& p) const {
        return {
            a0 * r.a0 - a1 * r.a1 - p * (a2 * r.a2 + a3 * r.a3),
            a0 * r.a1 + a1 * r.a0 + p * (a2 * r.a3 - a3 * r.a2),
            a0 * r.a2 - a1 * r.a3 + a2 * r.a0 + a3 * r.a1,
            a0 * r.a3 + a1 * r.a2 - a2 * r.a1 + a3 * r.a0
        };
    }
};

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

class ParameterGenerator {
public:
    // Генерация случайного простого p = 3 mod 4 заданной битности
    static mpz_class generate_prime_3mod4(unsigned int bits, gmp_randclass& rng) {
        mpz_class p;
        while (true) {
            p = rng.get_z_bits(bits);
            p |= (mpz_class(1) << (bits - 1)); // старший бит
            p |= 3; // гарантируем нечетность и младшие биты
            if (p % 4 != 3) p += 2;
            while (p % 4 != 3 || mpz_probab_prime_p(p.get_mpz_t(), 25) == 0) {
                p += 4;
            }
            if (mpz_sizeinbase(p.get_mpz_t(), 2) == bits) return p;
        }
    }
};

class KLPTCore {
public:
    // Полноценный RepresentInteger для любого размера M > p вплоть до 256/512 бит
    static std::optional<Quat> represent_integer(const mpz_class& M, const mpz_class& p, gmp_randclass& rng, unsigned long max_trials = 200000) {
        if (M <= p) return std::nullopt;
        mpz_class ratio = M / p;
        mpz_class bound;
        mpz_sqrt(bound.get_mpz_t(), ratio.get_mpz_t());
        if (bound < 2) bound = 2;

        for (unsigned long iter = 0; iter < max_trials; ++iter) {
            mpz_class z = rng.get_z_range(bound);
            mpz_class w = rng.get_z_range(bound);

            mpz_class p_part = p * (z * z + w * w);
            if (p_part >= M) continue;

            mpz_class M_prime = M - p_part;
            if (M_prime % 4 != 1) continue;
            if (mpz_probab_prime_p(M_prime.get_mpz_t(), 15) == 0) continue;

            auto xy = Cornacchia::solve(mpz_class(1), M_prime);
            if (!xy.has_value()) continue;

            auto [x, y] = *xy;
            Quat q{x, y, z, w};
            if (q.norm(p) == M) return q;
        }
        return std::nullopt;
    }
};

void run_benchmark() {
    std::cout << "================================================================================\n";
    std::cout << "           SQISign KLPT-Engine: Scalability & Benchmark Suite                   \n";
    std::cout << "================================================================================\n";
    std::cout << std::left << std::setw(12) << "Security" 
              << std::setw(16) << "Prime p (bits)" 
              << std::setw(20) << "Target M (bits)" 
              << std::setw(18) << "Time (ms)" 
              << "Status" << "\n";
    std::cout << "--------------------------------------------------------------------------------\n";

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(1337);

    std::vector<unsigned int> bit_sizes = {32, 64, 128, 256};

    for (unsigned int bits : bit_sizes) {
        mpz_class p = ParameterGenerator::generate_prime_3mod4(bits, rng);
        // В SQISign целевая норма M берется кратно больше характеристики p
        mpz_class M = p * (mpz_class(1) << 32) + 1;
        while (M % 4 != 1 || mpz_probab_prime_p(M.get_mpz_t(), 20) == 0) M += 2;

        auto t_start = std::chrono::high_resolution_clock::now();
        auto sol = KLPTCore::represent_integer(M, p, rng);
        auto t_end = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

        std::string sec_level = (bits == 256) ? "NIST-1 (PQC)" : (bits == 128 ? "Mid-Level" : "Toy-Level");

        std::cout << std::left << std::setw(12) << sec_level
                  << std::setw(16) << bits 
                  << std::setw(20) << mpz_sizeinbase(M.get_mpz_t(), 2) 
                  << std::setw(18) << std::fixed << std::setprecision(3) << ms;

        if (sol.has_value() && sol->norm(p) == M) {
            std::cout << "[VERIFIED]" << "\n";
        } else {
            std::cout << "[FAILED]" << "\n";
        }
    }
    std::cout << "================================================================================\n";
}

int main(int argc, char* argv[]) {
    if (argc >= 2 && std::string(argv[1]) == "--benchmark") {
        run_benchmark();
        return 0;
    }

    mpz_class p("431");
    mpz_class N("97");
    unsigned int e = 32;

    if (argc >= 2) p = mpz_class(argv[1]);
    if (argc >= 3) N = mpz_class(argv[2]);
    if (argc >= 4) e = std::stoul(argv[3]);

    std::cout << "==================================================================\n";
    std::cout << "        SQISign KLPT-Engine: C++20 Quaternion Navigation Core     \n";
    std::cout << "==================================================================\n";
    std::cout << "[*] Parameters:\n";
    std::cout << "    Characteristic (p) : " << p << " (p = 3 mod 4)\n";
    std::cout << "    Ideal Norm (N)     : " << N << "\n";
    std::cout << "    Smooth Target (2^e): 2^" << e << "\n";

    mpz_class neg_p = (-p % N + N) % N;
    mpz_class X0;
    if (!Cornacchia::mod_sqrt(neg_p, N, X0)) {
        std::cerr << "[ERROR] -p is not a QR mod N! Select another prime N.\n";
        return 1;
    }

    Quat alpha{X0, 0, 1, 0};
    mpz_class n_alpha = alpha.norm(p);
    std::cout << "\n[+] Ideal Generator: alpha = " << X0 << " + j\n";
    std::cout << "    Norm check: N(alpha) = " << n_alpha << " (" << n_alpha / N << " * N)\n";

    mpz_class two_pow_e = mpz_class(1) << e;
    mpz_class target_norm = N * two_pow_e;
    std::cout << "[*] Target Norm N(gamma) = N * 2^" << e << " = " << target_norm << "\n";

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(static_cast<unsigned long>(time(nullptr)));

    mpz_class max_bound;
    mpz_class ratio = target_norm / p;
    mpz_sqrt(max_bound.get_mpz_t(), ratio.get_mpz_t());
    max_bound /= 2;

    auto t_start = std::chrono::high_resolution_clock::now();
    bool found = false;
    Quat gamma;

    for (unsigned long trial = 0; trial < 1000000; ++trial) {
        mpz_class d0 = rng.get_z_range(max_bound);
        mpz_class d1 = rng.get_z_range(max_bound);

        mpz_class p_part = p * (d0 * d0 + d1 * d1);
        if (p_part >= target_norm) continue;

        mpz_class M_prime = target_norm - p_part;
        if (M_prime % 4 != 1) continue;
        if (mpz_probab_prime_p(M_prime.get_mpz_t(), 15) == 0) continue;

        auto xy = Cornacchia::solve(mpz_class(1), M_prime);
        if (!xy.has_value()) continue;

        auto [x, y] = *xy;

        mpz_class target_c0 = (d0 * X0) % N;
        mpz_class target_c1 = (d1 * X0) % N;

        std::pair<mpz_class, mpz_class> candidates[8] = {
            { (x % N + N) % N,  (y % N + N) % N},
            {(-x % N + N) % N,  (y % N + N) % N},
            { (x % N + N) % N, (-y % N + N) % N},
            {(-x % N + N) % N, (-y % N + N) % N},
            { (y % N + N) % N,  (x % N + N) % N},
            {(-y % N + N) % N,  (x % N + N) % N},
            { (y % N + N) % N, (-x % N + N) % N},
            {(-y % N + N) % N, (-x % N + N) % N}
        };

        mpz_class c0_val[8] = { x, -x,  x, -x,  y, -y,  y, -y };
        mpz_class c1_val[8] = { y,  y, -y, -y,  x,  x, -x, -x };

        for (int k = 0; k < 8; ++k) {
            if (candidates[k].first == target_c0 && candidates[k].second == target_c1) {
                gamma = Quat{c0_val[k], c1_val[k], d0, d1};
                found = true;
                break;
            }
        }

        if (found) break;
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

    if (!found) {
        std::cerr << "[FAIL] Path not found within search bound.\n";
        return 1;
    }

    std::cout << "\n[+] Navigation Solved in " << duration << " ms!\n";
    std::cout << "    gamma = " << gamma.a0 << " + " 
                               << gamma.a1 << "*i + " 
                               << gamma.a2 << "*j + " 
                               << gamma.a3 << "*k\n";
    std::cout << "    Norm(gamma) = " << gamma.norm(p) << "\n";

    Quat prod = gamma.mul(alpha.conj(), p);
    if (prod.a0 % N == 0 && prod.a1 % N == 0 && prod.a2 % N == 0 && prod.a3 % N == 0) {
        std::cout << "[SUCCESS] Ideal membership: gamma * alpha* = 0 (mod " << N << ")\n";
        std::cout << "[SUCCESS] Equivalent smooth ideal J has degree 2^" << e << ".\n";
        std::cout << "          Path length on E_0 is exactly " << e << " 2-isogeny steps.\n";
    } else {
        std::cerr << "[FAIL] Ideal membership violated!\n";
        return 1;
    }

    std::cout << "==================================================================\n";
    return 0;
}
