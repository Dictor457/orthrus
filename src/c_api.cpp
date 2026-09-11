#include "orthrus/orthrus.h"
#include "orthrus/quaternion.hpp"
#include "orthrus/order.hpp"
#include "orthrus/cornacchia.hpp"
#include "orthrus/ideal.hpp"
#include "orthrus/klpt.hpp"
#include <cstring>

extern "C" {

const char* orthrus_version(void) {
    return "0.3.0-enterprise";
}

orthrus_status_t orthrus_klpt_navigate(
    const char* prime_p_str,
    const char* ideal_norm_N_str,
    unsigned int target_e,
    orthrus_quat_t* out_gamma
) {
    if (!prime_p_str || !ideal_norm_N_str || !out_gamma) {
        return ORTHRUS_ERR_INVALID_PARAM;
    }

    try {
        mpz_class p(prime_p_str);
        mpz_class N(ideal_norm_N_str);

        if (p % 4 != 3 || N <= 1) {
            return ORTHRUS_ERR_INVALID_PARAM;
        }

        mpz_class neg_p = (-p % N + N) % N;
        mpz_class X0;
        if (!orthrus::Cornacchia::mod_sqrt(neg_p, N, X0)) {
            return ORTHRUS_ERR_INVALID_PARAM;
        }

        orthrus::IdealHNF I(N, X0);
        auto sol = orthrus::KLPTNavigator::find_smooth_generator(I, target_e, p);

        if (!sol.has_value()) {
            return ORTHRUS_ERR_NOT_FOUND;
        }

        orthrus::Quaternion scaled = sol->to_scaled_quaternion();

        std::string s_a0 = scaled.a0.get_str();
        std::string s_a1 = scaled.a1.get_str();
        std::string s_a2 = scaled.a2.get_str();
        std::string s_a3 = scaled.a3.get_str();

        std::strncpy(out_gamma->a0, s_a0.c_str(), sizeof(out_gamma->a0) - 1);
        std::strncpy(out_gamma->a1, s_a1.c_str(), sizeof(out_gamma->a1) - 1);
        std::strncpy(out_gamma->a2, s_a2.c_str(), sizeof(out_gamma->a2) - 1);
        std::strncpy(out_gamma->a3, s_a3.c_str(), sizeof(out_gamma->a3) - 1);

        out_gamma->a0[sizeof(out_gamma->a0) - 1] = '\0';
        out_gamma->a1[sizeof(out_gamma->a1) - 1] = '\0';
        out_gamma->a2[sizeof(out_gamma->a2) - 1] = '\0';
        out_gamma->a3[sizeof(out_gamma->a3) - 1] = '\0';

        return ORTHRUS_OK;
    } catch (...) {
        return ORTHRUS_ERR_INTERNAL;
    }
}

} // extern "C"
