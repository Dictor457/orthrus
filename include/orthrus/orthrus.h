#ifndef ORTHRUS_H
#define ORTHRUS_H

#ifdef __cplusplus
extern "C" {
#endif

// Коды статуса C-API
typedef enum {
    ORTHRUS_OK = 0,
    ORTHRUS_ERR_INVALID_PARAM = -1,
    ORTHRUS_ERR_NOT_FOUND = -2,
    ORTHRUS_ERR_INTERNAL = -3
} orthrus_status_t;

// Кватернион gamma в представлении: (a0 + a1*i + a2*j + a3*k) / 2
typedef struct {
    char a0[512];
    char a1[512];
    char a2[512];
    char a3[512];
} orthrus_quat_t;

// Основная функция C-API: навигация по идеалу к степени 2^e
orthrus_status_t orthrus_klpt_navigate(
    const char* prime_p_str,
    const char* ideal_norm_N_str,
    unsigned int target_e,
    orthrus_quat_t* out_gamma
);

const char* orthrus_version(void);

#ifdef __cplusplus
}
#endif

#endif // ORTHRUS_H
