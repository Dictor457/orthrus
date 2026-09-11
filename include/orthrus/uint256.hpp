#pragma once

#include <cstdint>
#include <array>
#include <string>
#include <iostream>
#include <compare>
#include <iomanip>
#include <sstream>
#include <x86intrin.h>

namespace orthrus {

// 256-битное беззнаковое число с фиксированным размером 32 байта
// Полностью стековое размещение: 0 аллокаций памяти, 0 malloc/free
struct alignas(32) uint256 {
    uint64_t limbs[4]{0, 0, 0, 0}; // Little-endian: limbs[0] - младшие 64 бита

    constexpr uint256() = default;
    constexpr uint256(uint64_t low) : limbs{low, 0, 0, 0} {}
    constexpr uint256(uint64_t l0, uint64_t l1, uint64_t l2, uint64_t l3) : limbs{l0, l1, l2, l3} {}

    // Аппаратное сложение с флагом переноса (Zen 3 ADCX/ADOX)
    [[nodiscard]] uint256 operator+(const uint256& b) const {
        uint256 r;
        unsigned char carry = 0;
        carry = _addcarry_u64(carry, limbs[0], b.limbs[0], (unsigned long long*)&r.limbs[0]);
        carry = _addcarry_u64(carry, limbs[1], b.limbs[1], (unsigned long long*)&r.limbs[1]);
        carry = _addcarry_u64(carry, limbs[2], b.limbs[2], (unsigned long long*)&r.limbs[2]);
        _addcarry_u64(carry, limbs[3], b.limbs[3], (unsigned long long*)&r.limbs[3]);
        return r;
    }

    // Аппаратное вычитание с заёмом
    [[nodiscard]] uint256 operator-(const uint256& b) const {
        uint256 r;
        unsigned char borrow = 0;
        borrow = _subborrow_u64(borrow, limbs[0], b.limbs[0], (unsigned long long*)&r.limbs[0]);
        borrow = _subborrow_u64(borrow, limbs[1], b.limbs[1], (unsigned long long*)&r.limbs[1]);
        borrow = _subborrow_u64(borrow, limbs[2], b.limbs[2], (unsigned long long*)&r.limbs[2]);
        _subborrow_u64(borrow, limbs[3], b.limbs[3], (unsigned long long*)&r.limbs[3]);
        return r;
    }

    // Сравнение в C++20 через spaceship-оператор
    bool operator==(const uint256& b) const = default;

    std::strong_ordering operator<=>(const uint256& b) const {
        for (int i = 3; i >= 0; --i) {
            if (limbs[i] < b.limbs[i]) return std::strong_ordering::less;
            if (limbs[i] > b.limbs[i]) return std::strong_ordering::greater;
        }
        return std::strong_ordering::equal;
    }

    // Побитовые сдвиги
    [[nodiscard]] uint256 operator<<(unsigned int shift) const {
        if (shift == 0) return *this;
        if (shift >= 256) return uint256(0);
        uint256 r;
        unsigned int word_shift = shift / 64;
        unsigned int bit_shift = shift % 64;
        for (unsigned int i = 0; i < 4; ++i) {
            if (i + word_shift < 4) {
                r.limbs[i + word_shift] |= (limbs[i] << bit_shift);
                if (bit_shift > 0 && i + word_shift + 1 < 4) {
                    r.limbs[i + word_shift + 1] |= (limbs[i] >> (64 - bit_shift));
                }
            }
        }
        return r;
    }

    [[nodiscard]] uint256 operator>>(unsigned int shift) const {
        if (shift == 0) return *this;
        if (shift >= 256) return uint256(0);
        uint256 r;
        unsigned int word_shift = shift / 64;
        unsigned int bit_shift = shift % 64;
        for (int i = 3; i >= 0; --i) {
            if (i - (int)word_shift >= 0) {
                r.limbs[i - word_shift] |= (limbs[i] >> bit_shift);
                if (bit_shift > 0 && i - (int)word_shift - 1 >= 0) {
                    r.limbs[i - word_shift - 1] |= (limbs[i] << (64 - bit_shift));
                }
            }
        }
        return r;
    }

    // Умножение 256x256 -> 256 на базе 128-битных перекрестных произведений (MULX)
    [[nodiscard]] uint256 operator*(const uint256& b) const {
        uint256 r;
        unsigned __int128 carry = 0;
        for (size_t i = 0; i < 4; ++i) {
            carry = 0;
            for (size_t j = 0; i + j < 4; ++j) {
                unsigned __int128 cur = (unsigned __int128)r.limbs[i + j] +
                                       (unsigned __int128)limbs[i] * (unsigned __int128)b.limbs[j] + carry;
                r.limbs[i + j] = (uint64_t)cur;
                carry = cur >> 64;
            }
        }
        return r;
    }

    // Деление и остаток
    static std::pair<uint256, uint256> divmod(const uint256& a, const uint256& b) {
        if (b == uint256(0)) throw std::runtime_error("uint256 division by zero");
        if (a < b) return {uint256(0), a};
        if (a == b) return {uint256(1), uint256(0)};

        uint256 quotient = 0;
        uint256 remainder = 0;

        for (int i = 255; i >= 0; --i) {
            remainder = remainder << 1;
            unsigned int word_idx = static_cast<unsigned int>(i / 64);
            unsigned int bit_idx = static_cast<unsigned int>(i % 64);
            if ((a.limbs[word_idx] >> bit_idx) & 1ULL) {
                remainder.limbs[0] |= 1ULL;
            }
            if (remainder >= b) {
                remainder = remainder - b;
                quotient.limbs[word_idx] |= (1ULL << bit_idx);
            }
        }
        return {quotient, remainder};
    }

    [[nodiscard]] uint256 operator/(const uint256& b) const { return divmod(*this, b).first; }
    [[nodiscard]] uint256 operator%(const uint256& b) const { return divmod(*this, b).second; }

    [[nodiscard]] std::string to_hex() const {
        std::stringstream ss;
        ss << "0x";
        bool leading = true;
        for (int i = 3; i >= 0; --i) {
            if (leading && limbs[i] == 0 && i != 0) continue;
            if (!leading) {
                ss << std::hex << std::setw(16) << std::setfill('0') << limbs[i];
            } else {
                ss << std::hex << limbs[i];
                leading = false;
            }
        }
        return ss.str();
    }
};

} // namespace orthrus
