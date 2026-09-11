#pragma once

#include <concepts>

namespace orthrus {

template <typename T>
concept BigIntBackend = requires(T a, T b) {
    { a + b } -> std::same_as<T>;
    { a - b } -> std::same_as<T>;
    { a * b } -> std::same_as<T>;
    { a / b } -> std::same_as<T>;
    { a % b } -> std::same_as<T>;
    { a == b } -> std::convertible_to<bool>;
    { a < b }  -> std::convertible_to<bool>;
};

} // namespace orthrus
