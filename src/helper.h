#pragma once

template<typename T>
struct is_complex : std::false_type {};

template<typename U>
struct is_complex<std::complex<U>> : std::true_type {};

template<typename T>
concept is_complex_v = is_complex<T>::value;

template<typename T>
concept is_scalar_v = is_complex<T>::value || std::is_arithmetic_v<T>;

template<typename T>
constexpr T epsilon() {
    // float, double, long double
    if constexpr (std::is_floating_point_v<T>) {
        return std::numeric_limits<T>::epsilon();
    }
    // bool, short, char, int, long, long long
    else if constexpr (std::is_integral_v<T>) {
        return 1;
    }
    // fallback
    else {
        return static_cast<T>(1e-12);
    }
}

template<is_complex_v T>
constexpr T::value_type epsilon() {
    return epsilon<typename T::value_type>();
}

template<typename T>
bool compare(const T val, const T target) {
    return std::abs(val - target) < epsilon<T>();
}

template<typename T, typename OTHER_T>
bool compare(const T val, const OTHER_T target) {
    return std::abs(val - target) < epsilon<T>();
}
