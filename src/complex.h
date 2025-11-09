#pragma once

template<typename T>
struct IsComplex : std::false_type {};

template<typename U>
struct IsComplex<std::complex<U>> : std::true_type {};