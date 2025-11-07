#pragma once
#include <sstream>
#include <string>

template<typename T>
struct Complex {
    T real;
    T imaginary;

    Complex<T> complexConjugate() const {
        Complex<T> cc = *this;
        cc.imaginary *= -1;
        return cc;
    }

    [[nodiscard]] std::string toString() const {
        std::stringstream ss;
        ss.precision(2);

        ss << std::to_string(real) << " " << std::to_string(imaginary) << "i";

        return ss.str();
    }
};
