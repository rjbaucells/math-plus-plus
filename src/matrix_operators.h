#pragma once
#include "matrix.h"

// m = m
template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::operator=(const Matrix<COLUMNS, ROWS, T>& other) {
    if (this != &other) {
        memcpy(data, other.data, sizeof(T) * COLUMNS * ROWS);
    }

    return *this;
}

// m == m
template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::equals(const Matrix<COLUMNS, ROWS, T>& other) const {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            if (!compare(data[c][r], other.data[c][r]))
                return false;
        }
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::operator==(const Matrix<COLUMNS, ROWS, T>& other) const {
    return equals(other);
}

// m + m
template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::add(const Matrix<COLUMNS, ROWS, T>& other) const {
    Matrix<COLUMNS, ROWS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = data[c][r] + other.data[c][r];
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::operator+(const Matrix<COLUMNS, ROWS, T>& other) const {
    return add(other);
}

// m - m
template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::subtract(const Matrix<COLUMNS, ROWS, T>& other) const {
    Matrix<COLUMNS, ROWS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = data[c][r] - other.data[c][r];
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::operator-(const Matrix<COLUMNS, ROWS, T>& other) const {
    return subtract(other);
}

// m * m
template<int COLUMNS, int ROWS, is_scalar_v T>
template<int OTHER_COLUMNS>
Matrix<OTHER_COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::multiply(const Matrix<OTHER_COLUMNS, COLUMNS, T>& other) const {
    Matrix<OTHER_COLUMNS, ROWS, T> result;

    for (int c = 0; c < OTHER_COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            for (int x = 0; x < COLUMNS; x++) {
                result[c][r] += data[x][r] * other.data[c][x];
            }
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<int OTHER_COLUMNS>
Matrix<OTHER_COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::operator*(const Matrix<OTHER_COLUMNS, COLUMNS, T>& other) const {
    return multiply(other);
}

// m * v
template<int COLUMNS, int ROWS, is_scalar_v T>
Vector<COLUMNS, T> Matrix<COLUMNS, ROWS, T>::multiply(const Vector<COLUMNS, T>& other) const {
    Vector<COLUMNS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[r] += data[c][r] * other[c];
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Vector<COLUMNS, T> Matrix<COLUMNS, ROWS, T>::operator*(const Vector<COLUMNS, T>& other) const {
    return multiply(other);
}

// m * #
template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::multiply(const T val) const {
    Matrix<COLUMNS, ROWS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = data[c][r] * val;
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::operator*(const T val) const {
    return multiply(val);
}

// m / #
template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::divide(const T scalar) const {
    Matrix<COLUMNS, ROWS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = data[c][r] / scalar;
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::operator/(const T scalar) const {
    return divide(scalar);
}

// m += m
template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::addEquals(const Matrix<COLUMNS, ROWS, T>& other) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            data[c][r] += other.data[c][r];
        }
    }

    return *this;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::operator+=(const Matrix<COLUMNS, ROWS, T>& other) {
    return addEquals(other);
}

// m -= m
template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::subtractEquals(const Matrix<COLUMNS, ROWS, T>& other) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            data[c][r] -= other.data[c][r];
        }
    }

    return *this;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::operator-=(const Matrix<COLUMNS, ROWS, T>& other) {
    return subtractEquals(other);
}

// m *= #
template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::multiplyEquals(const T val) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            data[c][r] *= val;
        }
    }

    return *this;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::operator*=(const T val) {
    return multiplyEquals(val);
}

// m /= #
template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::divideEquals(const T scalar) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            data[c][r] /= scalar;
        }
    }

    return *this;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::operator/=(const T scalar) {
    return divideEquals(scalar);
}

// Operators for different types

// m = m
template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires std::convertible_to<OTHER_T, T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::operator=(const Matrix<COLUMNS, ROWS, OTHER_T>& other) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            data[c][r] = other.data[c][r];
        }
    }

    return *this;
}

// m == m
template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires std::equality_comparable_with<OTHER_T, T>
bool Matrix<COLUMNS, ROWS, T>::equals(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            if (!compare(data[c][r], other.data[c][r]))
                return false;
        }
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires std::equality_comparable_with<OTHER_T, T>
bool Matrix<COLUMNS, ROWS, T>::operator==(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const {
    return equals(other);
}

// m + m
template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires has_common_type<OTHER_T, T>
Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::add(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const {
    Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = data[c][r] + other.data[c][r];
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires has_common_type<OTHER_T, T>
Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::operator+(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const {
    return add(other);
}

// m - m
template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires has_common_type<OTHER_T, T>
Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::subtract(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const {
    Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = data[c][r] - other.data[c][r];
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires has_common_type<OTHER_T, T>
Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::operator-(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const {
    return subtract(other);
}

// m * m
template<int COLUMNS, int ROWS, is_scalar_v T>
template<int OTHER_COLUMNS, typename OTHER_T> requires has_common_type<OTHER_T, T>
Matrix<OTHER_COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::multiply(const Matrix<OTHER_COLUMNS, COLUMNS, OTHER_T>& other) const {
    Matrix<OTHER_COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> result;

    for (int c = 0; c < OTHER_COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            for (int x = 0; x < COLUMNS; x++) {
                result[c][r] += data[x][r] * other.data[c][x];
            }
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<int OTHER_COLUMNS, typename OTHER_T> requires has_common_type<OTHER_T, T>
Matrix<OTHER_COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::operator*(const Matrix<OTHER_COLUMNS, COLUMNS, OTHER_T>& other) const {
    return multiply(other);
}

// m * v
template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires has_common_type<OTHER_T, T>
Vector<COLUMNS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::multiply(const Vector<COLUMNS, OTHER_T>& other) const {
    Vector<COLUMNS, std::common_type_t<T, OTHER_T>> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[r] += data[c][r] * other[c];
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires has_common_type<OTHER_T, T>
Vector<COLUMNS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::operator*(const Vector<COLUMNS, OTHER_T>& other) const {
    return multiply(other);
}

// m * #
template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires has_common_type<OTHER_T, T>
Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::multiply(const OTHER_T val) const {
    Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = data[c][r] * val;
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires has_common_type<OTHER_T, T>
Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::operator*(const OTHER_T val) const {
    return multiply(val);
}

// m / #
template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires has_common_type<OTHER_T, T>
Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::divide(const OTHER_T scalar) const {
    Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = data[c][r] / scalar;
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires has_common_type<OTHER_T, T>
Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::operator/(const OTHER_T scalar) const {
    return divide(scalar);
}

// m += m
template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires std::convertible_to<OTHER_T, T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::addEquals(const Matrix<COLUMNS, ROWS, OTHER_T>& other) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            data[c][r] += other.data[c][r];
        }
    }

    return *this;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires std::convertible_to<OTHER_T, T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::operator+=(const Matrix<COLUMNS, ROWS, OTHER_T>& other) {
    return addEquals(other);
}

// m -= m
template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires std::convertible_to<OTHER_T, T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::subtractEquals(const Matrix<COLUMNS, ROWS, OTHER_T>& other) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            data[c][r] -= other.data[c][r];
        }
    }

    return *this;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires std::convertible_to<OTHER_T, T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::operator-=(const Matrix<COLUMNS, ROWS, OTHER_T>& other) {
    return subtractEquals(other);
}

// m *= #
template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires std::convertible_to<OTHER_T, T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::multiplyEquals(const OTHER_T val) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            data[c][r] *= val;
        }
    }

    return *this;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires std::convertible_to<OTHER_T, T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::operator*=(const OTHER_T val) {
    return multiplyEquals(val);
}

// m /= #
template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires std::convertible_to<OTHER_T, T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::divideEquals(const OTHER_T scalar) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            data[c][r] /= scalar;
        }
    }

    return *this;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T>
requires std::convertible_to<OTHER_T, T>
Matrix<COLUMNS, ROWS, T>& Matrix<COLUMNS, ROWS, T>::operator/=(const OTHER_T scalar) {
    return divideEquals(scalar);
}

template<int COLUMNS, int ROWS, is_scalar_v T>
T* Matrix<COLUMNS, ROWS, T>::operator[](const int index) {
    return &data[index][0];
}

template<int COLUMNS, int ROWS, is_scalar_v T>
const T* Matrix<COLUMNS, ROWS, T>::operator[](const int index) const {
    return &data[index][0];
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::operator-() const {
    Matrix<COLUMNS, ROWS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = -data[c][r];
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>::operator T*() {
    return &data[0][0];
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>::operator const T*() const {
    return &data[0][0];
}
