#pragma once
#include "matrix.h"

template<int COLUMNS, int ROWS, is_scalar_v T>
template<int ITER>
Matrix<COLUMNS, ROWS, T>::template LanczosAlgorithm<Matrix<ITER, ITER, T>, Matrix<ITER + 1, COLUMNS, T>> Matrix<COLUMNS, ROWS, T>::lanczosAlgorithm() const requires (isSquare) {
    if (!isHermitian())
        throw std::runtime_error("Cannot do Lanczos algorithm on non hermitian matrix");

    std::array<Vector<COLUMNS, T>, ITER + 1> q;

    Matrix < ITER, ITER, T > t;
    Matrix < ITER + 1, COLUMNS, T > qMatrix;

    qMatrix.setColumnVector(0, Vector<COLUMNS, T>::random().normalize());

    for (int m = 0; m < ITER; m++) {
        Vector<COLUMNS, T> v = multiply(q[m]);
        t[m][m] = q[m] * v;

        if (m == 0) {
            v -= t[m][m] * q[m];
        }
        else {
            v -= t[m][m - 1] * q[m - 1] - t[m][m] * q[m];
        }

        T vNorm = v.euclidianNorm();

        if (m != ITER - 1) {
            t[m][m + 1] = t[m + 1][m] = vNorm;
        }

        q[m + 1] = v / vNorm;

        qMatrix.setColumnVector(m + 1, q[m + 1]);
    }

    return {t, qMatrix};
}

template<int COLUMNS, int ROWS, is_scalar_v T>
T Matrix<COLUMNS, ROWS, T>::rayleighQuotient(const Vector<COLUMNS, T>& vec) const {
    return (vec * *this * vec) / (vec * vec);
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>::template PowerIteration<Vector<COLUMNS, T>, T> Matrix<COLUMNS, ROWS, T>::powerIteration(const int maxIterations, const T tolerance) const {
    Vector<COLUMNS, T> vec = Vector<COLUMNS, T>::random();
    T val = {};

    for (int i = 0; i < maxIterations; i++) {
        vec = multiply(vec).normalize();
        T nextVal = rayleighQuotient(vec);

        if (std::abs(val - nextVal) < tolerance) {
            return {vec, nextVal};
        }

        val = nextVal;
    }

    return {vec, val};
}
