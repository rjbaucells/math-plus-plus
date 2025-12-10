#pragma once
#include "matrix.h"

template<int COLUMNS, int ROWS, is_scalar_v T>
template<int ITER>
Matrix<COLUMNS, ROWS, T>::template LanczosAlgorithm<Matrix<ITER, ITER, T>, Matrix<ITER + 1, COLUMNS, T>> Matrix<COLUMNS, ROWS, T>::lanczosAlgorithm() const requires (isSquare) {
    if (!isHermitian())
        throw std::runtime_error("Cannot do Lanczos algorithm on non hermitian matrix");

    std::array<Vector<COLUMNS, T>, ITER + 1> q;

    Matrix<ITER, ITER, T> t;
    Matrix<ITER + 1, COLUMNS, T> qMatrix;

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
    return vec.dot(this->multiply(vec)) / vec.euclidianNormSquared();
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Vector<COLUMNS, T> Matrix<COLUMNS, ROWS, T>::inverseIteration(const int maxIterations, const T eigenVal, const Vector<COLUMNS, T>& startingVector, const T tolerance) const {
    Vector<COLUMNS, T> b_k = startingVector;
    Matrix<COLUMNS, ROWS, T> thisMinusEigenIdentityInverse = subtract(eigenVal * identity()).inverse();

    for (int k = 0; k < maxIterations; k++) {
        Vector<COLUMNS, T> b_k1 = (thisMinusEigenIdentityInverse * b_k).normalized();

        if (tolerance > 0 && (b_k - b_k1).euclidianNorm() < tolerance)
            return b_k1;

        b_k = b_k1;
    }

    return b_k;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>::template RayleighQuotientIteration<Vector<COLUMNS, T>, T> Matrix<COLUMNS, ROWS, T>::rayleighQuotientIteration(const int maxIterations, const std::optional<T> valueApproximation, const Vector<COLUMNS, T>& vectorApproximation, const T tolerance) const {
    Vector<COLUMNS, T> b_k = vectorApproximation;
    T u_k = valueApproximation.value_or(rayleighQuotient(b_k));

    for (int k = 0; k < maxIterations; k++) {
        Vector<COLUMNS, T> b_k1 = (subtract(u_k * identity()).inverse() * b_k).normalized();
        T u_k1 = rayleighQuotient(b_k1);

        if (tolerance > 0 && (b_k - b_k1).euclidianNorm() < tolerance && std::abs(u_k - u_k1) < tolerance)
            return {b_k1, u_k1};

        b_k = b_k1;
        u_k = u_k1;
    }

    return {b_k, u_k};
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>::template PowerIteration<Vector<COLUMNS, T>, T> Matrix<COLUMNS, ROWS, T>::powerIteration(const int maxIterations, const Vector<COLUMNS, T>& vectorApproximation, const T tolerance) const {
    Vector<COLUMNS, T> b_k = vectorApproximation;
    T u_k = {};

    for (int k = 0; k < maxIterations; k++) {
        Vector<COLUMNS, T> b_k1 = multiply(b_k).normalized();
        T u_k1 = rayleighQuotient(b_k1);

        if (tolerance < 0)
            continue;

        if ((b_k - b_k1).euclidianNorm() < tolerance || std::abs(u_k - u_k1) < tolerance)
            return {b_k1, u_k1};

        b_k = b_k1;
        u_k = u_k1;
    }

    return {b_k, u_k};
}
