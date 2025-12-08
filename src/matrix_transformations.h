#pragma once
#include "matrix.h"

template<int COLUMNS, int ROWS, is_scalar_v T>
template<int N>
Vector<N, T> Matrix<COLUMNS, ROWS, T>::applyHomogeneousTransformation(const Vector<N, T>& point) const requires (isSquare) {
    Vector<COLUMNS, T> resizedPoint;

    for (int i = 0; i < N; i++) {
        resizedPoint[i] = point[i];
    }

    Vector<COLUMNS, T> transformedPoint = multiply(resizedPoint);

    Vector<N, T> result;

    for (int i = 0; i < N; i++) {
        result[i] = transformedPoint[i];
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<int N, typename OTHER_T> requires has_common_type<OTHER_T, T>
Vector<N, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::applyHomogeneousTransformation(const Vector<N, OTHER_T>& point) const requires (isSquare) {
    Vector<COLUMNS, std::common_type_t<T, OTHER_T>> resizedPoint;

    for (int i = 0; i < N; i++) {
        resizedPoint[i] = point[i];
    }

    Vector<COLUMNS, std::common_type_t<T, OTHER_T>> transformedPoint = multiply(resizedPoint);

    Vector<N, std::common_type_t<T, OTHER_T>> result;

    for (int i = 0; i < N; i++) {
        result[i] = transformedPoint[i];
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::scalingMatrix(const Vector<COLUMNS, T>& factors) requires (isSquare) {
    Matrix<COLUMNS, ROWS, T> matrix;

    for (int c = 0; c < COLUMNS; c++) {
        matrix[c][c] = factors[c];
    }

    return matrix;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::shearMatrix(const int i, const int j, const T k) requires (isSquare) {
    Matrix<COLUMNS, ROWS, T> matrix = identity();

    matrix[j][i] = k;

    return matrix;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::squeezeMatrix(const int i, const int j, const T k) requires (isSquare) {
    Matrix<COLUMNS, ROWS, T> matrix = identity();

    matrix[i][i] = k;
    matrix[j][j] = 1 / k;

    return matrix;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::rotationMatrixAboutOrigin(const T rot, const RotationType rotationType) requires (isSquare && COLUMNS == 2) {
    T asRadians = convert(rotationType, RotationType::radians, rot);

    T sin = std::sin(asRadians);
    T cos = std::cos(asRadians);

    Matrix<COLUMNS, ROWS, T> r;

    r[0][0] = cos;
    r[1][0] = -sin;
    r[0][1] = sin;
    r[1][1] = cos;

    return r;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS + 1, ROWS + 1, T> Matrix<COLUMNS, ROWS, T>::rotationMatrixAboutPoint(const Vector<COLUMNS, T>& p, const T rot, const RotationType rotationType) requires (isSquare && COLUMNS == 2) {
    Matrix<COLUMNS, ROWS, T> rotationMatrix = rotationMatrixAboutOrigin(rot, rotationType);
    Vector<COLUMNS, T> translationVector = (identity() - rotationMatrix) * p;

    Matrix<COLUMNS + 1, ROWS + 1, T> result = Matrix<COLUMNS + 1, ROWS + 1, T>::identity();

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = rotationMatrix[c][r];
        }
    }

    for (int r = 0; r < ROWS; r++) {
        result[COLUMNS][r] = translationVector[r];
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::rotationMatrixAroundAxisThroughOrigin(const Vector<COLUMNS, T>& axis, const T rot, const RotationType rotationType) requires (isSquare && COLUMNS == 3) {
    T asRadians = convert(rotationType, RotationType::radians, rot);

    T sin = std::sin(asRadians);
    T cos = std::cos(asRadians);

    Matrix<COLUMNS, ROWS, T> r = identity() * cos + (1 - cos) * axis.outerProduct(axis) + axis.crossProductMatrix() * sin;

    return r;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS + 1, ROWS + 1, T> Matrix<COLUMNS, ROWS, T>::rotationMatrixAroundAxisNotThroughOrigin(const Vector<COLUMNS, T>& axis, const Vector<COLUMNS, T>& point, const T rot, const RotationType rotationType) requires (isSquare && COLUMNS == 3) {
    Vector<COLUMNS, T> u = axis.normalize();

    Matrix<COLUMNS, ROWS, T> rotationMatrix = rotationMatrixAroundAxisThroughOrigin(u, rot, rotationType);
    Vector<COLUMNS, T> translationVector = (identity() - rotationMatrix) * point;

    Matrix<COLUMNS + 1, ROWS + 1, T> result = Matrix<COLUMNS + 1, ROWS + 1, T>::identity();

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = rotationMatrix[c][r];
        }
    }

    for (int r = 0; r < ROWS; r++) {
        result[COLUMNS][r] = translationVector[r];
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::rotationMatrixInPlaneThroughOrigin(const Vector<COLUMNS, T>& v1, const Vector<COLUMNS, T>& v2, const T rot, const RotationType rotationType) requires (isSquare && COLUMNS >= 3) {
    T asRadians = convert(rotationType, RotationType::radians, rot);

    T sin = std::sin(asRadians);
    T cos = std::cos(asRadians);

    auto orthonormalized = Vector<COLUMNS, T>::orthonormalize({v1, v2});

    Vector<COLUMNS, T> u = orthonormalized[0];
    Vector<COLUMNS, T> v = orthonormalized[1];

    Matrix<COLUMNS, ROWS, T> r = identity() + (cos - 1) * (u.outerProduct(u) + v.outerProduct(v)) + sin * (v.outerProduct(u) - u.outerProduct(v));

    return r;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::rotationMatrixInPLaneNotThroughOrigin(const Vector<COLUMNS, T>& v1, const Vector<COLUMNS, T>& v2, const Vector<COLUMNS, T>& point, const T rot, const RotationType rotationType) requires (isSquare && COLUMNS >= 3) {
    auto orthonormalized = Vector<COLUMNS, T>::orthonormalize({v1, v2});

    Vector<COLUMNS, T> u = orthonormalized[0];
    Vector<COLUMNS, T> v = orthonormalized[1];

    Matrix<COLUMNS, ROWS, T> rotationMatrix = rotationMatrixInPlaneThroughOrigin(u, v, rot, rotationType);

    Vector<COLUMNS, T> translationVector = (identity() - rotationMatrix) * point;

    Matrix<COLUMNS + 1, ROWS + 1, T> result = Matrix<COLUMNS + 1, ROWS + 1, T>::identity();

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = rotationMatrix[c][r];
        }
    }

    for (int r = 0; r < ROWS; r++) {
        result[COLUMNS][r] = translationVector[r];
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::reflectionMatrixAlongAxisThroughOrigin(const Vector<COLUMNS, T>& axis) requires (isSquare) {
    return 2 * axis.outerProduct(axis) - identity();
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS + 1, ROWS + 1, T> Matrix<COLUMNS, ROWS, T>::reflectionMatrixAlongAxisNotThroughOrigin(const Vector<COLUMNS, T>& axis, const Vector<COLUMNS, T>& point) requires (isSquare) {
    Matrix<COLUMNS, ROWS, T> reflectionMatrix = reflectionMatrixAlongAxisThroughOrigin(axis);
    Vector<COLUMNS, T> translationVector = (identity() - 2 * axis.outerProduct(axis)) * point;

    Matrix<COLUMNS + 1, ROWS + 1, T> result = Matrix<COLUMNS + 1, ROWS + 1, T>::identity();

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = reflectionMatrix[c][r];
        }
    }

    for (int r = 0; r < ROWS; r++) {
        result[COLUMNS][r] = translationVector[r];
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS + 1, ROWS + 1, T> Matrix<COLUMNS, ROWS, T>::translationMatrix(const Vector<COLUMNS, T>& translation) requires (isSquare) {
    Matrix<COLUMNS + 1, ROWS + 1, T> matrix = Matrix<COLUMNS + 1, ROWS + 1, T>::identity();

    for (int r = 0; r < ROWS; r++) {
        matrix[COLUMNS][r] = translation[r];
    }

    return matrix;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::orthoMatrix(const T left, const T right, const T bottom, const T top, const T near, const T far) requires (isSquare && COLUMNS == 4) {
    // identity
    Matrix<COLUMNS, ROWS, T> transformation = Matrix<COLUMNS, ROWS, T>::identity();
    // transformation
    transformation.data[0][0] = 2 / (right - left);
    transformation.data[1][1] = 2 / (top - bottom);
    transformation.data[2][2] = -2 / (far - near);
    transformation.data[3][0] = -(right + left) / (right - left);
    transformation.data[3][1] = -(top + bottom) / (top - bottom);
    transformation.data[3][2] = -(far + near) / (far - near);
    // return
    return transformation;
}
