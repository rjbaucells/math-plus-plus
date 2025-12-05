#pragma once

#include <cassert>
#include <complex>
#include <regex>
#include <stdexcept>
#include <string>
#include "helper.h"
#include "vector.h"

template<int COLUMNS, int ROWS, is_scalar_v T>
struct Matrix;

template<typename T>
struct is_matrix : std::false_type {};

template<int COLUMNS, int ROWS, is_scalar_v T>
struct is_matrix<Matrix<COLUMNS, ROWS, T>> : std::true_type {};

template<typename T>
concept is_matrix_v = is_matrix<T>::value;

template<int COLUMNS, int ROWS, is_scalar_v T = float>
struct Matrix {
    static constexpr int columns = COLUMNS;
    static constexpr int rows = ROWS;

    static constexpr bool isSquare = ROWS == COLUMNS;
    static constexpr bool isComplex = is_complex_v<T>;

    typedef T value_type;

    T data[COLUMNS][ROWS] = {};

    Matrix() = default;

    constexpr Matrix(std::initializer_list<std::initializer_list<T>> initializerList);

    Matrix(const Matrix<COLUMNS, ROWS, T>& other);

    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix(const Matrix<COLUMNS, ROWS, OTHER_T>& other);

    // m = m
    Matrix<COLUMNS, ROWS, T>& operator=(const Matrix<COLUMNS, ROWS, T>& other);

    // m == m
    bool equals(const Matrix<COLUMNS, ROWS, T>& other) const;
    bool operator==(const Matrix<COLUMNS, ROWS, T>& other) const;

    // m + m
    Matrix<COLUMNS, ROWS, T> add(const Matrix<COLUMNS, ROWS, T>& other) const;
    Matrix<COLUMNS, ROWS, T> operator+(const Matrix<COLUMNS, ROWS, T>& other) const;

    // m - m
    Matrix<COLUMNS, ROWS, T> subtract(const Matrix<COLUMNS, ROWS, T>& other) const;
    Matrix<COLUMNS, ROWS, T> operator-(const Matrix<COLUMNS, ROWS, T>& other) const;

    // m * m
    template<int OTHER_COLUMNS>
    Matrix<OTHER_COLUMNS, ROWS, T> multiply(const Matrix<OTHER_COLUMNS, COLUMNS, T>& other) const;
    template<int OTHER_COLUMNS>
    Matrix<OTHER_COLUMNS, ROWS, T> operator*(const Matrix<OTHER_COLUMNS, COLUMNS, T>& other) const;

    // m * v
    Vector<COLUMNS, T> multiply(const Vector<COLUMNS, T>& other) const;
    Vector<COLUMNS, T> operator*(const Vector<COLUMNS, T>& other) const;

    // m * #
    Matrix<COLUMNS, ROWS, T> multiply(T val) const;
    Matrix<COLUMNS, ROWS, T> operator*(T val) const;

    // m / #
    Matrix<COLUMNS, ROWS, T> divide(T scalar) const;
    Matrix<COLUMNS, ROWS, T> operator/(T scalar) const;

    // m += m
    Matrix<COLUMNS, ROWS, T>& addEquals(const Matrix<COLUMNS, ROWS, T>& other);
    Matrix<COLUMNS, ROWS, T>& operator+=(const Matrix<COLUMNS, ROWS, T>& other);

    // m -= m
    Matrix<COLUMNS, ROWS, T>& subtractEquals(const Matrix<COLUMNS, ROWS, T>& other);
    Matrix<COLUMNS, ROWS, T> operator-=(const Matrix<COLUMNS, ROWS, T>& other);

    // m *= #
    Matrix<COLUMNS, ROWS, T>& multiplyEquals(T val);
    Matrix<COLUMNS, ROWS, T>& operator*=(T val);

    // m /= #
    Matrix<COLUMNS, ROWS, T>& divideEquals(T scalar);
    Matrix<COLUMNS, ROWS, T>& operator/=(T scalar);

    // m = m
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T>& operator=(const Matrix<COLUMNS, ROWS, OTHER_T>& other);

    // m == m
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    bool equals(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const;
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    bool operator==(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const;

    // m + m
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T> add(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const;
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T> operator+(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const;

    // m - m
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T> subtract(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const;
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T> operator-(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const;

    // m * m
    template<int OTHER_COLUMNS, typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<OTHER_COLUMNS, ROWS, T> multiply(const Matrix<OTHER_COLUMNS, COLUMNS, OTHER_T>& other) const;
    template<int OTHER_COLUMNS, typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<OTHER_COLUMNS, ROWS, T> operator*(const Matrix<OTHER_COLUMNS, COLUMNS, OTHER_T>& other) const;

    // m * v
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Vector<COLUMNS, T> multiply(const Vector<COLUMNS, OTHER_T>& other) const;
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Vector<COLUMNS, T> operator*(const Vector<COLUMNS, OTHER_T>& other) const;

    // m * #
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T> multiply(OTHER_T val) const;
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T> operator*(OTHER_T val) const;

    // m / #
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T> divide(OTHER_T scalar) const;
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T> operator/(OTHER_T scalar) const;

    // m += m
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T>& addEquals(const Matrix<COLUMNS, ROWS, OTHER_T>& other);
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T>& operator+=(const Matrix<COLUMNS, ROWS, OTHER_T>& other);

    // m -= m
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T>& subtractEquals(const Matrix<COLUMNS, ROWS, OTHER_T>& other);
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T> operator-=(const Matrix<COLUMNS, ROWS, OTHER_T>& other);

    // m *= #
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T>& multiplyEquals(OTHER_T val);
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T>& operator*=(OTHER_T val);

    // m /= #
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T>& divideEquals(OTHER_T scalar);
    template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Matrix<COLUMNS, ROWS, T>& operator/=(OTHER_T scalar);

    template<int N>
    Vector<N, T> applyHomogeneousTransformation(const Vector<N, T>& point) const requires (isSquare);
    template<int N, typename OTHER_T> requires std::convertible_to<OTHER_T, T>
    Vector<N, T> applyHomogeneousTransformation(const Vector<N, OTHER_T>& point) const requires (isSquare);

    T* operator[](int index);
    const T* operator[](int index) const;

    Matrix<COLUMNS, ROWS, T> operator-() const;

    Matrix<ROWS, COLUMNS, T> transpose() const;
    Matrix<ROWS, COLUMNS, T> conjugateTranspose() const;

    Matrix<COLUMNS, ROWS, T> inverse() const requires (isSquare);

    enum DeterminantAlgorithm {
        laplace,
        triangular,
        hessenberg
    };

    T determinant(DeterminantAlgorithm algorithm = laplace) const requires (isSquare);

private:
    T laplaceDeterminant() const requires (isSquare);
    T triangularDeterminant() const requires (isSquare);

public:
    static Matrix<COLUMNS, ROWS, T> scalingMatrix(const Vector<COLUMNS, T>& factors) requires (isSquare);
    static Matrix<COLUMNS, ROWS, T> shearMatrix(int i, int j, T k) requires (isSquare);
    static Matrix<COLUMNS, ROWS, T> squeezeMatrix(int i, int j, T k) requires (isSquare);

    static Matrix<COLUMNS, ROWS, T> rotationMatrixAboutOrigin(T rot, RotationType rotationType = RotationType::radians) requires (isSquare && COLUMNS == 2);
    static Matrix<COLUMNS + 1, ROWS + 1, T> rotationMatrixAboutPoint(const Vector<COLUMNS, T>& p, T rot, RotationType rotationType = RotationType::radians) requires (isSquare && COLUMNS == 2);
    static Matrix<COLUMNS, ROWS, T> rotationMatrixAroundAxisThroughOrigin(const Vector<COLUMNS, T>& axis, T rot, RotationType rotationType = RotationType::radians) requires (isSquare && COLUMNS == 3);
    static Matrix<COLUMNS + 1, ROWS + 1, T> rotationMatrixAroundAxisNotThroughOrigin(const Vector<COLUMNS, T>& axis, const Vector<COLUMNS, T>& point, T rot, RotationType rotationType = RotationType::radians) requires (isSquare && COLUMNS == 3);
    static Matrix<COLUMNS, ROWS, T> rotationMatrixInPlaneThroughOrigin(const Vector<COLUMNS, T>& v1, const Vector<COLUMNS, T>& v2, T rot, RotationType rotationType = RotationType::radians) requires (isSquare && COLUMNS >= 3);
    static Matrix<COLUMNS, ROWS, T> rotationMatrixInPLaneNotThroughOrigin(const Vector<COLUMNS, T>& v1, const Vector<COLUMNS, T>& v2, const Vector<COLUMNS, T>& point, T rot, RotationType rotationType = RotationType::radians) requires (isSquare && COLUMNS >= 3);

    static Matrix<COLUMNS, ROWS, T> reflectionMatrixAlongAxisThroughOrigin(const Vector<COLUMNS, T>& axis) requires (isSquare);
    static Matrix<COLUMNS + 1, ROWS + 1, T> reflectionMatrixAlongAxisNotThroughOrigin(const Vector<COLUMNS, T>& axis, const Vector<COLUMNS, T>& point) requires (isSquare);

    static Matrix<COLUMNS + 1, ROWS + 1, T> translationMatrix(const Vector<COLUMNS, T>& translation) requires (isSquare);

    static Matrix<COLUMNS, ROWS, T> orthoMatrix(T left, T right, T bottom, T top, T near, T far) requires (isSquare && COLUMNS == 4);

    std::string toString() const;
    std::string toLaTex() const;

    template<int NUM_COLUMNS_TO_REMOVE>
    Matrix<COLUMNS - NUM_COLUMNS_TO_REMOVE, ROWS, T> removeColumns(const std::array<int, NUM_COLUMNS_TO_REMOVE>& columnsToRemove) const;
    template<int NUM_ROWS_TO_REMOVE>
    Matrix<COLUMNS, ROWS - NUM_ROWS_TO_REMOVE, T> removeRows(const std::array<int, NUM_ROWS_TO_REMOVE>& rowsToRemove) const;
    template<int NUM_COLUMNS_TO_REMOVE, int NUM_ROWS_TO_REMOVE>
    Matrix<COLUMNS - NUM_COLUMNS_TO_REMOVE, ROWS - NUM_ROWS_TO_REMOVE, T> removeColumnsAndRows(const std::array<int, NUM_COLUMNS_TO_REMOVE>& columnsToRemove, const std::array<int, NUM_ROWS_TO_REMOVE>& rowsToRemove) const;

    Matrix<COLUMNS - 1, ROWS, T> removeColumn(int columnToRemove) const;
    Matrix<COLUMNS, ROWS - 1, T> removeRow(int rowToRemove) const;
    Matrix<COLUMNS - 1, ROWS - 1, T> removeColumnAndRow(int columnToRemove, int rowToRemove) const;

    Matrix<COLUMNS, ROWS, T> swapRows(int rowA, int rowB);
    Matrix<COLUMNS, ROWS, T> swapColumns(int columnA, int columnB);

    explicit operator const T*() const;

    explicit operator T*();

    static Matrix<COLUMNS, ROWS, T> identity() requires (isSquare);

    bool isRowEchelon(bool pivotMustBeOne = false) const;
    Matrix<COLUMNS, ROWS, T> toRowEchelon() const;

    bool isReducedRowEchelon() const;
    Matrix<COLUMNS, ROWS, T> toReducedRowEchelon() const;

    int rank() const;

    bool isSymmetrical() const requires (isSquare);
    bool isSkewSymmetrical() const requires (isSquare);

    bool isHermitian() const requires (isSquare);
    bool isSkewHermitian() const requires (isSquare);

    bool isPositiveDefinite() const;
    bool isPositiveSemiDefinite() const;
    bool isNegativeDefinite() const;
    bool isNegativeSemiDefinite() const;

    Vector<ROWS, T> getColumnVector(int i) const;
    std::array<Vector<ROWS>, COLUMNS> getColumnVectors() const;
    Vector<COLUMNS, T> getRowVector(int i) const;
    std::array<Vector<COLUMNS>, ROWS> getRowVectors() const;

    void setColumnVectors(const std::array<Vector<ROWS, T>, COLUMNS>& columnVectors);
    void setColumnVector(int i, const Vector<ROWS, T>& v);
    void setRowVectors(std::array<Vector<COLUMNS, T>, ROWS> rowVectors);
    void setRowVector(int i, const Vector<COLUMNS, T>& v);

    T trace() const requires (isSquare);

    bool isUnitary() const requires (isSquare);
    bool isSpecialUnitary() const requires (isSquare);

    bool isOrthogonal() const requires (!isComplex && isSquare);
    bool isSpecialOrthogonal() const requires (!isComplex && isSquare);

    bool isSemiOrthogonal() const requires (!isComplex && !isSquare);

    bool isUpperTriangleMatrix() const requires (isSquare);
    bool isLowerTriangleMatrix() const requires (isSquare);

    bool isDiagonalMatrix() const requires (isSquare);

    bool isUpperUnitriangularMatrix() const requires (isSquare);
    bool isLowerUnitriangularMatrix() const requires (isSquare);

    bool isStrictlyUpperTriangularMatrix() const requires (isSquare);
    bool isStrictlyLowerTriangularMatrix() const requires (isSquare);

    bool isFrobeniusMatrix() const requires (isSquare);

    template<typename L_TYPE, typename U_TYPE, typename P_TYPE>
    struct LUPDecomposition {
        L_TYPE l;
        U_TYPE u;
        P_TYPE p;
    };

    LUPDecomposition<Matrix<std::min(ROWS, COLUMNS), ROWS, T>, Matrix<COLUMNS, std::min(ROWS, COLUMNS), T>, Matrix<ROWS, ROWS, T>> lupDecomposition() const;

    template<typename L_TYPE, typename U_TYPE, typename P_TYPE, typename Q_TYPE>
    struct LUPQDecomposition {
        L_TYPE l;
        U_TYPE u;
        P_TYPE p;
        Q_TYPE q;
    };

    LUPQDecomposition<Matrix<std::min(ROWS, COLUMNS), ROWS, T>, Matrix<COLUMNS, std::min(ROWS, COLUMNS), T>, Matrix<ROWS, ROWS, T>, Matrix<COLUMNS, COLUMNS, T>> lupqDecomposition() const;

    template<typename L_TYPE, typename U_TYPE>
    struct LUDecomposition {
        L_TYPE l;
        U_TYPE u;
    };

    LUDecomposition<Matrix<std::min(ROWS, COLUMNS), ROWS, T>, Matrix<COLUMNS, std::min(ROWS, COLUMNS), T>> luDecomposition() const;

    template<typename L_TYPE, typename D_TYPE, typename U_TYPE>
    struct LDUDecomposition {
        L_TYPE l;
        D_TYPE d;
        U_TYPE u;
    };

    LDUDecomposition<Matrix<std::min(ROWS, COLUMNS), ROWS, T>, Matrix<std::min(ROWS, COLUMNS), std::min(ROWS, COLUMNS), T>, Matrix<COLUMNS, std::min(ROWS, COLUMNS), T>> lduDecomposition() const;

    template<typename L_TYPE, typename L_TRANSPOSE_TYPE>
    struct CholeskyDecomposition {
        L_TYPE l;
        L_TRANSPOSE_TYPE lTranspose;
    };

    CholeskyDecomposition<Matrix<COLUMNS, ROWS, T>, Matrix<ROWS, COLUMNS, T>> choleskyDecomposition(bool allowPositiveSemiDefinite = false) const requires (isSquare);

    template<typename L_TYPE, typename D_TYPE, typename L_TRANSPOSE_TYPE>
    struct LDLDecomposition {
        L_TYPE l;
        D_TYPE d;
        L_TRANSPOSE_TYPE lTranspose;
    };

    LDLDecomposition<Matrix<COLUMNS, ROWS, T>, Matrix<COLUMNS, ROWS, T>, Matrix<ROWS, COLUMNS, T>> ldlDecomposition() const requires (isSquare);

    template<typename Q_TYPE, typename R_TYPE>
    struct QRDecomposition {
        Q_TYPE q;
        R_TYPE r;
    };

    QRDecomposition<Matrix<ROWS, ROWS, T>, Matrix<COLUMNS, ROWS, T>> qrDecomposition() const requires (isSquare);

    T minorOfElement(int c, int r) const requires (isSquare);

    Matrix<COLUMNS, ROWS, T> minorMatrix() const requires (isSquare);

    T cofactorOfElement(int c, int r) const requires (isSquare);

    Matrix<COLUMNS, ROWS, T> cofactorMatrix() const requires (isSquare);

    Matrix<ROWS, COLUMNS, T> adjoint() const requires (isSquare);

    bool isUpperHessenberg() const requires (isSquare);
    bool isUnreducedUpperHessenberg() const requires (isSquare);

    bool isLowerHessenberg() const requires (isSquare);
    bool isUnreducedLowerHessenberg() const requires (isSquare);

    bool isTridiagonal() const requires (isSquare);

    template<typename T_TYPE, typename Q_TYPE>
    struct LanczosAlgorithm {
        T_TYPE t;
        Q_TYPE q;
    };

    template<int ITER>
    LanczosAlgorithm<Matrix<ITER, ITER, T>, Matrix<ITER + 1, COLUMNS, T>> lanczosAlgorithm() const requires (isSquare);

    /**
     * Used to get an eigen-value approximation from an eigen-vector approximation
     * @param vec eigen-vector approximation
     * @return corresponding eigen-value approximation for given vector @a vec
     */
    T rayleighQuotient(const Vector<COLUMNS, T>& vec) const;

    template<typename EIGENVECTOR_TYPE, typename EIGENVALUE_TYPE>
    struct PowerIteration {
        EIGENVECTOR_TYPE vector;
        EIGENVALUE_TYPE value;
    };

    PowerIteration<Vector<COLUMNS, T>, T> powerIteration(int maxIterations, T tolerance = 1e-12) const;
};

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> multiply(const T lhs, const Matrix<COLUMNS, ROWS, T>& rhs) {
    return rhs.multiply(lhs);
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> operator*(const T lhs, const Matrix<COLUMNS, ROWS, T>& rhs) {
    return multiply(lhs, rhs);
}

template<int COLUMNS, int ROWS, is_scalar_v T, typename OTHER_T> requires std::convertible_to<OTHER_T, T>
Matrix<COLUMNS, ROWS, T> multiply(const OTHER_T lhs, const Matrix<COLUMNS, ROWS, T>& rhs) {
    return rhs.multiply(lhs);
}

template<int COLUMNS, int ROWS, is_scalar_v T, typename OTHER_T> requires std::convertible_to<OTHER_T, T>
Matrix<COLUMNS, ROWS, T> operator*(const OTHER_T lhs, const Matrix<COLUMNS, ROWS, T>& rhs) {
    return multiply(lhs, rhs);
}

template<int N, is_scalar_v T>
Matrix<N, N, T> Vector<N, T>::crossProductMatrix() const requires (N == 3) {
    return {{0, -data[2], data[1]}, {data[2], 0, -data[0]}, {-data[1], data[0], 0}};
}

template<int N, is_scalar_v T>
template<int OTHER_N>
Matrix<OTHER_N, N, T> Vector<N, T>::outerProductMatrix(const Vector<OTHER_N, T>& v) const {
    Matrix<OTHER_N, N, T> result;

    for (int c = 0; c < result.columns; c++) {
        for (int r = 0; r < result.rows; r++) {
            if constexpr (isComplex)
                result[c][r] = data[r] * std::conj(v[c]);
            else
                result[c][r] = data[r] * v[c];
        }
    }

    return result;
}

template<int N, is_scalar_v T>
template<int OTHER_N, typename OTHER_T> requires std::convertible_to<OTHER_T, T>
Matrix<OTHER_N, N, T> Vector<N, T>::outerProductMatrix(const Vector<OTHER_N, OTHER_T>& v) const {
    Matrix<OTHER_N, N, T> result;

    for (int c = 0; c < result.columns; c++) {
        for (int r = 0; r < result.rows; r++) {
            if constexpr (isComplex)
                result[c][r] = data[r] * std::conj(v[c]);
            else
                result[c][r] = data[r] * v[c];
        }
    }

    return result;
}

template<int N, is_scalar_v T>
template<int COLUMNS>
Vector<COLUMNS, T> Vector<N, T>::multiply(const Matrix<COLUMNS, N, T>& m) const {
    Vector<COLUMNS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < N; r++) {
            result[c] += data[r] * m[c][r];
        }
    }

    return result;
}

template<int N, is_scalar_v T>
template<int COLUMNS>
Vector<COLUMNS, T> Vector<N, T>::operator*(const Matrix<COLUMNS, N, T>& m) const {
    return multiply(m);
}

template<int N, is_scalar_v T>
template<int COLUMNS, typename OTHER_T> requires std::convertible_to<OTHER_T, T>
Vector<COLUMNS, T> Vector<N, T>::multiply(const Matrix<COLUMNS, N, OTHER_T>& m) const {
    Vector<COLUMNS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < N; r++) {
            result[c] += data[r] * m[r][c];
        }
    }

    return result;
}

template<int N, is_scalar_v T>
template<int COLUMNS, typename OTHER_T> requires std::convertible_to<OTHER_T, T>
Vector<COLUMNS, T> Vector<N, T>::operator*(const Matrix<COLUMNS, N, OTHER_T>& m) const {
    return multiply(m);
}

// block matrix
template<int COLUMNS, int ROWS, int B_COLUMNS, int B_ROWS, is_scalar_v B_T>
struct Matrix<COLUMNS, ROWS, Matrix<B_COLUMNS, B_ROWS, B_T>> {
    static constexpr int columns = COLUMNS;
    static constexpr int rows = ROWS;
    static constexpr int bColumns = B_COLUMNS;
    static constexpr int bRows = B_ROWS;

    static constexpr bool isSquare = ROWS == COLUMNS;
    static constexpr bool bIsSquare = B_ROWS == B_COLUMNS;

    Matrix<B_COLUMNS, B_ROWS, B_T> data[COLUMNS][ROWS] = {};
};
