#pragma once

#include <cassert>
#include <complex>
#include <regex>
#include <sstream>
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

    static constexpr T epsilon = ::epsilon<T>();

    T data[COLUMNS][ROWS] = {};

    Matrix() = default;

    constexpr Matrix(std::initializer_list<std::initializer_list<T>> initializerList);

    Matrix(const Matrix<COLUMNS, ROWS, T>& other);

    template<is_convertable_to<T> OTHER_T>
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
    Matrix<COLUMNS, ROWS, T>& operator*=(T val) const;

    // m /= #
    Matrix<COLUMNS, ROWS, T>& divideEquals(T scalar);
    Matrix<COLUMNS, ROWS, T>& operator/=(T scalar);

    // m = m
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& operator=(const Matrix<COLUMNS, ROWS, OTHER_T>& other);

    // m == m
    template<is_convertable_to<T> OTHER_T>
    bool equals(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const;
    template<is_convertable_to<T> OTHER_T>
    bool operator==(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const;

    // m + m
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> add(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const;
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> operator+(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const;

    // m - m
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> subtract(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const;
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> operator-(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const;

    // m * m
    template<int OTHER_COLUMNS, is_convertable_to<T> OTHER_T>
    Matrix<OTHER_COLUMNS, ROWS, T> multiply(const Matrix<OTHER_COLUMNS, COLUMNS, OTHER_T>& other) const;
    template<int OTHER_COLUMNS, is_convertable_to<T> OTHER_T>
    Matrix<OTHER_COLUMNS, ROWS, T> operator*(const Matrix<OTHER_COLUMNS, COLUMNS, OTHER_T>& other) const;

    // m * v
    template<is_convertable_to<T> OTHER_T>
    Vector<COLUMNS, T> multiply(const Vector<COLUMNS, OTHER_T>& other) const;
    template<is_convertable_to<T> OTHER_T>
    Vector<COLUMNS, T> operator*(const Vector<COLUMNS, OTHER_T>& other) const;

    // m * #
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> multiply(OTHER_T val) const;
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> operator*(OTHER_T val) const;

    // m / #
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> divide(OTHER_T scalar) const;
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> operator/(OTHER_T scalar) const;

    // m += m
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& addEquals(const Matrix<COLUMNS, ROWS, OTHER_T>& other);
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& operator+=(const Matrix<COLUMNS, ROWS, OTHER_T>& other);

    // m -= m
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& subtractEquals(const Matrix<COLUMNS, ROWS, OTHER_T>& other);
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> operator-=(const Matrix<COLUMNS, ROWS, OTHER_T>& other);

    // m *= #
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& multiplyEquals(OTHER_T val);
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& operator*=(OTHER_T val) const;

    // m /= #
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& divideEquals(OTHER_T scalar);
    template<is_convertable_to<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& operator/=(OTHER_T scalar);

    template<int N>
    Vector<N, T> applyHomogeneousTransformation(const Vector<N, T>& point) const requires (isSquare);

    template<int N, is_convertable_to<T> OTHER_T>
    Vector<N, T> applyHomogeneousTransformation(const Vector<N, OTHER_T>& point) const requires (isSquare);

    T* operator[](int index);
    const T* operator[](int index) const;

    Matrix<ROWS, COLUMNS, T> transpose() const;

    Matrix<ROWS, COLUMNS, T> conjugateTranspose() const;

    Matrix<COLUMNS, ROWS, T> inverse() const requires (isSquare);

    enum DeterminantAlgorithm {
        laplace,
        triangular,
        hessenberg
    };

    T determinant(const DeterminantAlgorithm algorithm = laplace) const requires (isSquare);

private:
    T laplaceDeterminant() const requires (isSquare);

    T triangularDeterminant() const requires (isSquare);

public:
    static Matrix<COLUMNS, ROWS, T> scalingMatrix(const Vector<COLUMNS, T>& factors) requires (isSquare);

    static Matrix<COLUMNS, ROWS, T> shearMatrix(const int i, const int j, const T k) requires (isSquare);

    static Matrix<COLUMNS, ROWS, T> squeezeMatrix(const int i, const int j, const T k) requires (isSquare);

    static Matrix<COLUMNS, ROWS, T> rotationMatrixAboutOrigin(const T rot, const RotationType rotationType = RotationType::radians) requires (isSquare && COLUMNS == 2);

    static Matrix<COLUMNS + 1, ROWS + 1, T> rotationMatrixAboutPoint(const Vector<COLUMNS, T>& p, const T rot, const RotationType rotationType = RotationType::radians) requires (isSquare && COLUMNS == 2);

    static Matrix<COLUMNS, ROWS, T> rotationMatrixAroundAxisThroughOrigin(const Vector<COLUMNS, T>& axis, const T rot, const RotationType rotationType = RotationType::radians) requires (isSquare && COLUMNS == 3);

    static Matrix<COLUMNS + 1, ROWS + 1, T> rotationMatrixAroundAxisNotThroughOrigin(const Vector<COLUMNS, T>& axis, const Vector<COLUMNS, T>& point, const T rot, const RotationType rotationType = RotationType::radians) requires (isSquare && COLUMNS == 3);

    static Matrix<COLUMNS, ROWS, T> rotationMatrixInPlaneThroughOrigin(const Vector<COLUMNS, T>& v1, const Vector<COLUMNS, T>& v2, const T rot, const RotationType rotationType = RotationType::radians) requires (isSquare && COLUMNS >= 3);

    static Matrix<COLUMNS, ROWS, T> rotationMatrixInPLaneNotThroughOrigin(const Vector<COLUMNS, T>& v1, const Vector<COLUMNS, T>& v2, const Vector<COLUMNS, T>& point, const T rot, const RotationType rotationType = RotationType::radians) requires (isSquare && COLUMNS >= 3);

    static Matrix<COLUMNS, ROWS, T> reflectionMatrixAlongAxisThroughOrigin(const Vector<COLUMNS, T>& axis) requires (isSquare);

    static Matrix<COLUMNS + 1, ROWS + 1, T> reflectionMatrixAlongAxisNotThroughOrigin(const Vector<COLUMNS, T>& axis, const Vector<COLUMNS, T>& point) requires (isSquare);

    static Matrix<COLUMNS + 1, ROWS + 1, T> translationMatrix(const Vector<COLUMNS, T>& translation) requires (isSquare);

    static Matrix<COLUMNS, ROWS, T> orthoMatrix(const T left, const T right, const T bottom, const T top, const T near, const T far) requires (isSquare && COLUMNS == 4);

    std::string toString() const;
    std::string toLaTex() const;

    template<int NUM_COLUMNS_TO_REMOVE>
    Matrix<COLUMNS - NUM_COLUMNS_TO_REMOVE, ROWS, T> removeColumns(const std::array<int, NUM_COLUMNS_TO_REMOVE>& columnsToRemove) const;

    template<int NUM_ROWS_TO_REMOVE>
    Matrix<COLUMNS, ROWS - NUM_ROWS_TO_REMOVE, T> removeRows(const std::array<int, NUM_ROWS_TO_REMOVE>& rowsToRemove) const;

    template<int NUM_COLUMNS_TO_REMOVE, int NUM_ROWS_TO_REMOVE>
    Matrix<COLUMNS - NUM_COLUMNS_TO_REMOVE, ROWS - NUM_ROWS_TO_REMOVE, T> removeColumnsAndRows(const std::array<int, NUM_COLUMNS_TO_REMOVE>& columnsToRemove, const std::array<int, NUM_ROWS_TO_REMOVE>& rowsToRemove) const;

    Matrix<COLUMNS, ROWS, T> swapRows(const int rowA, const int rowB);

    Matrix<COLUMNS, ROWS, T> swapColumns(const int columnA, const int columnB);

    explicit operator const T*() const;

    explicit operator T*();

    static Matrix<COLUMNS, ROWS, T> identity() requires (isSquare);

    bool isRowEchelon(bool pivotMustBeOne = false) const {
        int lastPivotColumn = -1;

        for (int r = 0; r < ROWS; r++) {
            bool foundNonZero = false;

            for (int c = 0; c < COLUMNS; c++) {
                // this is a pivot
                if (!foundNonZero && !compare(data[c][r], 0)) {
                    // this is to the left of the last pivot
                    if (c < lastPivotColumn)
                        return false;

                    // the pivot must be one and it isn't
                    if (!compare(data[c][r], 1) && pivotMustBeOne)
                        return false;

                    lastPivotColumn = c;
                }

                if (!compare(data[c][r], 0))
                    foundNonZero = true;
            }

            // this entire row was zeros, and this wasn't the last row
            if (!foundNonZero || r != ROWS - 1) {
                return false;
            }
        }

        return true;
    }

    Matrix<COLUMNS, ROWS, T> toRowEchelon() const {
        Matrix<COLUMNS, ROWS, T> m = *this;

        for (int c = 0; c < std::min(ROWS, COLUMNS); c++) {
            // handle row swaps
            {
                int rowIndex = -1;
                const T pivot = m[c][c];
                T rowValue = std::abs(pivot);

                // iterate through rows of this column. Looking for the biggest boi
                for (int r = c + 1; r < ROWS; r++) {
                    T curValue = std::abs(m[c][r]);

                    if (curValue > rowValue) {
                        rowIndex = r;
                        rowValue = curValue;
                    }
                }

                // we found nothing so we skip this column
                if (rowIndex == -1 && compare(pivot, 0)) {
                    continue;
                }

                if (rowIndex != -1 && rowIndex != c) {
                    // swap u and p rows like normal
                    m.swapRows(c, rowIndex);
                }
            }

            const T pivot = m[c][c];

            // iterate through things beneath that pivot in the matrix
            for (int r = c + 1; r < ROWS; r++) {
                T val = m[c][r];

                T multiplierToPivotRow = val / pivot;

                // do this row minus other row times multiplier
                for (int i = c; i < COLUMNS; i++) {
                    m[i][r] += -multiplierToPivotRow * m[i][c];
                }
            }
        }

        return m;
    }

    bool isReducedRowEchelon() const {
        int lastPivotColumn = -1;

        for (int r = 0; r < ROWS; r++) {
            bool foundNonZero = false;

            for (int c = 0; c < COLUMNS; c++) {
                // this is a pivot
                if (!foundNonZero && !compare(data[c][r], 0)) {
                    // this is to the left of the last pivot
                    if (c < lastPivotColumn)
                        return false;

                    // the pivot isn't 1
                    if (!compare(data[c][r], 1))
                        return false;

                    // check that no other number in that column is a nonzero value
                    for (int i = 0; i < ROWS; i++) {
                        if (!compare(data[c][i], 0) && !compare(data[c][i], data[c][r]))
                            return false;
                    }

                    lastPivotColumn = c;
                }

                if (!compare(data[c][r], 0))
                    foundNonZero = true;
            }

            // this entire row was zeros, and this wasn't the last row
            if (!foundNonZero || r != ROWS - 1) {
                return false;
            }
        }

        return true;
    }

    Matrix<COLUMNS, ROWS, T> toReducedRowEchelon() const {
        Matrix<COLUMNS, ROWS, T> m = *this;

        for (int c = 0; c < std::min(ROWS, COLUMNS); c++) {
            // handle row swaps
            {
                int rowIndex = -1;
                const T pivot = m[c][c];
                T rowValue = std::abs(pivot);

                // iterate through rows of this column. Looking for the biggest boi
                for (int r = c + 1; r < ROWS; r++) {
                    T curValue = std::abs(m[c][r]);

                    if (curValue > rowValue) {
                        rowIndex = r;
                        rowValue = curValue;
                    }
                }

                // we found nothing so we skip this column
                if (rowIndex == -1 && compare(pivot, 0)) {
                    continue;
                }

                if (rowIndex != -1 && rowIndex != c) {
                    // swap u and p rows like normal
                    m.swapRows(c, rowIndex);
                }
            }

            {
                // normalize pivot row to 1
                T value = m[c][c];
                for (int cc = c + 1; cc < COLUMNS; cc++) {
                    m[cc][c] /= value;
                }
                m[c][c] = 1;
            }

            const T pivot = m[c][c];

            // iterate through things beneath that pivot in the matrix
            for (int r = c + 1; r < ROWS; r++) {
                T val = m[c][r];

                T multiplierToPivotRow = val / pivot;

                // do this row minus other row times multiplier
                for (int i = c; i < COLUMNS; i++) {
                    m[i][r] += -multiplierToPivotRow * m[i][c];
                }
            }
        }

        return m;
    }

    int rank() const {
        Matrix<COLUMNS, ROWS, T> ref = toRowEchelon();
        int result = 0;

        for (int r = 0; r < ROWS; r++) {
            if (!compare(ref[0][r], 0)) {
                result++;
            }
            else {
                bool nonZero = false;

                for (int c = 0; c < COLUMNS; c++) {
                    if (!compare(ref[c][r], 0)) {
                        nonZero = true;
                        break;
                    }
                }

                if (nonZero)
                    result++;
            }
        }

        return result;
    }

    bool isSymmetrical() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                if (!compare(data[c][r], data[r][c]))
                    return false;
            }
        }

        return true;
    }

    bool isSkewSymmetrical() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                if (!compare(data[c][r], -data[r][c]))
                    return false;
            }
        }

        return true;
    }

    bool isHermitian() const requires (isSquare) {
        if constexpr (!isComplex)
            return isSymmetrical();

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                if (!compare(data[c][r], std::conj(data[r][c]))) {
                    return false;
                }
            }
        }

        return true;
    }

    bool isSkewHermitian() const requires (isSquare) {
        if constexpr (!isComplex)
            return isSkewSymmetrical();

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                if (!compare(data[c][r], -std::conj(data[r][c])))
                    return false;
            }
        }

        return true;
    }

    bool isPositiveDefinite() const requires (isSquare) {
        // TODO:
    }

    bool isPositiveSemiDefinite() const requires (isSquare) {
        // TODO:
    }

    bool isNegativeDefinite() const requires (isSquare) {
        // TODO:
    }

    bool isNegativeSemiDefinite() const requires (isSquare) {
        // TODO:
    }

    Vector<ROWS, T> getColumnVector(const int i) const {
        Vector<ROWS, T> v;

        for (int j = 0; j < ROWS; j++) {
            v[j] = data[i][j];
        }

        return v;
    }

    std::array<Vector<ROWS>, COLUMNS> getColumnVectors() const {
        std::array<Vector<ROWS>, COLUMNS> vecs;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                vecs[c][r] = data[c][r];
            }
        }

        return vecs;
    }

    Vector<COLUMNS, T> getRowVector(const int i) const {
        Vector<COLUMNS, T> v;

        for (int j = 0; j < COLUMNS; j++) {
            v[j] = data[j][i];
        }

        return v;
    }

    std::array<Vector<COLUMNS>, ROWS> getRowVectors() const {
        std::array<Vector<COLUMNS>, ROWS> vecs;

        for (int r = 0; r < ROWS; r++) {
            for (int c = 0; c < COLUMNS; c++) {
                vecs[r][c] = data[c][r];
            }
        }

        return vecs;
    }

    void setColumnVectors(const std::array<Vector<ROWS, T>, COLUMNS>& columnVectors) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                data[c][r] = columnVectors[c][r];
            }
        }
    }

    void setColumnVector(const int i, const Vector<ROWS, T>& v) {
        for (int j = 0; j < ROWS; j++) {
            data[i][j] = v[j];
        }
    }

    void setRowVectors(const std::array<Vector<COLUMNS, T>, ROWS> rowVectors) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                data[c][r] = rowVectors[r][c];
            }
        }
    }

    void setRowVector(const int i, const Vector<COLUMNS, T>& v) {
        for (int j = 0; j < COLUMNS; j++) {
            data[j][i] = v[j];
        }
    }

    T trace() const requires (isSquare) {
        T sum = {};

        for (int c = 0; c < COLUMNS; c++) {
            sum += data[c][c];
        }

        return sum;
    }

    bool isUnitary() const requires (isSquare) {
        // TODO:
    }

    bool isSpecialUnitary() const requires (isSquare) {
        // TODO:
    }

    bool isOrthogonal() const requires (!isComplex && isSquare) {
        // TODO:
    }

    bool isSpecialOrthogonal() const requires (!isComplex && isSquare) {
        // TODO:
    }

    bool isUpperTriangleMatrix() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            // dont worry about out of bounds, loop wont even run if c + 1 is too big
            for (int r = c + 1; r < ROWS; r++) {
                if (!compare(data[c][r], 0))
                    return false;
            }
        }

        return true;
    }

    bool isLowerTriangleMatrix() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < c; r++) {
                if (!compare(data[c][r], 0))
                    return false;
            }
        }

        return true;
    }

    bool isDiagonalMatrix() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                if (c == r)
                    continue;

                if (!compare(data[c][r], 0))
                    return false;
            }
        }

        return true;
    }

    bool isUpperUnitriangularMatrix() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = c; r < ROWS; r++) {
                if (c == r) {
                    if (!compare(data[c][r], 1)) {
                        return false;
                    }
                }
                else if (!compare(data[c][r], 0)) {
                    return false;
                }
            }
        }

        return true;
    }

    bool isLowerUnitriangularMatrix() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r <= c; r++) {
                if (c == r) {
                    if (!compare(data[c][r], 1)) {
                        return false;
                    }
                }
                else if (!compare(data[c][r], 0))
                    return false;
            }
        }

        return true;
    }

    bool isStrictlyUpperTriangularMatrix() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = c; r < ROWS; r++) {
                if (!compare(data[c][r], 0)) {
                    return false;
                }
            }
        }

        return true;
    }

    bool isStrictlyLowerTriangularMatrix() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r <= c; r++) {
                if (!compare(data[c][r], 0))
                    return false;
            }
        }

        return true;
    }

    bool isFrobeniusMatrix() const requires (isSquare) {
        int columnWithNonZero = -1;
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                // not a one on diagonal
                if (c == r && !compare(data[c][r], 1)) {
                    return false;
                }
                if (!compare(data[c][r], 0)) {
                    // its above
                    if (r <= c)
                        return false;

                    // too many columns with nonzero
                    if (columnWithNonZero != c)
                        return false;

                    columnWithNonZero = c;
                }
            }
        }

        return true;
    }

#pragma region Decomposiitons
    template<typename L_TYPE, typename U_TYPE, typename P_TYPE>
    struct LUPDecomposition {
        L_TYPE l;
        U_TYPE u;
        P_TYPE p;
    };

    LUPDecomposition<Matrix<std::min(ROWS, COLUMNS), ROWS, T>, Matrix<COLUMNS, std::min(ROWS, COLUMNS), T>, Matrix<ROWS, ROWS, T>> lupDecomposition() const {
        Matrix<std::min(ROWS, COLUMNS), ROWS, T> l = Matrix<ROWS, ROWS, T>::identity();
        Matrix<COLUMNS, std::min(ROWS, COLUMNS), T> u = *this;
        Matrix<ROWS, ROWS, T> p = Matrix<ROWS, ROWS, T>::identity();

        // use std::min so we never access the pivot that's out of the matrix
        for (int c = 0; c < std::min(ROWS, COLUMNS); c++) {
            // handle row swaps
            {
                int rowIndex = -1;
                const T pivot = u[c][c];
                T rowValue = std::abs(pivot);

                // iterate through rows of this column. Looking for the biggest boi
                for (int r = c + 1; r < ROWS; r++) {
                    T curValue = std::abs(u[c][r]);

                    if (curValue > rowValue) {
                        rowIndex = r;
                        rowValue = curValue;
                    }
                }

                // we found nothing but we needed to find something
                if (rowIndex == -1 && compare(pivot, 0)) {
                    throw std::runtime_error("Cannot LUP decompose singular matrix");
                }

                if (rowIndex != -1) {
                    // swap u and p rows like normal
                    u.swapRows(c, rowIndex);
                    p.swapRows(c, rowIndex);

                    // swap rows of l before column c
                    for (int i = 0; i < c; ++i) {
                        T temp = l[i][c];
                        l[i][c] = l[i][rowIndex];
                        l[i][rowIndex] = temp;
                    }
                }
            }

            T pivot = u[c][c];

            // iterate through things beneath that pivot in the matrix
            for (int r = c + 1; r < ROWS; r++) {
                T val = u[c][r];

                T multiplierToPivotRow = val / pivot;

                l[c][r] = multiplierToPivotRow;

                // do this row minus other row times multiplier
                for (int i = c; i < COLUMNS; i++) {
                    u[i][r] += -multiplierToPivotRow * u[i][c];
                }
            }
        }

        return {l, u, p};
    }

    template<typename L_TYPE, typename U_TYPE, typename P_TYPE, typename Q_TYPE>
    struct LUPQDecomposition {
        L_TYPE l;
        U_TYPE u;
        P_TYPE p;
        Q_TYPE q;
    };

    LUPQDecomposition<Matrix<std::min(ROWS, COLUMNS), ROWS, T>, Matrix<COLUMNS, std::min(ROWS, COLUMNS), T>, Matrix<ROWS, ROWS, T>, Matrix<COLUMNS, COLUMNS, T>> lupqDecomposition() const {
        Matrix<std::min(ROWS, COLUMNS), ROWS, T> l = Matrix<ROWS, ROWS, T>::identity();
        Matrix<COLUMNS, std::min(ROWS, COLUMNS), T> u = *this;
        Matrix<ROWS, ROWS, T> p = Matrix<ROWS, ROWS, T>::identity();
        Matrix<COLUMNS, COLUMNS, T> q = Matrix<COLUMNS, COLUMNS, T>::identity();

        // use std::min so we never access the pivot that's out of the matrix
        for (int c = 0; c < std::min(ROWS, COLUMNS); c++) {
            // handle row AND column swaps
            {
                int rowIndex = -1;
                int columnIndex = -1;

                const T pivot = u[c][c];

                T maxValue = std::abs(pivot);

                for (int i = c; i < COLUMNS; i++) {
                    for (int j = c; j < ROWS; j++) {
                        T curValue = std::abs(u[i][j]);

                        if (curValue > maxValue) {
                            maxValue = curValue;
                            columnIndex = i;
                            rowIndex = j;
                        }
                    }
                }

                if (rowIndex == -1 && columnIndex == -1 && compare(pivot, 0)) {
                    throw std::runtime_error("Cannot LUPQ decompose singular matrix");
                }

                if (rowIndex != c) {
                    u.swapRows(c, rowIndex);
                    p.swapRows(c, rowIndex);
                    // swap rows of L before column c
                    for (int i = 0; i < c; ++i) {
                        T temp = l[i][c];
                        l[i][c] = l[i][rowIndex];
                        l[i][rowIndex] = temp;
                    }
                }

                if (columnIndex != c) {
                    u.swapColumns(c, columnIndex);
                    q.swapColumns(c, columnIndex);
                }
            }

            T pivot = u[c][c];

            // iterate through things beneath that pivot in the matrix
            for (int r = c + 1; r < ROWS; r++) {
                T val = u[c][r];

                T multiplierToPivotRow = val / pivot;

                l[c][r] = multiplierToPivotRow;

                // do this row minus other row times multiplier
                for (int i = c; i < COLUMNS; i++) {
                    u[i][r] += -multiplierToPivotRow * u[i][c];
                }
            }
        }

        return {l, u, p, q};
    }

    template<typename L_TYPE, typename U_TYPE>
    struct LUDecomposition {
        L_TYPE l;
        U_TYPE u;
    };

    LUDecomposition<Matrix<std::min(ROWS, COLUMNS), ROWS, T>, Matrix<COLUMNS, std::min(ROWS, COLUMNS), T>> luDecomposition() const {
        Matrix<std::min(ROWS, COLUMNS), ROWS, T> l = Matrix<ROWS, ROWS, T>::identity();
        Matrix<COLUMNS, std::min(ROWS, COLUMNS), T> u = *this;

        // use std::min so we never access the pivot that's out of the matrix
        for (int c = 0; c < std::min(ROWS, COLUMNS); c++) {
            T pivot = u[c][c];

            if (compare(pivot, 0))
                throw std::runtime_error("Cannot LU decompose matrix due to zero pivot, try LUP or LUPQ");

            // iterate through things beneath that pivot in the matrix
            for (int r = c + 1; r < ROWS; r++) {
                T val = u[c][r];

                T multiplierToPivotRow = val / pivot;

                l[c][r] = multiplierToPivotRow;

                // do this row minus other row times multiplier
                for (int i = c; i < COLUMNS; i++) {
                    u[i][r] += -multiplierToPivotRow * u[i][c];
                }
            }
        }

        return {l, u};
    }

    template<typename L_TYPE, typename D_TYPE, typename U_TYPE>
    struct LDUDecomposition {
        L_TYPE l;
        D_TYPE d;
        U_TYPE u;
    };

    LDUDecomposition<Matrix<std::min(ROWS, COLUMNS), ROWS, T>, Matrix<std::min(ROWS, COLUMNS), std::min(ROWS, COLUMNS), T>, Matrix<COLUMNS, std::min(ROWS, COLUMNS), T>> lduDecomposition() const {
        Matrix<std::min(ROWS, COLUMNS), ROWS, T> l = Matrix<ROWS, ROWS, T>::identity();
        Matrix<std::min(ROWS, COLUMNS), std::min(ROWS, COLUMNS), T> d = Matrix<std::min(ROWS, COLUMNS), std::min(ROWS, COLUMNS), T>::identity();
        Matrix<COLUMNS, std::min(ROWS, COLUMNS), T> u = *this;

        // use std::min so we never access the pivot that's out of the matrix
        for (int c = 0; c < std::min(ROWS, COLUMNS); c++) {
            T pivot = u[c][c];

            if (compare(pivot, 0))
                throw std::runtime_error("Cannot LDU decompose matrix due to zero pivot");

            // iterate through things beneath that pivot in the matrix
            for (int r = c + 1; r < ROWS; r++) {
                T val = u[c][r];

                T multiplierToPivotRow = val / pivot;

                l[c][r] = multiplierToPivotRow;

                // do this row minus other row times multiplier
                for (int i = c; i < COLUMNS; i++) {
                    u[i][r] += -multiplierToPivotRow * u[i][c];
                }
            }
        }

        for (int c = 0; c < COLUMNS; c++) {
            d[c][c] = u[c][c];

            for (int r = 0; r < ROWS; r++) {
                u[c][r] /= d[c][c];
            }
        }

        return {l, d, u};
    }

    template<typename L_TYPE, typename L_TRANSPOSE_TYPE>
    struct CholeskyDecomposition {
        L_TYPE l;
        L_TRANSPOSE_TYPE lTranspose;
    };

    CholeskyDecomposition<Matrix<COLUMNS, ROWS, T>, Matrix<ROWS, COLUMNS, T>> choleskyDecomposition(const bool allowPositiveSemiDefinite = false) const requires (isSquare) {
        if (!isHermitian()) {
            throw std::runtime_error("Cannot find Cholesky Decomposition of non hermitian/symmetric matrix");
        }

        Matrix<COLUMNS, ROWS, T> l;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r <= c; r++) {
                if (r == c) {
                    T value = data[c][c];

                    for (int k = 0; k < c; k++) {
                        // TODO: Find a way to replace this if
                        if constexpr (isComplex) {
                            value -= l[k][c] * std::conj(l[k][c]);
                        }
                        else {
                            value -= std::pow(l[k][c], 2);
                        }
                    }

                    if (value < 0) {
                        throw std::runtime_error("Cannot cholesky decompose non positive definite matrix");
                    }

                    if (compare(value, 0) && !allowPositiveSemiDefinite) {
                        throw std::runtime_error("Cannot cholesky decompose non semi-positive definite matrix");
                    }

                    l[c][c] = std::sqrt(value);
                }
                else if (r > c) {
                    T value = data[c][r];

                    for (int k = 0; k < c; k++) {
                        value -= l[k][r] * l[k][c];
                    }

                    l[c][r] = value / l[c][c];
                }
            }
        }

        return {l, l.conjugateTranspose()};
    }

    template<typename L_TYPE, typename D_TYPE, typename L_TRANSPOSE_TYPE>
    struct LDLDecomposition {
        L_TYPE l;
        D_TYPE d;
        L_TRANSPOSE_TYPE lTranspose;
    };

    LDLDecomposition<Matrix<COLUMNS, ROWS, T>, Matrix<COLUMNS, ROWS, T>, Matrix<ROWS, COLUMNS, T>> ldlDecomposition() const requires (isSquare) {
        Matrix<COLUMNS, ROWS, T> l;
        Matrix<COLUMNS, ROWS, T> d;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r <= c; r++) {
                if constexpr (!isComplex) {
                    if (r == c) {
                        T value = data[c][c];

                        for (int k = 0; k < c; k++) {
                            value -= std::pow(l[k][c], 2) * d[k][k];
                        }

                        d[c][c] = value;
                    }
                    else if (r > c) {
                        T value = data[c][r];

                        for (int k = 0; k < c; k++) {
                            value -= l[k][r] * l[k][c] * d[k][k];
                        }

                        l[c][r] = (1 / d[c][c]) * value;
                    }
                }
                else {
                    if (r == c) {
                        T value = data[c][c];

                        for (int k = 0; k < c; k++) {
                            T complexConjugate = l[k][c];
                            complexConjugate.imag *= -1;
                            value -= l[k][c] * complexConjugate * d[k][k];
                        }

                        d[c][c] = value;
                    }
                    else if (r > c) {
                        T value = data[c][r];

                        for (int k = 0; k < c; k++) {
                            T complexConjugate = l[k][c];
                            complexConjugate.imag *= -1;

                            value -= l[k][r] * complexConjugate * d[k][k];
                        }

                        l[c][r] = (1 / d[c][c]) * value;
                    }
                }
            }
        }

        return {l, d, l.conjugateTranspose()};
    }

    template<typename Q_TYPE, typename R_TYPE>
    struct QRDecomposition {
        Q_TYPE q;
        R_TYPE r;
    };

    QRDecomposition<Matrix<ROWS, ROWS, T>, Matrix<COLUMNS, ROWS, T>> qrDecomposition() const requires (isSquare) {
        std::array<Vector<ROWS>, COLUMNS> a = getColumnVectors();
        std::array<Vector<ROWS>, COLUMNS> u = {};
        Matrix<ROWS, ROWS, T> q;

        for (int k = 0; k < COLUMNS; k++) {
            u[k] = a[k];

            for (int j = 0; j < k; j++) {
                u[k] -= u[j].projection(a[k]);
            }

            T uMagnitude = u[k].magnitude();

            for (int i = 0; i < ROWS; i++) {
                q[k][i] = u[k][i] / uMagnitude;
            }
        }

        Matrix<COLUMNS, ROWS, T> r = q.transpose() * *this;

        return {q, r};
    }

#pragma endregion

    T minorOfElement(const int c, const int r) const requires (isSquare) {
        return removeColumnsAndRows({c}, {r}).determinant();
    }

    Matrix<COLUMNS, ROWS, T> minorMatrix() const requires (isSquare) {
        Matrix<COLUMNS, ROWS, T> result;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                result[c][r] = minorOfElement(c, r);
            }
        }

        return result;
    }

    T cofactorOfElement(const int c, const int r) const requires (isSquare) {
        return minorOfElement(c, r) * std::pow(-1, c + r);
    }

    Matrix<COLUMNS, ROWS, T> cofactorMatrix() const requires (isSquare) {
        Matrix<COLUMNS, ROWS, T> result;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                result[c][r] = cofactorOfElement(c, r);
            }
        }

        return result;
    }

    Matrix<ROWS, COLUMNS, T> adjoint() const requires (isSquare) {
        return cofactorMatrix().transpose();
    }

    bool isUpperHessenberg() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = c + 2; r < ROWS; r++) {
                if (!compare(data[c][r], 0))
                    return false;
            }
        }

        return true;
    }

    bool isUnreducedUpperHessenberg() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = c + 1; r < ROWS; r++) {
                // subdiagonal
                if (r == c + 1) {
                    if (compare(data[c][r], 0)) {
                        return false;
                    }
                }
                else {
                    if (!compare(data[c][r], 0)) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    bool isLowerHessenberg() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = c - 2; r >= 0; r--) {
                if (!compare(data[c][r], 0))
                    return false;
            }
        }

        return true;
    }

    bool isUnreducedLowerHessenberg() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = c - 1; r >= 0; r--) {
                // supra-diagonal
                if (r == c - 1) {
                    if (compare(data[c][r], 0)) {
                        return false;
                    }
                }
                else {
                    if (!compare(data[c][r], 0)) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    bool isTridiagonal() const requires (isSquare) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                // diagonal, subdiagonal, supradiagonal
                if (c == r || c == r + 1 || c == r - 1)
                    continue;

                if (!compare(data[c][r], 0))
                    return false;
            }
        }

        return true;
    }

#pragma region eigen stuffs
    template<typename T_TYPE, typename Q_TYPE>
    struct LanczosAlgorithm {
        T_TYPE t;
        Q_TYPE q;
    };

    template<int ITER>
    LanczosAlgorithm<Matrix<ITER, ITER, T>, Matrix<ITER + 1, COLUMNS, T>> lanczosAlgorithm() const requires (isSquare) {
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

    /**
     * Used to get an eigen-value approximation from an eigen-vector approximation
     * @param vec eigen-vector approximation
     * @return corresponding eigen-value approximation for given vector @a vec
     */
    T rayleighQuotient(const Vector<COLUMNS, T>& vec) const {
        return (vec * *this * vec) / (vec * vec);
    }

    template<typename EIGENVECTOR_TYPE, typename EIGENVALUE_TYPE>
    struct PowerIteration {
        EIGENVECTOR_TYPE vector;
        EIGENVALUE_TYPE value;
    };

    PowerIteration<Vector<COLUMNS, T>, T> powerIteration(const int maxIterations, const T tolerance = epsilon) const {
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
#pragma endregion
};

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> multiply(const T lhs, const Matrix<COLUMNS, ROWS, T>& rhs) {
    return rhs.multiply(lhs);
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> operator*(const T lhs, const Matrix<COLUMNS, ROWS, T>& rhs) {
    return multiply(lhs, rhs);
}

template<int COLUMNS, int ROWS, is_scalar_v T, is_convertable_to<T> OTHER_T>
Matrix<COLUMNS, ROWS, T> multiply(const OTHER_T lhs, const Matrix<COLUMNS, ROWS, T>& rhs) {
    return rhs.multiply(lhs);
}

template<int COLUMNS, int ROWS, is_scalar_v T, is_convertable_to<T> OTHER_T>
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
template<int OTHER_N, is_convertable_to<T> OTHER_T>
Matrix<OTHER_N, N, T> Vector<N, T>::outerProductMatrix(const Vector<OTHER_N, OTHER_T>& v) const {
    Matrix<OTHER_N, N, T> result;

    for (int c = 0; c < result.columns; c++) {
        for (int r = 0; r < result.rows; r++) {
            if (isComplex)
                result[c][r] = data[r] * std::conj(v[c]);
            else
                result[c][r] = data[r] * v[c];
        }
    }

    return result;
}

template<int N, is_scalar_v T>
template<int ROWS>
Vector<N, T> Vector<N, T>::multiply(const Matrix<N, ROWS, T>& m) const {
    Vector<N, T> result;

    for (int c = 0; c < N; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[r] += data[c] * m.data[c][r];
        }
    }

    return result;
}

template<int N, is_scalar_v T>
template<int ROWS>
Vector<N, T> Vector<N, T>::operator*(const Matrix<N, ROWS, T>& m) const {
    return multiply(m);
}

template<int N, is_scalar_v T>
template<int ROWS, is_convertable_to<T> OTHER_T>
Vector<N, T> Vector<N, T>::multiply(const Matrix<N, ROWS, OTHER_T>& m) const {
    Vector<N, T> result;

    for (int c = 0; c < N; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[r] += data[c] * m.data[c][r];
        }
    }

    return result;
}

template<int N, is_scalar_v T>
template<int ROWS, is_convertable_to<T> OTHER_T>
Vector<N, T> Vector<N, T>::operator*(const Matrix<N, ROWS, OTHER_T>& m) const {
    return multiply(m);
}

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
