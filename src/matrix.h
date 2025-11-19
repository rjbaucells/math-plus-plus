#pragma once
#include <cassert>
#include <complex>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include "complex.h"
#include "vector.h"

template<int COLUMNS, int ROWS, typename T = float>
struct Matrix {
    static constexpr int rows = ROWS;
    static constexpr int columns = COLUMNS;
    static constexpr bool square = ROWS == COLUMNS;

    T data[COLUMNS][ROWS] = {};

    Matrix() = default;

    Matrix(std::initializer_list<std::initializer_list<T>> initializerList) {
        if (initializerList.size() != ROWS) {
            throw std::runtime_error("Incorrect number of rows in initializer list");
        }

        int r = 0;
        for (const auto& row : initializerList) {
            if (row.size() != COLUMNS) {
                throw std::runtime_error("Incorrect number of columns in initializer list");
            }

            int c = 0;

            for (const auto element : row) {
                data[c][r] = element;
                c++;
            }

            r++;
        }
    }

    // copy constructor with same type
    Matrix(const Matrix<COLUMNS, ROWS, T>& other) {
        memcpy(data, other.data, sizeof(T) * ROWS * COLUMNS);
    }

    // copy constructor with different type
    template<IsConvertableTo<T> OTHER_T>
    Matrix(const Matrix<COLUMNS, ROWS, OTHER_T>& other) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                data[c][r] = other.data[c][r];
            }
        }
    }

#pragma region Same Type Operators
    // copy assignment operator with same type
    Matrix<COLUMNS, ROWS, T>& operator=(const Matrix<COLUMNS, ROWS, T>& other) {
        if (this != &other) {
            memcpy(data, other.data, sizeof(T) * COLUMNS * ROWS);
        }

        return *this;
    }

    Matrix<COLUMNS, ROWS, T> add(const Matrix<COLUMNS, ROWS, T>& other) const {
        Matrix<COLUMNS, ROWS, T> result;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                result[c][r] = data[c][r] + other.data[c][r];
            }
        }

        return result;
    }

    Matrix<COLUMNS, ROWS, T> operator+(const Matrix<COLUMNS, ROWS, T>& other) const {
        return add(other);
    }

    Matrix<COLUMNS, ROWS, T>& addEquals(const Matrix<COLUMNS, ROWS, T>& other) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                data[c][r] += other.data[c][r];
            }
        }

        return *this;
    }

    Matrix<COLUMNS, ROWS, T>& operator+=(const Matrix<COLUMNS, ROWS, T>& other) {
        return addEquals(other);
    }

    Matrix<COLUMNS, ROWS, T> subtract(const Matrix<COLUMNS, ROWS, T>& other) const {
        Matrix<COLUMNS, ROWS, T> result;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                result[c][r] = data[c][r] - other.data[c][r];
            }
        }

        return result;
    }

    Matrix<COLUMNS, ROWS, T> operator-(const Matrix<COLUMNS, ROWS, T>& other) const {
        return subtract(other);
    }

    Matrix<COLUMNS, ROWS, T>& subtractEquals(const Matrix<COLUMNS, ROWS, T>& other) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                data[c][r] -= other.data[c][r];
            }
        }

        return *this;
    }

    Matrix<COLUMNS, ROWS, T> operator-=(const Matrix<COLUMNS, ROWS, T>& other) {
        return subtractEquals(other);
    }

    template<int OTHER_COLUMNS>
    Matrix<OTHER_COLUMNS, ROWS, T> multiply(const Matrix<OTHER_COLUMNS, COLUMNS, T>& other) const {
        Matrix<OTHER_COLUMNS, ROWS, T> result;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                for (int x = 0; x < COLUMNS; x++) {
                    result[c][r] += data[x][r] * other.data[c][x];
                }
            }
        }

        return result;
    }

    template<int OTHER_COLUMNS>
    Matrix<OTHER_COLUMNS, ROWS, T> operator*(const Matrix<OTHER_COLUMNS, COLUMNS, T>& other) const {
        return multiply(other);
    }

    Matrix<COLUMNS, ROWS, T> multiply(const T& val) const {
        Matrix<COLUMNS, ROWS, T> result;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                result[c][r] = data[c][r] * val;
            }
        }

        return result;
    }

    Matrix<COLUMNS, ROWS, T> operator*(const T& val) const {
        return multiply(val);
    }

    Matrix<COLUMNS, ROWS, T>& multiplyEquals(const T& val) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                data[c][r] *= val;
            }
        }

        return *this;
    }

    Matrix<COLUMNS, ROWS, T>& operator*=(const T& val) const {
        return multiplyEquals(val);
    }

    bool compare(const Matrix<COLUMNS, ROWS, T>& other) const {
        if constexpr (std::is_floating_point_v<T>) {
            T epsilon = std::numeric_limits<T>::epsilon();
            for (int c = 0; c < COLUMNS; c++) {
                for (int r = 0; r < ROWS; r++) {
                    if (std::abs(other.data[c][r] - data[c][r]) > epsilon)
                        return false;
                }
            }

            return true;
        }
        else {
            for (int c = 0; c < COLUMNS; c++) {
                for (int r = 0; r < ROWS; r++) {
                    if (other.data[c][r] == data[c][r])
                        return false;
                }
            }

            return false;
        }
    }

    bool operator==(const Matrix<COLUMNS, ROWS, T>& other) const {
        return compare(other);
    }

#pragma endregion
#pragma region Different Type Operators
    // copy assignment operator with different type
    template<IsConvertableTo<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& operator=(const Matrix<COLUMNS, ROWS, OTHER_T>& other) {
        if (*this != other) {
            for (int c = 0; c < COLUMNS; c++) {
                for (int r = 0; r < ROWS; r++) {
                    data[c][r] = other.data[c][r];
                }
            }
        }

        return *this;
    }

    template<IsConvertableTo<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> add(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const {
        Matrix<COLUMNS, ROWS, T> result;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                result[c][r] = data[c][r] + other.data[c][r];
            }
        }

        return result;
    }

    template<IsConvertableTo<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> operator+(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const {
        return add(other);
    }

    template<IsConvertableTo<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& addEquals(const Matrix<COLUMNS, ROWS, OTHER_T>& other) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                data[c][r] += other.data[c][r];
            }
        }

        return *this;
    }

    template<IsConvertableTo<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& operator+=(const Matrix<COLUMNS, ROWS, OTHER_T>& other) {
        return addEquals(other);
    }

    template<IsConvertableTo<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> subtract(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const {
        Matrix<COLUMNS, ROWS, T> result;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                result[c][r] = data[c][r] - other.data[c][r];
            }
        }

        return result;
    }

    template<IsConvertableTo<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> operator-(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const {
        return subtract(other);
    }

    template<IsConvertableTo<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& subtractEquals(const Matrix<COLUMNS, ROWS, OTHER_T>& other) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                data[c][r] -= other.data[c][r];
            }
        }

        return *this;
    }

    template<IsConvertableTo<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& operator-=(const Matrix<COLUMNS, ROWS, OTHER_T>& other) {
        return subtractEquals(other);
    }

    template<int OTHER_COLUMNS, IsConvertableTo<T> OTHER_T>
    Matrix<OTHER_COLUMNS, ROWS, T> multiply(const Matrix<OTHER_COLUMNS, COLUMNS, OTHER_T>& other) const {
        Matrix<OTHER_COLUMNS, ROWS, T> result;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                for (int x = 0; x < COLUMNS; x++) {
                    result[c][r] += data[x][r] * other.data[c][x];
                }
            }
        }

        return result;
    }

    template<int OTHER_COLUMNS, IsConvertableTo<T> OTHER_T>
    Matrix<OTHER_COLUMNS, ROWS, T> operator*(const Matrix<OTHER_COLUMNS, COLUMNS, OTHER_T>& other) const {
        return multiply(other);
    }

    template<IsConvertableTo<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> multiply(const OTHER_T& val) const {
        Matrix<COLUMNS, ROWS, T> result;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                result[c][r] = data[c][r] * val;
            }
        }

        return result;
    }

    template<IsConvertableTo<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T> operator*(const OTHER_T& val) const {
        return multiply(val);
    }

    template<IsConvertableTo<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& multiplyEquals(const OTHER_T& val) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                data[c][r] *= val;
            }
        }

        return *this;
    }

    template<IsConvertableTo<T> OTHER_T>
    Matrix<COLUMNS, ROWS, T>& operator*=(const OTHER_T& val) const {
        return multiplyEquals(val);
    }

    template<IsConvertableTo<T> OTHER_T>
    bool compare(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const {
        if constexpr (std::is_floating_point_v<T>) {
            T epsilon = std::numeric_limits<T>::epsilon();
            for (int c = 0; c < COLUMNS; c++) {
                for (int r = 0; r < ROWS; r++) {
                    if (std::abs(static_cast<T>(other.data[c][r]) - data[c][r]) > epsilon)
                        return false;
                }
            }

            return true;
        }
        else {
            for (int c = 0; c < COLUMNS; c++) {
                for (int r = 0; r < ROWS; r++) {
                    if (static_cast<T>(other.data[c][r]) == data[c][r])
                        return false;
                }
            }

            return false;
        }
    }

    template<IsConvertableTo<T> OTHER_T>
    bool operator==(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const {
        return compare(other);
    }

#pragma endregion

    // subscript
    T* operator[](const int index) {
        return &data[index][0];
    }

    // const subscript
    const T* operator[](const int index) const {
        return &data[index][0];
    }

    Matrix<ROWS, COLUMNS, T> transpose() const requires (!IsComplex<T>::value) {
        Matrix<ROWS, COLUMNS, T> result;

        for (int c = 0; c < ROWS; c++) {
            for (int r = 0; r < COLUMNS; r++) {
                result[c][r] = data[r][c];
            }
        }

        return result;
    }

    Matrix<ROWS, COLUMNS, T> conjugateTranspose() const requires (IsComplex<T>::value) {
        Matrix<ROWS, COLUMNS, T> result;

        for (int c = 0; c < ROWS; c++) {
            for (int r = 0; r < COLUMNS; r++) {
                result[c][r].real = data[r][c].real;
                result[c][r].real = data[r][c].imag;
            }
        }

        return result;
    }

    Matrix<COLUMNS, ROWS, T> inverse() const requires (square) {
        // its a one by one, we can just return 1 / value
        if constexpr (ROWS == 1) {
            if (data[0][0] == 0) {
                throw std::runtime_error("Cannot find inverse of singular matrix");
            }

            Matrix<1, 1, T> result;

            result[0][0] = 1 / data[0][0];
            return result;
        }

        // its a two by two, we can do the special fast thing
        if constexpr (ROWS == 2) {
            T det = determinant();

            if (det == 0) {
                throw std::runtime_error("Cannot find inverse of singular matrix");
            }

            Matrix<2, 2, T> result;

            result[0][0] = data[1][1] / det;
            result[0][1] = -data[0][1] / det;
            result[1][0] = -data[1][0] / det;
            result[1][1] = data[0][0] / det;

            return result;
        }

        const Matrix<COLUMNS, ROWS, T> identityMatrix = Matrix<COLUMNS, ROWS, T>::identity();
        Matrix<COLUMNS * 2, ROWS, T> temp;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                temp[c][r] = data[c][r];
            }
        }

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                temp[c + COLUMNS][r] = identityMatrix[c][r];
            }
        }

        for (int c = 0; c < COLUMNS; c++) {
            // if the pivot is zero
            if (temp[c][c] == 0) {
                int rowIndex = -1;
                for (int r = c + 1; r < ROWS; r++) {
                    if (temp[c][r] != 0 && (rowIndex == -1 || std::abs(temp[c][r]) > std::abs(temp[c][rowIndex]))) {
                        rowIndex = r;
                    }
                }

                if (rowIndex == -1) {
                    throw std::runtime_error("Cannot find inverse of singular matrix");
                }

                // since c (column) = r (row), the dest row is at temp[row (aka, c}] the src row is the memory address of the start of the biggest row so &temp[c][biggestRow]
                temp = temp.swapRows(c, rowIndex);
            }

            // normalize pivot row to 1
            T value = temp[c][c];
            for (int cc = c + 1; cc < temp.columns; cc++) {
                temp[cc][c] /= value;
            }
            temp[c][c] = 1;

            for (int r = 0; r < ROWS; r++) {
                // the current value is on the main diagonal. (expected to be a 1). Should be good
                if (r == c) {
                    // we dont want to affect the pivot
                    continue;
                }
                // it isnt on the main diagonal, (expected to be a 0)
                T k = temp[c][r];
                for (int cc = c + 1; cc < temp.columns; cc++) {
                    temp[cc][r] += -k * temp[cc][c];
                }
                temp[c][r] = 0;
            }
        }

        Matrix<COLUMNS, ROWS, T> result;

        for (int c = 0; c < COLUMNS; ++c) {
            for (int r = 0; r < ROWS; r++) {
                result[c][r] = temp[c + COLUMNS][r];
            }
        }

        return result;
    }

    T determinant() const requires (square) {
        if constexpr (ROWS == 1) {
            return data[0][0];
        }
        else if constexpr (ROWS == 2) {
            return data[0][0] * data[1][1] - data[0][1] * data[1][0];
        }
        else {
            T result = 0;
            int sign = 1;

            for (int c = 0; c < COLUMNS; c++) {
                Matrix<COLUMNS - 1, ROWS - 1, T> insideMatrix = getSubMatrix(0, c);

                result += sign * data[c][0] * insideMatrix.determinant();
                sign *= -1;
            }

            return result;
        }
    }

#pragma region 4x4 stuffs
    static Matrix<4, 4, T> ortho(const T left, const T right, const T bottom, const T top, const T near, const T far) {
        // identity
        Matrix<4, 4, T> transformation = Matrix<4, 4, T>::identity();
        // transformation
        transformation.data[0][0] = 2.0f / (right - left);
        transformation.data[1][1] = 2.0f / (top - bottom);
        transformation.data[2][2] = -2.0f / (far - near);
        transformation.data[3][0] = -(right + left) / (right - left);
        transformation.data[3][1] = -(top + bottom) / (top - bottom);
        transformation.data[3][2] = -(far + near) / (far - near);
        // return
        return transformation;
    }

#pragma endregion
    Vector<COLUMNS, T> multiply(const Vector<COLUMNS, T>& other) {
        Vector<COLUMNS, T> result;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                result[r] += data[c][r] * other[c];
            }
        }

        return result;
    }

    template<IsConvertableTo<T> OTHER_T>
    Vector<COLUMNS, T> operator*(const Vector<COLUMNS, OTHER_T>& other) {
        return multiply(other);
    }

    template<IsConvertableTo<T> OTHER_T>
    Vector<COLUMNS, T> multiply(const Vector<COLUMNS, OTHER_T>& other) {
        Vector<COLUMNS, T> result;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                result[r] += data[c][r] * other[c];
            }
        }

        return result;
    }

    Vector<COLUMNS, T> operator*(const Vector<COLUMNS, T>& other) {
        return multiply(other);
    }

    std::string toString() const {
        std::stringstream ss;
        ss.precision(2);

        ss << "[";

        for (int r = 0; r < ROWS; r++) {
            ss << "[";
            for (int c = 0; c < COLUMNS; c++) {
                ss << data[c][r];

                if (c < COLUMNS - 1)
                    ss << ", ";
            }

            ss << "]";

            if (r < ROWS - 1)
                ss << ", ";
        }

        ss << "]";
        return ss.str();
    }

    std::string toLaTex() const {
        std::stringstream ss;
        ss.precision(2);

        ss << "\\begin{bmatrix}";

        for (int r = 0; r < ROWS; r++) {
            for (int c = 0; c < COLUMNS; c++) {
                ss << data[c][r];

                if (c < COLUMNS - 1)
                    ss << " & ";
            }

            if (r < ROWS - 1)
                ss << "\\\\";
        }

        ss << "\\end{bmatrix}";
        return ss.str();
    }

    Matrix<COLUMNS - 1, ROWS - 1, T> getSubMatrix(const int rowToRemove, const int columnToRemove) const {
        Matrix<COLUMNS - 1, ROWS - 1, T> subMatrix;
        int subMatrixR = 0;
        for (int r = 0; r < ROWS; r++) {
            if (r == rowToRemove)
                continue;

            int subMatrixC = 0;
            for (int c = 0; c < COLUMNS; c++) {
                if (c == columnToRemove)
                    continue;

                subMatrix.data[subMatrixC][subMatrixR] = data[c][r];
                subMatrixC++;
            }
            subMatrixR++;
        }

        return subMatrix;
    }

    template<int NUM_COLUMNS_TO_REMOVE, int NUM_ROWS_TO_REMOVE>
    Matrix<COLUMNS - NUM_COLUMNS_TO_REMOVE, ROWS - NUM_ROWS_TO_REMOVE, T> getSubMatrix(const std::array<int, NUM_COLUMNS_TO_REMOVE>& columnsToRemove, const std::array<int, NUM_ROWS_TO_REMOVE>& rowsToRemove) const {
        Matrix<COLUMNS - NUM_COLUMNS_TO_REMOVE, ROWS - NUM_ROWS_TO_REMOVE, T> subMatrix;

        int subMatrixC = 0;
        for (int c = 0; c < COLUMNS; c++) {
            if (std::find(columnsToRemove.begin(), columnsToRemove.end(), c))
                continue;

            int subMatrixR = 0;
            for (int r = 0; r < ROWS; r++) {
                if (std::find(rowsToRemove.begin(), rowsToRemove.end(), r))
                    continue;

                subMatrix[subMatrixC][subMatrixR] = data[c][r];
                subMatrixR++;
            }

            subMatrixC++;
        }

        return subMatrix;
    }

    Matrix<COLUMNS, ROWS, T> swapRows(const int rowA, const int rowB) {
        Matrix<COLUMNS, ROWS, T> m = *this;

        T temp[COLUMNS] = {};

        for (int c = 0; c < COLUMNS; c++) {
            temp[c] = m[c][rowA];
        }

        for (int c = 0; c < COLUMNS; c++) {
            m[c][rowA] = m[c][rowB];
        }

        for (int c = 0; c < COLUMNS; c++) {
            m[c][rowB] = temp[c];
        }

        return m;
    }

    Matrix<COLUMNS, ROWS, T> swapColumns(const int columnA, const int columnB) {
        Matrix<COLUMNS, ROWS, T> m = *this;

        T temp[ROWS] = {};

        for (int r = 0; r < ROWS; r++) {
            temp[r] = m[columnA][r];
        }

        for (int r = 0; r < COLUMNS; r++) {
            m[columnA][r] = m[columnB][r];
        }

        for (int r = 0; r < COLUMNS; r++) {
            m[columnB][r] = temp[r];
        }

        return m;
    }

    explicit operator const T*() const {
        return &data[0][0];
    }

    explicit operator T*() {
        return &data[0][0];
    }

    static Matrix<COLUMNS, ROWS, T> identity() requires (square) {
        Matrix<COLUMNS, ROWS, T> result;

        for (int i = 0; i < COLUMNS; i++) {
            result[i][i] = 1;
        }

        return result;
    }

    bool isRowEchelon(bool pivotMustBeOne = false) const {
        int lastPivotColumn = -1;

        for (int r = 0; r < ROWS; r++) {
            bool foundNonZero = false;

            for (int c = 0; c < COLUMNS; c++) {
                // this is a pivot
                if (!foundNonZero && data[c][r] != 0) {
                    // this is to the left of the last pivot
                    if (c < lastPivotColumn)
                        return false;

                    // the pivot must be one and it isn't
                    if (data[c][r] != 1 && pivotMustBeOne)
                        return false;

                    lastPivotColumn = c;
                }

                if (data[c][r] != 0)
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
                if (rowIndex == -1 && pivot == 0) {
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
                if (!foundNonZero && data[c][r] != 0) {
                    // this is to the left of the last pivot
                    if (c < lastPivotColumn)
                        return false;

                    // the pivot isn't 1
                    if (data[c][r] != 1)
                        return false;

                    // check that no other number in that column is a nonzero value
                    for (int i = 0; i < ROWS; i++) {
                        if (data[c][i] != 0 && data[c][i] != data[c][r])
                            return false;
                    }

                    lastPivotColumn = c;
                }

                if (data[c][r] != 0)
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
                if (rowIndex == -1 && pivot == 0) {
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
            if (ref[0][r] != 0) {
                result++;
            }
            else {
                bool nonZero = false;

                for (int c = 0; c < COLUMNS; c++) {
                    if (ref[c][r] != 0) {
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

    bool symmetrical() const requires (!IsComplex<T>::value && square) {
        return *this == transpose();
    }

    bool skewSymmetrical() const requires (!IsComplex<T>::value && square) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                if (data[r][c] != -data[c][r])
                    return false;
            }
        }

        return true;
    }

    bool hermitian() const requires (IsComplex<T>::value && square) {
        auto ts = conjugateTranspose();

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                if (data[c][r] != ts[c][r]) {
                    return false;
                }
            }
        }

        return true;
    }

    bool skewHermitian() const requires (IsComplex<T>::value && square) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                T complexConjugate = data[c][r];
                complexConjugate.imag *= -1;

                if (data[r][c] != -complexConjugate)
                    return false;
            }
        }

        return true;
    }

    bool positiveDefinite() const requires (square) {
        for (int c = 0; c < COLUMNS; c++) {
            Vector<ROWS> x;

            for (int r = 0; r < ROWS; r++) {
                x[r] = data[c][r];
            }

            if (!(x.componentDot(x) > 0)) {
                return false;
            }
        }

        return true;
    }

    bool positiveSemiDefinite() const requires (square) {
        for (int c = 0; c < COLUMNS; c++) {
            Vector<ROWS> x;

            for (int r = 0; r < ROWS; r++) {
                x[r] = data[c][r];
            }

            if (!(x.componentDot(x) >= 0)) {
                return false;
            }
        }

        return true;
    }

    bool negativeDefinite() const requires (square) {
        for (int c = 0; c < COLUMNS; c++) {
            Vector<ROWS> x;

            for (int r = 0; r < ROWS; r++) {
                x[r] = data[c][r];
            }

            if (!(x.componentDot(x) < 0)) {
                return false;
            }
        }

        return true;
    }

    bool negativeSemiDefinite() const requires (square) {
        for (int c = 0; c < COLUMNS; c++) {
            Vector<ROWS> x;

            for (int r = 0; r < ROWS; r++) {
                x[r] = data[c][r];
            }

            if (!(x.componentDot(x) <= 0)) {
                return false;
            }
        }

        return true;
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

    std::array<Vector<COLUMNS>, ROWS> getRowVectors() const {
        std::array<Vector<COLUMNS>, ROWS> vecs;

        for (int r = 0; r < ROWS; r++) {
            for (int c = 0; c < COLUMNS; c++) {
                vecs[r][c] = data[c][r];
            }
        }

        return vecs;
    }

    T trace() const requires (square) {
        T sum = {};

        for (int c = 0; c < COLUMNS; c++) {
            sum += data[c][c];
        }

        return sum;
    }

    bool unitary() const requires (IsComplex<T>::value && square) {
        auto iden = identity();
        auto tr = conjugateTranspose();

        if (multiply(tr) != iden)
            return false;

        if (tr.multiply(*this) != iden)
            return false;

        return true;
    }

    bool specialUnitary() const requires (IsComplex<T>::value && square) {
        auto iden = identity();
        auto tr = conjugateTranspose();

        if (multiply(tr) != iden)
            return false;

        if (tr.multiply(*this) != iden)
            return false;

        return determinant() == 1;
    }

    bool orthogonal() const requires (!IsComplex<T>::value && square) {
        auto iden = identity();
        auto tr = transpose();

        if (multiply(tr) != iden)
            return false;

        if (tr.multiply(*this) != iden)
            return false;

        return true;
    }

    bool upperTriangleMatrix() const requires (square) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < c; r++) {
                if (data[c][r] != 0)
                    return false;
            }
        }

        return true;
    }

    bool lowerTriangleMatrix() const requires (square) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = c; r < ROWS; r++) {
                if (data[c][r] != 0)
                    return false;
            }
        }

        return true;
    }

    bool diagonalMatrix() const requires (square) {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                if (c == r)
                    continue;

                if (data[c][r] != 0)
                    return false;
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
                if (rowIndex == -1 && pivot == 0) {
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

                if (rowIndex == -1 && columnIndex == -1 && pivot == 0) {
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

            if (pivot == 0)
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

            if (pivot == 0)
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

    CholeskyDecomposition<Matrix<COLUMNS, ROWS, T>, Matrix<ROWS, COLUMNS, T>> choleskyDecomposition(bool allowPositiveSemiDefinite = false) const requires (square) {
        if ((IsComplex<T>::value && !hermitian()) || (!IsComplex<T>::value && !symmetrical())) {
            throw std::runtime_error("Cannot find Cholesky Decomposition of non hermitian/symmetric matrix");
        }

        Matrix<COLUMNS, ROWS, T> l;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r <= c; r++) {
                if constexpr (!IsComplex<T>::value) {
                    if (r == c) {
                        T value = data[c][c];

                        for (int k = 0; k < c; k++) {
                            value -= std::pow(l[k][c], 2);
                        }

                        if (value < 0) {
                            throw std::runtime_error("Cannot cholesky decompose non positive definite matrix");
                        }

                        if (value == 0 && !allowPositiveSemiDefinite) {
                            throw std::runtime_error("Cannot cholesky decompose non semi-positive definite matrix");
                        }

                        l[c][c] = std::sqrt(value);
                    }
                    else if (r > c) {
                        T value = data[c][r];

                        for (int k = 0; k < c; k++) {
                            value -= l[k][r] * l[k][c];
                        }

                        l[c][r] = (1 / l[c][c]) * value;
                    }
                }
                else {
                    if (r == c) {
                        T value = data[c][c];

                        for (int k = 0; k < c; k++) {
                            T complexConjugate = l[k][c];
                            complexConjugate.imag *= -1;
                            value -= l[k][c] * complexConjugate;
                        }

                        if (value < 0) {
                            throw std::runtime_error("Cannot cholesky decompose non positive definite matrix");
                        }

                        if (value == 0 && !allowPositiveSemiDefinite) {
                            throw std::runtime_error("Cannot cholesky decompose non semi-positive definite matrix");
                        }

                        l[c][c] = std::sqrt(value);
                    }
                    else if (r > c) {
                        T value = data[c][r];

                        for (int k = 0; k < c; k++) {
                            T complexConjugate = l[k][c];
                            complexConjugate.imag *= -1;
                            value -= l[k][r] * complexConjugate;
                        }

                        l[c][r] = (1 / l[c][c]) * value;
                    }
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

    LDLDecomposition<Matrix<COLUMNS, ROWS, T>, Matrix<COLUMNS, ROWS, T>, Matrix<ROWS, COLUMNS, T>> ldlDecomposition() const requires (square) {
        Matrix<COLUMNS, ROWS, T> l;
        Matrix<COLUMNS, ROWS, T> d;

        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r <= c; r++) {
                if constexpr (!IsComplex<T>::value) {
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

    QRDecomposition<Matrix<ROWS, ROWS, T>, Matrix<COLUMNS, ROWS, T>> qrDecomposition() const requires (square) {
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

    T rayleighQuotient(const Vector<COLUMNS, T>& x) {
        if constexpr (IsComplex<T>::value) {
            return (x.conjugate() * *this * x) / (x.conjugate() * x);
        }
        else {
            return (x * *this * x) / x * x;
        }
    }

    Matrix<COLUMNS, ROWS, T> cofactor() const requires (square) {
        if constexpr (COLUMNS == 1) {
            return 1;
        }
        else if constexpr (COLUMNS == 2) {
            return {{data[1][1], -data[1][0]}, {-data[0][1], data[0][0]}};
        }
        else if constexpr (COLUMNS >= 3) {
            Matrix<COLUMNS, ROWS, T> result;

            int sign = 1;
            for (int c = 0; c < COLUMNS; c++) {
                for (int r = 0; r < ROWS; r++) {
                    result[c][r] = getSubMatrix(c, r).determinant() * sign;
                    sign *= -1;
                }
            }

            return result;
        }
    }

    Matrix<COLUMNS, ROWS, T> adjugate() const requires (square) {
        return cofactor().transpose();
    }

    T getEigenValueGivenEigenVector(const Vector<COLUMNS, T>& v) {

    }

    template<typename EIGENVALUE_TYPE, typename EIGENVECTOR_TYPE>
    struct BiggestEigenValueAndVector {
        EIGENVALUE_TYPE eigenValue;
        EIGENVECTOR_TYPE eigenVector;
    };

    BiggestEigenValueAndVector<T, Vector<COLUMNS, T>> biggestEigenValueAndVector(const T tolerance = 0.01, const int maxIterations = 100) const {
        // approximation or just a random vector
        Vector<COLUMNS, T> b;
        for (int c = 0; c < COLUMNS; c++) {
            b[c] = c;
        }

        for (int i = 0; i < maxIterations; i++) {
            Vector<COLUMNS> b_i = multiply(b);

            // if we have reached the tolerance
            for (int c = 0; c < COLUMNS; c++) {
                if (b_i[c] - b[c] < tolerance)
                    return b_i;
            }

            // re-normalize
            Vector<COLUMNS> normB_i = b_i.normalize();
            b = b_i / normB_i;
        }

        return {getEigenValueGivenEigenVector(b), b};
    }
};

template<int COLUMNS, int ROWS, typename T>
Vector<COLUMNS, T> multiply(const Vector<COLUMNS, T>& v, const Matrix<COLUMNS, ROWS, T>& m) {
    Vector<COLUMNS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[r] += v[c] * m.data[c][r];
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, typename T, typename OTHER_T>
Vector<COLUMNS, T> multiply(const Vector<COLUMNS, OTHER_T>& v, const Matrix<COLUMNS, ROWS, T>& m) {
    Vector<COLUMNS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[r] += v[c] * m.data[c][r];
        }
    }

    return result;
}