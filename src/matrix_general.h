#pragma once
#include "matrix.h"
#include <cstring>

template<int COLUMNS, int ROWS, is_scalar_v T>
constexpr Matrix<COLUMNS, ROWS, T>::Matrix(std::initializer_list<std::initializer_list<T>> initializerList) {
    int r = 0;
    for (const auto& row : initializerList) {
        int c = 0;

        for (const auto element : row) {
            data[c][r] = element;
            c++;
        }

        r++;
    }
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T>::Matrix(const Matrix<COLUMNS, ROWS, T>& other) {
    memcpy(data, other.data, sizeof(T) * ROWS * COLUMNS);
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T> requires std::convertible_to<OTHER_T, T>
Matrix<COLUMNS, ROWS, T>::Matrix(const Matrix<COLUMNS, ROWS, OTHER_T>& other) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            data[c][r] = other.data[c][r];
        }
    }
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<ROWS, COLUMNS, T> Matrix<COLUMNS, ROWS, T>::transpose() const {
    Matrix<ROWS, COLUMNS, T> result;

    for (int c = 0; c < ROWS; c++) {
        for (int r = 0; r < COLUMNS; r++) {
            result[c][r] = data[r][c];
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<ROWS, COLUMNS, T> Matrix<COLUMNS, ROWS, T>::conjugateTranspose() const {
    if constexpr (!isComplex) {
        return transpose();
    }
    else {
        Matrix<ROWS, COLUMNS, T> result;

        for (int c = 0; c < ROWS; c++) {
            for (int r = 0; r < COLUMNS; r++) {
                result[c][r] = std::conj(data[r][c]);
            }
        }

        return result;
    }
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::inverse() const requires (isSquare) {
    if constexpr (ROWS == 1) { // its a one by one, we can just return 1 / value
        if (compare(data[0][0], 0)) {
            throw std::runtime_error("Cannot find inverse of singular matrix");
        }

        Matrix<1, 1, T> result;

        result[0][0] = 1 / data[0][0];
        return result;
    }
    else if constexpr (ROWS == 2) { // its a two by two, we can do the special fast thing
        T det = determinant();

        if (compare(det, 0)) {
            throw std::runtime_error("Cannot find inverse of singular matrix");
        }

        Matrix<2, 2, T> result;

        result[0][0] = data[1][1] / det;
        result[0][1] = -data[0][1] / det;
        result[1][0] = -data[1][0] / det;
        result[1][1] = data[0][0] / det;

        return result;
    }
    else {
        Matrix<COLUMNS, ROWS, T> left = *this;
        Matrix<COLUMNS, ROWS, T> right = identity();

        for (int c = 0; c < COLUMNS; c++) {
            // handle row swaps for performance
            int rowIndex = -1;
            T value = left[c][c];

            // only check for rows to swap with below main diagonal
            for (int i = c; i < ROWS; i++) {
                T norm = std::norm(left[c][i]);

                if (norm > value) {
                    value = norm;
                    rowIndex = i;
                }
            }

            // if we didnt swap, but needed to have
            if (rowIndex == -1 && compare(value, 0))
                throw std::runtime_error("Cannot find inverse of singular matrix");

            // swap if we can
            if (rowIndex != -1) {
                left = left.swapRows(c, rowIndex);
                right = right.swapRows(c, rowIndex);
            }

            {
                // normalize pivot row
                T pivot = left[c][c];
                for (int cc = c; cc < COLUMNS; cc++) {
                    left[cc][c] /= pivot;
                    right[cc][c] /= pivot;
                }
            }

            for (int r = 0; r < ROWS; r++) {
                if (c == r)
                    continue;

                T multiplierToPivotRow = left[c][r];
                for (int i = c; i < COLUMNS; i++) {
                    left[i][r] += -multiplierToPivotRow * left[i][c];
                    right[i][r] += -multiplierToPivotRow * left[i][c];
                }
            }
        }

        return right;
    }
}

enum DeterminantAlgorithm {
    laplace,
    triangular,
    hessenberg
};

template<int COLUMNS, int ROWS, is_scalar_v T>
T Matrix<COLUMNS, ROWS, T>::determinant(const DeterminantAlgorithm algorithm) const requires (isSquare) {
    if constexpr (COLUMNS == 1) {
        return data[0][0];
    }
    else if constexpr (COLUMNS == 2) {
        // ad - bc
        return data[0][0] * data[1][1] - data[0][1] * data[1][0];
    }
    else if constexpr (COLUMNS == 3) {
        // aei + bfg + cdh - ceg - bdi - afh
        return (data[0][0] * data[1][1] * data[2][2]) + (data[1][0] * data[2][1] * data[0][2]) + (data[2][0] * data[0][1] * data[1][2]) - (data[2][0] * data[1][1] * data[0][2]) - (data[1][0] * data[0][1] * data[2][2]) - (data[0][0] * data[2][1] * data[1][2]);
    }
    else {
        switch (algorithm) {
            case triangular:
                return triangularDeterminant();
            case tridiagonal:
                return tridiagonalDeterminant();
            default:
                return laplaceDeterminant();
        }
    }
}

template<int COLUMNS, int ROWS, is_scalar_v T>
T Matrix<COLUMNS, ROWS, T>::laplaceDeterminant() const requires (isSquare) {
    T result = {};
    int sign = 1;

    for (int c = 0; c < COLUMNS; c++) {
        Matrix<COLUMNS - 1, ROWS - 1, T> insideMatrix = removeColumnAndRow(c, 0);

        result += sign * data[c][0] * insideMatrix.determinant();
        sign *= -1;
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
T Matrix<COLUMNS, ROWS, T>::triangularDeterminant() const requires (isSquare) {
    T result = {};

    for (int i = 0; i < COLUMNS; i++) {
        result *= data[i][i];
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
T Matrix<COLUMNS, ROWS, T>::tridiagonalDeterminant() const requires (isSquare) {

}


template<int COLUMNS, int ROWS, is_scalar_v T>
std::string Matrix<COLUMNS, ROWS, T>::toString() const {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
std::string Matrix<COLUMNS, ROWS, T>::toLaTex() const {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
template<int NUM_COLUMNS_TO_REMOVE>
Matrix<COLUMNS - NUM_COLUMNS_TO_REMOVE, ROWS, T> Matrix<COLUMNS, ROWS, T>::removeColumns(const std::array<int, NUM_COLUMNS_TO_REMOVE>& columnsToRemove) const {
    Matrix<COLUMNS - NUM_COLUMNS_TO_REMOVE, ROWS, T> m;

    for (int c = 0; c < COLUMNS; c++) {
        if (std::find(columnsToRemove.begin(), columnsToRemove.end(), c))
            continue;

        for (int r = 0; r < ROWS; r++) {
            m[c][r] = data[c][r];
        }
    }

    return m;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<int NUM_ROWS_TO_REMOVE>
Matrix<COLUMNS, ROWS - NUM_ROWS_TO_REMOVE, T> Matrix<COLUMNS, ROWS, T>::removeRows(const std::array<int, NUM_ROWS_TO_REMOVE>& rowsToRemove) const {
    Matrix<COLUMNS, ROWS - NUM_ROWS_TO_REMOVE, T> m;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            if (std::find(rowsToRemove.begin(), rowsToRemove.end(), r))
                continue;

            m[c][r] = data[c][r];
        }
    }

    return m;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<int NUM_COLUMNS_TO_REMOVE, int NUM_ROWS_TO_REMOVE>
Matrix<COLUMNS - NUM_COLUMNS_TO_REMOVE, ROWS - NUM_ROWS_TO_REMOVE, T> Matrix<COLUMNS, ROWS, T>::removeColumnsAndRows(const std::array<int, NUM_COLUMNS_TO_REMOVE>& columnsToRemove, const std::array<int, NUM_ROWS_TO_REMOVE>& rowsToRemove) const {
    Matrix<COLUMNS - NUM_COLUMNS_TO_REMOVE, ROWS - NUM_ROWS_TO_REMOVE, T> m;

    for (int c = 0; c < COLUMNS; c++) {
        if (std::find(columnsToRemove.begin(), columnsToRemove.end(), c))
            continue;

        for (int r = 0; r < ROWS; r++) {
            if (std::find(rowsToRemove.begin(), rowsToRemove.end(), r))
                continue;

            m[c][r] = data[c][r];
        }
    }

    return m;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS - 1, ROWS, T> Matrix<COLUMNS, ROWS, T>::removeColumn(const int columnToRemove) const {
    Matrix<COLUMNS - 1, ROWS, T> m;

    for (int c = 0; c < COLUMNS; c++) {
        if (c == columnToRemove)
            continue;

        for (int r = 0; r < ROWS; r++) {
            m[c][r] = data[c][r];
        }
    }

    return m;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS - 1, T> Matrix<COLUMNS, ROWS, T>::removeRow(const int rowToRemove) const {
    Matrix<COLUMNS, ROWS - 1, T> m;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            if (r == rowToRemove)
                continue;

            m[c][r] = data[c][r];
        }
    }

    return m;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS - 1, ROWS - 1, T> Matrix<COLUMNS, ROWS, T>::removeColumnAndRow(const int columnToRemove, const int rowToRemove) const {
    Matrix<COLUMNS - 1, ROWS - 1, T> m;

    for (int c = 0; c < COLUMNS; c++) {
        if (c == columnToRemove)
            continue;

        for (int r = 0; r < ROWS; r++) {
            if (r == rowToRemove)
                continue;

            m[c][r] = data[c][r];
        }
    }

    return m;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::swapRows(const int rowA, const int rowB) {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::swapColumns(const int columnA, const int columnB) {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
constexpr Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::identity() requires (isSquare) {
    Matrix<COLUMNS, ROWS, T> result;

    for (int i = 0; i < COLUMNS; i++) {
        result[i][i] = 1;
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::toRowEchelon() const {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::toReducedRowEchelon() const {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
int Matrix<COLUMNS, ROWS, T>::rank() const {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
Vector<ROWS, T> Matrix<COLUMNS, ROWS, T>::getColumnVector(const int i) const {
    Vector<ROWS, T> v;

    for (int j = 0; j < ROWS; j++) {
        v[j] = data[i][j];
    }

    return v;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
std::array<Vector<ROWS>, COLUMNS> Matrix<COLUMNS, ROWS, T>::getColumnVectors() const {
    std::array<Vector<ROWS>, COLUMNS> vecs;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            vecs[c][r] = data[c][r];
        }
    }

    return vecs;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Vector<COLUMNS, T> Matrix<COLUMNS, ROWS, T>::getRowVector(const int i) const {
    Vector<COLUMNS, T> v;

    for (int j = 0; j < COLUMNS; j++) {
        v[j] = data[j][i];
    }

    return v;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
std::array<Vector<COLUMNS>, ROWS> Matrix<COLUMNS, ROWS, T>::getRowVectors() const {
    std::array<Vector<COLUMNS>, ROWS> vecs;

    for (int r = 0; r < ROWS; r++) {
        for (int c = 0; c < COLUMNS; c++) {
            vecs[r][c] = data[c][r];
        }
    }

    return vecs;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
void Matrix<COLUMNS, ROWS, T>::setColumnVectors(const std::array<Vector<ROWS, T>, COLUMNS>& columnVectors) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            data[c][r] = columnVectors[c][r];
        }
    }
}

template<int COLUMNS, int ROWS, is_scalar_v T>
void Matrix<COLUMNS, ROWS, T>::setColumnVector(const int i, const Vector<ROWS, T>& v) {
    for (int j = 0; j < ROWS; j++) {
        data[i][j] = v[j];
    }
}

template<int COLUMNS, int ROWS, is_scalar_v T>
void Matrix<COLUMNS, ROWS, T>::setRowVectors(const std::array<Vector<COLUMNS, T>, ROWS> rowVectors) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            data[c][r] = rowVectors[r][c];
        }
    }
}

template<int COLUMNS, int ROWS, is_scalar_v T>
void Matrix<COLUMNS, ROWS, T>::setRowVector(const int i, const Vector<COLUMNS, T>& v) {
    for (int j = 0; j < COLUMNS; j++) {
        data[j][i] = v[j];
    }
}

template<int COLUMNS, int ROWS, is_scalar_v T>
T Matrix<COLUMNS, ROWS, T>::trace() const requires (isSquare) {
    T sum = {};

    for (int c = 0; c < COLUMNS; c++) {
        sum += data[c][c];
    }

    return sum;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
T Matrix<COLUMNS, ROWS, T>::minorOfElement(const int c, const int r) const requires (isSquare) {
    return removeColumnAndRow(c, r).determinant();
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::minorMatrix() const requires (isSquare) {
    Matrix<COLUMNS, ROWS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = minorOfElement(c, r);
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
T Matrix<COLUMNS, ROWS, T>::cofactorOfElement(const int c, const int r) const requires (isSquare) {
    return minorOfElement(c, r) * std::pow(-1, c + r);
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::cofactorMatrix() const requires (isSquare) {
    Matrix<COLUMNS, ROWS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = cofactorOfElement(c, r);
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<ROWS, COLUMNS, T> Matrix<COLUMNS, ROWS, T>::adjoint() const requires (isSquare) {
    return cofactorMatrix().transpose();
}

template<int COLUMNS, int ROWS, is_scalar_v T>
Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::hadamardProduct(const Matrix<COLUMNS, ROWS, T>& other) const {
    Matrix<COLUMNS, ROWS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = data[c][r] * other[c][r];
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<typename OTHER_T> requires has_common_type<OTHER_T, T>
Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::hadamardProduct(const Matrix<COLUMNS, ROWS, OTHER_T>& other) const {
    Matrix<COLUMNS, ROWS, std::common_type_t<T, OTHER_T>> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            result[c][r] = data[c][r] * other[c][r];
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<int OTHER_COLUMNS, int OTHER_ROWS>
Matrix<COLUMNS * OTHER_COLUMNS, ROWS * OTHER_ROWS, T> Matrix<COLUMNS, ROWS, T>::kroneckerProduct(const Matrix<OTHER_COLUMNS, OTHER_ROWS, T>& other) const {
    Matrix<COLUMNS * OTHER_COLUMNS, ROWS * OTHER_ROWS, T> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            const T val = data[c][r];

            for (int cc = 0; cc < OTHER_COLUMNS; cc++) {
                for (int rr = 0; rr < OTHER_ROWS; rr++) {
                    result[c * OTHER_COLUMNS + cc][r * OTHER_ROWS + rr] = val * other[cc][rr];
                }
            }
        }
    }

    return result;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
template<int OTHER_COLUMNS, int OTHER_ROWS, typename OTHER_T> requires has_common_type<OTHER_T, T>
Matrix<COLUMNS * OTHER_COLUMNS, ROWS * OTHER_ROWS, std::common_type_t<T, OTHER_T>> Matrix<COLUMNS, ROWS, T>::kroneckerProduct(const Matrix<OTHER_COLUMNS, OTHER_ROWS, OTHER_T>& other) const {
    Matrix<COLUMNS * OTHER_COLUMNS, ROWS * OTHER_ROWS, std::common_type_t<T, OTHER_T>> result;

    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            const T val = data[c][r];

            for (int cc = 0; cc < OTHER_COLUMNS; cc++) {
                for (int rr = 0; rr < OTHER_ROWS; rr++) {
                    result[c * OTHER_COLUMNS + cc][r * OTHER_ROWS + rr] = val * other[cc][rr];
                }
            }
        }
    }

    return result;
}
