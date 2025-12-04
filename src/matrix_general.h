#pragma once

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
template<is_convertable_to<T> OTHER_T>
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

    Matrix<ROWS, COLUMNS, T> result;

    for (int c = 0; c < ROWS; c++) {
        for (int r = 0; r < COLUMNS; r++) {
            result[c][r] = std::conj(data[r][c]);
        }
    }

    return result;
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

    switch (algorithm) {
        case triangular:
            return triangularDeterminant();
        case hessenberg:
            return 0;
        default:
            return laplaceDeterminant();
    }
}

template<int COLUMNS, int ROWS, is_scalar_v T>
T Matrix<COLUMNS, ROWS, T>::laplaceDeterminant() const requires (isSquare) {
    T result = {};
    int sign = 1;

    for (int c = 0; c < COLUMNS; c++) {
        Matrix<COLUMNS - 1, ROWS - 1, T> insideMatrix = removeColumnsAndRows({c}, {0});

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
    Matrix<COLUMNS, ROWS - NUM_ROWS_TO_REMOVE, T> m;

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
static Matrix<COLUMNS, ROWS, T> Matrix<COLUMNS, ROWS, T>::identity() requires (isSquare) {
    Matrix<COLUMNS, ROWS, T> result;

    for (int i = 0; i < COLUMNS; i++) {
        result[i][i] = 1;
    }

    return result;
}
