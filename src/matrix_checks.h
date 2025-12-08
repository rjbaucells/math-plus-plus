#pragma once
#include <complex>

#include "matrix.h"

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isRowEchelon(bool pivotMustBeOne) const {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isReducedRowEchelon() const {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isSymmetrical() const requires (isSquare) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            if (!compare(data[c][r], data[r][c]))
                return false;
        }
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isSkewSymmetrical() const requires (isSquare) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < ROWS; r++) {
            if (!compare(data[c][r], -data[r][c]))
                return false;
        }
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isHermitian() const requires (isSquare) {
    if constexpr (!isComplex)
        return isSymmetrical();
    else {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                if (!compare(data[c][r], std::conj(data[r][c]))) {
                    return false;
                }
            }
        }

        return true;
    }
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isSkewHermitian() const requires (isSquare) {
    if constexpr (!isComplex)
        return isSkewSymmetrical();
    else {
        for (int c = 0; c < COLUMNS; c++) {
            for (int r = 0; r < ROWS; r++) {
                if (!compare(data[c][r], -std::conj(data[r][c])))
                    return false;
            }
        }

        return true;
    }
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isPositiveDefinite() const {
    Matrix<COLUMNS, ROWS, T> ref = toRowEchelon();

    // pivots of ref are signs of eigenvalues
    for (int c = 0; c < COLUMNS; c++) {
        if (ref[c][c] <= 0)
            return false;
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isPositiveSemiDefinite() const {
    Matrix<COLUMNS, ROWS, T> ref = toRowEchelon();

    // pivots of ref are signs of eigenvalues
    for (int c = 0; c < COLUMNS; c++) {
        if (ref[c][c] < 0)
            return false;
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isNegativeDefinite() const {
    Matrix<COLUMNS, ROWS, T> ref = toRowEchelon();

    // pivots of ref are signs of eigenvalues
    for (int c = 0; c < COLUMNS; c++) {
        if (ref[c][c] >= 0)
            return false;
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isNegativeSemiDefinite() const {
    Matrix<COLUMNS, ROWS, T> ref = toRowEchelon();

    // pivots of ref are signs of eigenvalues
    for (int c = 0; c < COLUMNS; c++) {
        if (ref[c][c] > 0)
            return false;
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isUnitary() const requires (isSquare) {
    return conjugateTranspose() * *this == identity();
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isSpecialUnitary() const requires (isSquare) {
    return conjugateTranspose() * *this == identity() && compare(determinant(), 1);
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isOrthogonal() const requires (!isComplex && isSquare) {
    return transpose() * *this == identity();
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isSpecialOrthogonal() const requires (!isComplex && isSquare) {
    return transpose() * *this == identity() && compare(determinant(), 1);
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isSemiOrthogonal() const requires (!isComplex && !isSquare) {
    if constexpr (COLUMNS > ROWS) { // wide
        return multiply(transpose()) == identity();
    }
    else { // tall
        return transpose() * *this == identity();
    }
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isUpperTriangleMatrix() const requires (isSquare) {
    for (int c = 0; c < COLUMNS; c++) {
        // dont worry about out of bounds, loop wont even run if c + 1 is too big
        for (int r = c + 1; r < ROWS; r++) {
            if (!compare(data[c][r], 0))
                return false;
        }
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isLowerTriangleMatrix() const requires (isSquare) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r < c; r++) {
            if (!compare(data[c][r], 0))
                return false;
        }
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isDiagonalMatrix() const requires (isSquare) {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isUpperUnitriangularMatrix() const requires (isSquare) {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isLowerUnitriangularMatrix() const requires (isSquare) {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isStrictlyUpperTriangularMatrix() const requires (isSquare) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = c; r < ROWS; r++) {
            if (!compare(data[c][r], 0)) {
                return false;
            }
        }
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isStrictlyLowerTriangularMatrix() const requires (isSquare) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = 0; r <= c; r++) {
            if (!compare(data[c][r], 0))
                return false;
        }
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isFrobeniusMatrix() const requires (isSquare) {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isUpperHessenberg() const requires (isSquare) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = c + 2; r < ROWS; r++) {
            if (!compare(data[c][r], 0))
                return false;
        }
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isUnreducedUpperHessenberg() const requires (isSquare) {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isLowerHessenberg() const requires (isSquare) {
    for (int c = 0; c < COLUMNS; c++) {
        for (int r = c - 2; r >= 0; r--) {
            if (!compare(data[c][r], 0))
                return false;
        }
    }

    return true;
}

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isUnreducedLowerHessenberg() const requires (isSquare) {
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

template<int COLUMNS, int ROWS, is_scalar_v T>
bool Matrix<COLUMNS, ROWS, T>::isTridiagonal() const requires (isSquare) {
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