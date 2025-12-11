#include <gtest/gtest.h>
#include "math++/math.h"

TEST(MatrixGeneral, inverse_expected) {
    // arrange
    constexpr Matrix<3, 3> a = {{1, 2, -1}, {-2, 0, 1}, {1, -1, 0}};
    constexpr Matrix<3, 3> expected = {{1, 1, 2}, {1, 1, 1}, {2, 3, 4}};
    // act
    const Matrix<3, 3> inverse = a.inverse();
    // assert
    ASSERT_TRUE(inverse == expected);
}

TEST(MatrixGeneral, inverse_random) {
    // arrange
    const Matrix<3, 3> a = Matrix<3, 3>::random();
    const Matrix<3, 3> identity = Matrix<3, 3>::identity();
    // act
    const Matrix<3, 3> inverse = a.inverse();
    // assert
    const Matrix<3, 3> aInverse = a * inverse;
    const Matrix<3, 3> inverseA = inverse * a;
    ASSERT_TRUE(aInverse.equals(identity, 1e-6));
    ASSERT_TRUE(inverseA.equals(identity, 1e-6));
}

TEST(MatrixGeneral, forward_substitution) {
    // arrange
    constexpr Matrix<2, 2> l = {{1, 0}, {2, 3}};
    constexpr Vector<2> b = {4, 11};
    constexpr Vector<2> expected = {4, 1};
    // act
    const Vector<2> x = l.forwardSubstitution(b);
    // assert
    ASSERT_TRUE(x == expected);
}

TEST(MatrixGeneral, backward_substitution) {
    // arrange
    constexpr Matrix<3, 3> u = {{1, -2, 1}, {0, 1, 6}, {0, 0, 1}};
    constexpr Vector<3> b = {4, -1, 2};
    constexpr Vector<3> expected = {-24, -13, 2};
    // act
    const Vector<3> x = u.backwardsSubstitution(b);
    // assert
    ASSERT_TRUE(x == expected);
}

// TEST(MatrixGeneral, solve_linear_system_inverse) {
//     // arrange
//     constexpr Matrix<3, 3> a = {{1, 1, 1}, {0, 2, 5}, {2, 5, -1}};
//     constexpr Vector<3> b = {6, -4, 27};
//     constexpr Vector<3> expected = {5, 3, -2};
//     // act
//     const Vector<3> x = a.solveLinearSystem(b, Matrix<3, 3>::LinearSystemAlgorithm::inverse);
//     // assert
//     ASSERT_TRUE(x == expected);
// }
//
// TEST(MatrixGeneral, solve_linear_system_lu) {
//     // arrange
//     constexpr Matrix<3, 3> a = {{1, 1, 1}, {0, 2, 5}, {2, 5, -1}};
//     constexpr Vector<3> b = {6, -4, 27};
//     constexpr Vector<3> expected = {5, 3, -2};
//     // act
//     const Vector<3> x = a.solveLinearSystem(b, Matrix<3, 3>::LinearSystemAlgorithm::lu_factorization);
//     // assert
//     ASSERT_TRUE(x == expected);
// }