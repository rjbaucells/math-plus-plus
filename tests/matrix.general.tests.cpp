#include <gtest/gtest.h>
#include "math++/math.h"

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