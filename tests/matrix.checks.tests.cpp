#include <gtest/gtest.h>
#include "math++/math.h"

TEST(MatrixChecks, positive_definite) {
    // Arrange
    constexpr Matrix<2, 2> m = {{2, -1}, {-1, 2}};
    // Act / Assert
    ASSERT_TRUE(m.isPositiveDefinite());
    ASSERT_TRUE(m.isPositiveSemiDefinite());
    ASSERT_FALSE(m.isNegativeDefinite());
    ASSERT_FALSE(m.isNegativeSemiDefinite());
}

TEST(MatrixChecks, positive_semi_definite) {
    // Arrange
    constexpr Matrix<2, 2> m = {{1, -1}, {-1, 1}};
    // Act / Assert
    ASSERT_TRUE(m.isPositiveSemiDefinite());
    ASSERT_FALSE(m.isPositiveDefinite());
    ASSERT_FALSE(m.isNegativeDefinite());
    ASSERT_FALSE(m.isNegativeSemiDefinite());
}

TEST(MatrixChecks, negative_definite) {
    // Arrange
    constexpr Matrix<2, 2> m = {{-2, -1}, {-1, -2}};
    // Act / Assert
    ASSERT_TRUE(m.isNegativeDefinite());
    ASSERT_TRUE(m.isNegativeSemiDefinite());
    ASSERT_FALSE(m.isPositiveDefinite());
    ASSERT_FALSE(m.isPositiveSemiDefinite());
}

TEST(MatrixChecks, negative_semi_definite) {
    // Arrange
    constexpr Matrix<2, 2> m = {{-1, 0}, {0, 0}};
    // Act / Assert
    ASSERT_TRUE(m.isNegativeSemiDefinite());
    ASSERT_FALSE(m.isNegativeDefinite());
    ASSERT_FALSE(m.isPositiveDefinite());
    ASSERT_FALSE(m.isPositiveSemiDefinite());
}

TEST(MatrixChecks, indefinite) {
    // Arrange
    constexpr Matrix<2, 2> m = {{1, 2}, {2, -3}};
    // Act / Assert
    ASSERT_FALSE(m.isPositiveDefinite());
    ASSERT_FALSE(m.isPositiveSemiDefinite());
    ASSERT_FALSE(m.isNegativeDefinite());
    ASSERT_FALSE(m.isNegativeSemiDefinite());
}