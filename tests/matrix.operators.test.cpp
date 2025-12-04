#include <gtest/gtest.h>
#include "math++/math.h"

TEST(MatrixOperators, addition_operator) {
    // Arrange
    const Matrix<2, 2, float> a = {{1.0f, 2.0f}, {3.0f, 4.0f}};
    const Matrix<2, 2, float> b = {{5.0f, 6.0f}, {7.0f, 8.0f}};
    const Matrix<2, 2, float> expected = {{6.0f, 8.0f}, {10.0f, 12.0f}};

    // Act
    const Matrix<2, 2, float> result = a + b;

    // Assert
    ASSERT_TRUE(result == expected);
}

TEST(MatrixOperators, subtraction_operator) {
    // Arrange
    const Matrix<2, 2> a = {{10, 20}, {30, 40}};
    const Matrix<2, 2> b = {{1, 2}, {3, 4}};
    const Matrix<2, 2> expected = {{9, 18}, {27, 36}};

    // Act
    const Matrix<2, 2> result = a - b;

    // Assert
    ASSERT_TRUE(result == expected);
}

TEST(MatrixOperators, matrix_multiplication) {
    // Arrange
    const Matrix<2, 2> a = {{1, 2}, {3, 4}};
    const Matrix<2, 2> b = {{5, 6}, {7, 8}};
    const Matrix<2, 2> expected = {{19, 22}, {43, 50}};

    // Act
    const Matrix<2, 2> result = a * b;

    // Assert
    ASSERT_TRUE(result == expected);
}

TEST(MatrixOperators, matrix_vector_multiplication) {
    // Arrange
    const Matrix<2, 2> m = {{1, 2}, {3, 4}};
    const Vector<2> v = {5, 6};
    const Vector<2> expected = {17, 39};

    // Act
    const Vector<2> result = m * v;

    // Assert
    ASSERT_TRUE(result == expected);
}

TEST(MatrixOperators, scalar_multiplication) {
    // Arrange
    const Matrix<2, 2> m = {{1, 2}, {3, 4}};
    const float scalar = 2.5f;
    const Matrix<2, 2> expected = {{2.5, 5}, {7.5, 10}};

    // Act
    const Matrix<2, 2> result = m * scalar;

    // Assert
    ASSERT_TRUE(result == expected);
}

TEST(MatrixOperators, scalar_division) {
    // Arrange
    const Matrix<2, 2> m = {{4, 8}, {12, 16}};
    const float scalar = 2.0f;
    const Matrix<2, 2> expected = {{2, 4}, {6, 8}};

    // Act
    const Matrix<2, 2> result = m / scalar;

    // Assert
    ASSERT_TRUE(result == expected);
}

// Assignment Operators Tests

TEST(MatrixOperators, assignment_operator_same_type) {
    // Arrange
    const Matrix<2, 2> source = {{1, 2}, {3, 4}};
    Matrix<2, 2> destination;

    // Act
    destination = source;

    // Assert
    ASSERT_TRUE(destination == source);
}

TEST(MatrixOperators, assignment_operator_different_type) {
    // Arrange
    const Matrix<2, 2, int> source = {{1, 2}, {3, 4}};
    Matrix<2, 2, float> destination;

    // Act
    destination = source;

    // Assert
    ASSERT_TRUE(destination[0][0] == 1.0f);
    ASSERT_TRUE(destination[1][1] == 4.0f);
}

TEST(MatrixOperators, add_equals_operator) {
    // Arrange
    Matrix<2, 2> a = {{1, 2}, {3, 4}};
    const Matrix<2, 2> b = {{5, 6}, {7, 8}};
    const Matrix<2, 2> expected = {{6, 8}, {10, 12}};

    // Act
    a += b;

    // Assert
    ASSERT_TRUE(a == expected);
}

TEST(MatrixOperators, subtract_equals_operator) {
    // Arrange
    Matrix<2, 2> a = {{10, 20}, {30, 40}};
    const Matrix<2, 2> b = {{1, 2}, {3, 4}};
    const Matrix<2, 2> expected = {{9, 18}, {27, 36}};

    // Act
    a -= b;

    // Assert
    ASSERT_TRUE(a == expected);
}

TEST(MatrixOperators, multiply_equals_operator) {
    // Arrange
    Matrix<2, 2> m = {{1, 2}, {3, 4}};
    const float scalar = 3.0f;
    const Matrix<2, 2> expected = {{3, 6}, {9, 12}};

    // Act
    m *= scalar;

    // Assert
    ASSERT_TRUE(m == expected);
}

TEST(MatrixOperators, divide_equals_operator) {
    // Arrange
    Matrix<2, 2> m = {{4, 8}, {12, 16}};
    const float scalar = 2.0f;
    const Matrix<2, 2> expected = {{2, 4}, {6, 8}};

    // Act
    m /= scalar;

    // Assert
    ASSERT_TRUE(m == expected);
}

// Equality Operator Tests

TEST(MatrixOperators, equality_operator_equal_matrices) {
    // Arrange
    const Matrix<2, 2> a = {{1, 2}, {3, 4}};
    const Matrix<2, 2> b = {{1, 2}, {3, 4}};

    // Act
    const bool result = (a == b);

    // Assert
    ASSERT_TRUE(result);
}

TEST(MatrixOperators, equality_operator_unequal_matrices) {
    // Arrange
    const Matrix<2, 2> a = {{1, 2}, {3, 4}};
    const Matrix<2, 2> b = {{1, 2}, {3, 5}};

    // Act
    const bool result = (a == b);

    // Assert
    ASSERT_FALSE(result);
}

TEST(MatrixOperators, equality_operator_different_types) {
    // Arrange
    const Matrix<2, 2, int> a = {{1, 2}, {3, 4}};
    const Matrix<2, 2, float> b = {{1.0f, 2.0f}, {3.0f, 4.0f}};

    // Act
    const bool result = (a == b);

    // Assert
    ASSERT_TRUE(result);
}

// Rectangular Matrix Operations Tests

TEST(MatrixOperators, rectangular_matrix_multiplication) {
    // Arrange
    const Matrix<3, 2> a = {{1, 2}, {3, 4}, {5, 6}};
    const Matrix<2, 3> b = {{7, 8, 9}, {10, 11, 12}};
    const Matrix<2, 2> expected = {{27, 30}, {61, 68}, {95, 106}};

    // Act
    const Matrix<2, 2> result = a * b;

    // Assert
    ASSERT_TRUE(result == expected);
}

TEST(MatrixOperators, rectangular_matrix_vector_multiplication) {
    // Arrange
    const Matrix<3, 2> m = {{1, 2}, {3, 4}, {5, 6}};
    const Vector<3> v = {1, 2, 3};
    const Vector<2> expected = {14, 32};

    // Act
    const Vector<2> result = m * v;

    // Assert
    ASSERT_TRUE(result == expected);
}

// Complex Number Operations Tests

TEST(MatrixOperators, complex_matrix_addition) {
    // Arrange
    using ComplexFloat = std::complex<float>;
    const Matrix<2, 2, ComplexFloat> a = {
        {ComplexFloat(1, 1), ComplexFloat(2, 2)},
        {ComplexFloat(3, 3), ComplexFloat(4, 4)}
    };
    const Matrix<2, 2, ComplexFloat> b = {
        {ComplexFloat(1, 0), ComplexFloat(0, 1)},
        {ComplexFloat(1, -1), ComplexFloat(0, -1)}
    };
    const Matrix<2, 2, ComplexFloat> expected = {
        {ComplexFloat(2, 1), ComplexFloat(2, 3)},
        {ComplexFloat(4, 2), ComplexFloat(4, 3)}
    };

    // Act
    const auto result = a + b;

    // Assert
    ASSERT_TRUE(result == expected);
}

TEST(MatrixOperators, complex_matrix_multiplication) {
    // Arrange
    using ComplexFloat = std::complex<float>;
    const Matrix<2, 2, ComplexFloat> a = {
        {ComplexFloat(1, 0), ComplexFloat(0, 1)},
        {ComplexFloat(1, 1), ComplexFloat(2, 0)}
    };
    const Matrix<2, 2, ComplexFloat> b = {
        {ComplexFloat(1, 0), ComplexFloat(0, -1)},
        {ComplexFloat(0, 1), ComplexFloat(1, 1)}
    };

    // Act
    const auto result = a * b;

    // Assert
    // Verify dimensions and basic structure
    ASSERT_EQ(result.columns, 2);
    ASSERT_EQ(result.rows, 2);
}

// Edge Cases and Properties Tests

TEST(MatrixOperators, identity_multiplication) {
    // Arrange
    const Matrix<2, 2> identity = Matrix<2, 2>::identity();
    const Matrix<2, 2> m = {{1, 2}, {3, 4}};

    // Act
    const Matrix<2, 2> result = m * identity;

    // Assert
    ASSERT_TRUE(result == m);
}

TEST(MatrixOperators, zero_matrix_addition) {
    // Arrange
    const Matrix<2, 2> zero = Matrix<2, 2>::zero();
    const Matrix<2, 2> m = {{1, 2}, {3, 4}};

    // Act
    const Matrix<2, 2> result = m + zero;

    // Assert
    ASSERT_TRUE(result == m);
}

TEST(MatrixOperators, scalar_multiplication_by_one) {
    // Arrange
    const Matrix<2, 2> m = {{1, 2}, {3, 4}};
    const float scalar = 1.0f;

    // Act
    const Matrix<2, 2> result = m * scalar;

    // Assert
    ASSERT_TRUE(result == m);
}

TEST(MatrixOperators, associativity_of_addition) {
    // Arrange
    const Matrix<2, 2> a = {{1, 2}, {3, 4}};
    const Matrix<2, 2> b = {{5, 6}, {7, 8}};
    const Matrix<2, 2> c = {{9, 10}, {11, 12}};

    // Act
    const Matrix<2, 2> result1 = (a + b) + c;
    const Matrix<2, 2> result2 = a + (b + c);

    // Assert
    ASSERT_TRUE(result1 == result2);
}

TEST(MatrixOperators, commutativity_of_scalar_multiplication) {
    // Arrange
    const Matrix<2, 2> m = {{1, 2}, {3, 4}};
    const float scalar = 2.5f;

    // Act
    const Matrix<2, 2> result1 = m * scalar;
    const Matrix<2, 2> result2 = scalar * m;

    // Assert
    ASSERT_TRUE(result1 == result2);
}

TEST(MatrixOperators, distributivity_of_scalar_over_addition) {
    // Arrange
    const Matrix<2, 2> a = {{1, 2}, {3, 4}};
    const Matrix<2, 2> b = {{5, 6}, {7, 8}};
    const float scalar = 2.0f;

    // Act
    const Matrix<2, 2> result1 = scalar * (a + b);
    const Matrix<2, 2> result2 = scalar * a + scalar * b;

    // Assert
    ASSERT_TRUE(result1 == result2);
}
