#include <gtest/gtest.h>
#include "math++/math.h"

TEST(Matrix, power_iteration) {
    const Matrix<2, 2> m = {{6, 4}, {7, 1}};

    auto [eigenVec, eigenVal] = m.powerIteration(10000);

    ASSERT_TRUE(eigenVec * eigenVal == m * eigenVec);
}

TEST(Matrix, lanczos_algorithm) {
    const Matrix<2, 2> m = {{1, 4}, {2, 3}};

    auto [t, q] = m.lanczosAlgorithm<15>();
}

TEST(Matrix, should_compile) {
    const Matrix<2, 2> m = {{1, 0}, {0, 1}};
    const Vector<2> v = {2, 2};

    m * v;
    v * m;
    m * m;
    v * v;
    m * 2;
    m / 2;
    2 * m;
    v * 2;
    v / 2;
    2 * v;
}