#include <gtest/gtest.h>
#include "math++/math.h"

// float with float
TEST(Helper, compare_f_f_true) {
    // arrange
    constexpr float a = 1.0f;
    constexpr float b = 1.0f;

    // act
    const bool result = compare(a, b);

    // assert
    EXPECT_TRUE(result);
}

TEST(Helper, compare_f_f_false) {
    // arrange
    constexpr float a = 1.0f;
    constexpr float b = 1.1f;

    // act
    const bool result = compare(a, b);

    // assert
    EXPECT_FALSE(result);
}

// int with int
TEST(Helper, compare_i_i_true) {
    // arrange
    constexpr int a = 5;
    constexpr int b = 5;

    // act
    const bool result = compare(a, b);

    // assert
    EXPECT_TRUE(result);
}

TEST(Helper, compare_i_i_false) {
    // arrange
    constexpr int a = 5;
    constexpr int b = 7;

    // act
    const bool result = compare(a, b);

    // assert
    EXPECT_FALSE(result);
}

// complex<float> with complex<float>
TEST(Helper, compare_cf_cf_true) {
    // arrange
    constexpr std::complex<float> a{1.0f, 2.0f};
    constexpr std::complex<float> b{1.0f, 2.0f};

    // act
    const bool result = compare(a, b);

    // assert
    EXPECT_TRUE(result);
}

TEST(Helper, compare_cf_cf_false) {
    // arrange
    constexpr std::complex<float> a{1.0f, 2.0f};
    constexpr std::complex<float> b{1.0f, 2.1f};

    // act
    const bool result = compare(a, b);

    // assert
    EXPECT_FALSE(result);
}

// float with double
TEST(Helper, compare_f_d_true) {
    // arrange
    constexpr float a = 1.0f;
    constexpr double b = 1.0;

    // act
    const bool result = compare(a, b);

    // assert
    EXPECT_TRUE(result);
}

TEST(Helper, compare_f_d_false) {
    // arrange
    constexpr float a = 1.0f;
    constexpr double b = 2.0;

    // act
    const bool result = compare(a, b);

    // assert
    EXPECT_FALSE(result);
}

// int with long
TEST(Helper, compare_i_l_true) {
    // arrange
    constexpr int a = 42;
    constexpr long b = 42;

    // act
    const bool result = compare(a, b);

    // assert
    EXPECT_TRUE(result);
}

TEST(Helper, compare_i_l_false) {
    // arrange
    constexpr int a = 42;
    constexpr long b = 99;

    // act
    const bool result = compare(a, b);

    // assert
    EXPECT_FALSE(result);
}

// complex<float> with complex<double>
TEST(Helper, compare_cf_cd_true) {
    // arrange
    constexpr std::complex<float> a{1.0f, 2.0f};
    constexpr std::complex<double> b{1.0, 2.0};

    // act
    const bool result = compare(a, b);

    // assert
    EXPECT_TRUE(result);
}

TEST(Helper, compare_cf_cd_false) {
    // arrange
    std::complex<float> a{1.0f, 2.0f};
    std::complex<double> b{1.5, 2.0};

    // act
    const bool result = compare(a, b);

    // assert
    EXPECT_FALSE(result);
}

// float with float, custom precision
TEST(Helper, compare_f_f_custom_true) {
    // arrange
    constexpr float a = 1.0f;
    constexpr float b = 1.05f;
    constexpr float precision = 0.1f;

    // act
    const bool result = compare(a, b, precision);

    // assert
    EXPECT_TRUE(result);
}

TEST(Helper, compare_f_f_custom_false) {
    // arrange
    constexpr float a = 1.0f;
    constexpr float b = 1.2f;
    constexpr float precision = 0.1f;

    // act
    const bool result = compare(a, b, precision);

    // assert
    EXPECT_FALSE(result);
}

// float with double, custom precision
TEST(Helper, compare_f_d_custom_true) {
    // arrange
    constexpr float a = 1.0f;
    constexpr double b = 1.05;
    constexpr double precision = 0.1;

    // act
    const bool result = compare(a, b, precision);

    // assert
    EXPECT_TRUE(result);
}

TEST(Helper, compare_f_d_custom_false) {
    // arrange
    constexpr float a = 1.0f;
    constexpr double b = 1.2;
    constexpr double precision = 0.1;

    // act
    const bool result = compare(a, b, precision);

    // assert
    EXPECT_FALSE(result);
}

// complex<float> with complex<float>, custom precision
TEST(Helper, compare_cf_cf_custom_true) {
    // arrange
    constexpr std::complex<float> a{1.0f, 2.0f};
    constexpr std::complex<float> b{1.05f, 2.05f};
    constexpr float precision = 0.1f;

    // act
    const bool result = compare(a, b, precision);

    // assert
    EXPECT_TRUE(result);
}

TEST(Helper, compare_cf_cf_custom_false) {
    // arrange
    constexpr std::complex<float> a{1.0f, 2.0f};
    constexpr std::complex<float> b{1.2f, 2.0f};
    constexpr float precision = 0.1f;

    // act
    const bool result = compare(a, b, precision);

    // assert
    EXPECT_FALSE(result);
}

// complex<float> with complex<double>, custom precision
TEST(Helper, compare_cf_cd_custom_true) {
    // arrange
    constexpr std::complex<float> a{1.0f, 2.0f};
    constexpr std::complex<double> b{1.05, 2.05};
    constexpr double precision = 0.1;

    // act
    const bool result = compare(a, b, precision);

    // assert
    EXPECT_TRUE(result);
}

TEST(Helper, compare_cf_cd_custom_false) {
    // arrange
    constexpr std::complex<float> a{1.0f, 2.0f};
    constexpr std::complex<double> b{1.2, 2.0};
    constexpr double precision = 0.1;

    // act
    const bool result = compare(a, b, precision);

    // assert
    EXPECT_FALSE(result);
}
