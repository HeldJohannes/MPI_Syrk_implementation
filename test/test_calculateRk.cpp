#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include "two_d_syrk.h"
#include "utils.h"

extern void calculate_R_k(int *r_k, int k, int c);
extern int calculate_D_k(int k, int c);

class CalculateRkTest : public ::testing::Test {
protected:
    static constexpr int c3 = 3; // c ist eine Primzahl (z. B. 3, sodass P = c(c+1) = 12)
    static constexpr int P3 = c3 * (c3 + 1);
};

// Definition der statischen Konstanten außerhalb der Klasse
constexpr int CalculateRkTest::c3;
constexpr int CalculateRkTest::P3;

// Erwartete Werte aus dem Paper für c = 3
std::vector<std::vector<int>> expected_r_k_c3 = {
    {0,3,6}, {0,4,7}, {0,5,8},
    {1,3,7}, {1,4,8}, {1,5,6},
    {2,3,8}, {2,4,6}, {2,5,7},
    {0,1,2}, {3,4,5}, {6,7,8}
};

// Erwartete Werte für D_k aus dem Paper für c = 3
std::vector<int> expected_d_k_c3 = {
    -1, -1, -1,
    1, 4, 5,
    2, 6, 7,
    0, 3, 8
};

// Testfälle für verschiedene Werte von i
TEST_F(CalculateRkTest, HandlesCorrectComputation_R_k_c3) {
    for (int k = 0; k < P3; ++k) {
        intArray r_k;
        allocate_int_array(&r_k, c3);
        
        calculate_R_k(r_k.data, k, c3);

        // Prüfen, ob das Ergebnis mit der erwarteten Menge übereinstimmt
        for (int j = 0; j < c3; ++j) {
            EXPECT_EQ(r_k.data[j], expected_r_k_c3[k][j])
                << "Fehlermeldung bei k = " << k << ", Index " << j;
        }

        free_int_array(r_k);
    }
}

// Testfälle für verschiedene Werte von i
TEST_F(CalculateRkTest, HandlesCorrectComputation_D_k_c3) {
    for (int k = 0; k < P3; ++k) {
        int d_k = calculate_D_k(k, c3);

        // Prüfen, ob das Ergebnis mit der erwarteten Menge übereinstimmt
        EXPECT_EQ(d_k, expected_d_k_c3[k])
            << "Fehlermeldung bei k = " << k;
    }
}

