#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include "two_d_syrk.h"
#include "utils.h"

extern int calculate_Q_i(int *q_i, int m_bloc_i, int c);

class CalculateQITest : public ::testing::Test {
protected:
    static constexpr int c3 = 3; // c ist eine Primzahl (z. B. 3, sodass P = c(c+1) = 12)
    static constexpr int P3 = c3 * (c3 + 1);

    static constexpr int c2 = 2;
    static constexpr int P2 = c2 * (c2 + 1);
};

// Definition der statischen Konstanten außerhalb der Klasse
constexpr int CalculateQITest::c3;
constexpr int CalculateQITest::P3;

constexpr int CalculateQITest::c2;
constexpr int CalculateQITest::P2;

// Erwartete Werte aus dem Bild für c = 3
std::vector<std::vector<int>> expected_q_i_c3 = {
    {0, 1, 2, 9}, {3, 4, 5, 9}, {6, 7, 8, 9},
    {0, 3, 6, 10}, {1, 4, 7, 10}, {2, 5, 8, 10},
    {0, 5, 7, 11}, {1, 3, 8, 11}, {2, 4, 6, 11}
};

// Expected values for c = 2
std::vector<std::vector<int>> expected_q_i_c2 = {
    {0, 1, 4}, {2, 3, 4},
    {0, 2, 5}, {1, 3, 5}
};

// Testfälle für verschiedene Werte von i
TEST_F(CalculateQITest, HandlesCorrectComputation_c3) {
    for (int i = 0; i < expected_q_i_c3.size(); ++i) {
        intArray q_i;
        q_i.data = (int *) malloc((c3 + 1) * sizeof(int));
        q_i.length = c3 + 1;
        calculate_Q_i(q_i, i, c3);

        // Prüfen, ob das Ergebnis mit der erwarteten Menge übereinstimmt
        for (int j = 0; j <= c3; ++j) {
            EXPECT_EQ(q_i.data[j], expected_q_i_c3[i][j])
                << "Fehlermeldung bei i = " << i << ", Index " << j;
        }
    }
}

// Testfälle für verschiedene Werte von i
TEST_F(CalculateQITest, HandlesCorrectComputation_c2) {
    for (int i = 0; i < expected_q_i_c2.size(); ++i) {
        intArray q_i;
        q_i.data = (int *) malloc((c2 + 1) * sizeof(int));
        q_i.length = c2 + 1;
        calculate_Q_i(q_i, i, c2);

        // Prüfen, ob das Ergebnis mit der erwarteten Menge übereinstimmt
        for (int j = 0; j <= c2; ++j) {
            EXPECT_EQ(q_i.data[j], expected_q_i_c2[i][j])
                << "Fehlermeldung bei i = " << i << ", Index " << j;
        }
    }
}

TEST_F(CalculateQITest, HandlesCorrectComputationForNegativeValues) {
    // Testen, ob die Funktion korrekt mit negativen Werten umgeht
    intArray q_i;
    q_i.data = (int *) malloc((c3 + 1) * sizeof(int));
    q_i.length = c3 + 1;
    EXPECT_EQ(calculate_Q_i(q_i, -1, c3), -1);
    EXPECT_EQ(calculate_Q_i(q_i, c3 * c3, c3), -1);
}

TEST_F(CalculateQITest, HandlesCorrectComputationAndCount) {
    GTEST_SKIP() << "Skip for now";
    
    int counts[c3 * c3] = {0};  // Array für die Anzahl der Elemente

    intArray q_i;
    q_i.data = (int *) malloc((c3 + 1) * sizeof(int)); // Array für Ergebnisse
    q_i.length = c3 + 1;

    for (int i = 0; i < expected_q_i_c3.size(); ++i) {

        calculate_Q_i(q_i, i, 3);

        // Prüfen, ob die Anzahl der Elemente korrekt ist
        for (int j = 0; j <= c3; ++j) {
            ASSERT_LT(q_i.data[j], P3)
                << "Fehlermeldung bei Index " << j;
            if (q_i.data[j] != -1) {
                counts[q_i.data[j]] += 1;
            }
        }
    }

    // Prüfen, ob die Anzahl der Elemente korrekt ist
    for (int i = 0; i < c3 * c3; ++i) {
        EXPECT_EQ(counts[i], c3)
            << "Fehlermeldung bei Index " << i;
    }
}
