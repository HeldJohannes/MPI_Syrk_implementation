#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include "two_d_syrk.h"

extern void calculate_Q_i(int *q_i, int m_bloc_i, int c);

class CalculateQITest : public ::testing::Test {
protected:
    static constexpr int c = 3; // c ist eine Primzahl (z. B. 3, sodass P = c(c+1) = 12)
    static constexpr int P = c * (c + 1);
};

// Erwartete Werte aus dem Bild für c = 3
std::vector<std::vector<int>> expected_q_i = {
    {0, 1, 2, 9},
    {3, 4, 5, 9},
    {6, 7, 8, 9},
    {0, 3, 6, 10},
    {1, 4, 7, 10},
    {2, 5, 8, 10},
    {0, 5, 7, 11},
    {1, 3, 8, 11},
    {2, 4, 6, 11}
};

// Testfälle für verschiedene Werte von m_bloc_i
TEST_F(CalculateQITest, HandlesCorrectComputation) {
    for (int i = 0; i < expected_q_i.size(); ++i) {
        int q_i[c + 1] = {0};  // Array für Ergebnisse
        calculate_Q_i(q_i, i, c);

        // Prüfen, ob das Ergebnis mit der erwarteten Menge übereinstimmt
        for (int j = 0; j <= c; ++j) {
            EXPECT_EQ(q_i[j], expected_q_i[i][j])
                << "Fehlermeldung bei m_bloc_i = " << i << ", Index " << j;
        }
    }
}
