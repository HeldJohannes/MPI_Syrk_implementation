#include <gtest/gtest.h>
#include <vector>
#include <cstring>  // For memset
#include "two_d_syrk.h"
#include "utils.h"

extern int calculate_Q_i(intArray q_i, int m_bloc_i, int c);
extern void accumulate_B_into_A(run_config *s, int k, floatArray A, floatArray B, intArray R_k, intArray Q_i, floatMatrix input);

// Sample data from the problem statement
static constexpr int C = 3;  // Example value, replace with the actual one
static constexpr int BLOCK_HEIGHT = 4;  // Example value
static constexpr int BLOCK_LENGTH = 3;  // Example value
static constexpr int M = 36;  // Example matrix height
static constexpr int N = 12;  // Example matrix width
static constexpr int k = 11;  // Example rank (this is the rank of the processor)

class AccumulateBIntoATest : public ::testing::Test {
    protected:
    run_config config;
    floatArray A, B;
    intArray R_k, Q_i;
    floatMatrix input;

    void SetUp() override {
        // Configure run_config
        config.m = M;
        config.n = N;
        config.c = C;

        // Allocate memory for A (output matrix)
        A.length = C * (C + 1) * BLOCK_HEIGHT * BLOCK_LENGTH;
        A.data = new float[A.length];
        std::memset(A.data, 0, A.length * sizeof(float));

        // Allocate memory for B (input matrix)
        // this is the matrix after the MPI_ALLTOALL operation
        B.length = C * (C + 1) * BLOCK_HEIGHT * BLOCK_LENGTH;
        B.data = new float[B.length]{
        // [0][1][2][3][4][5][6][7][8][9][10][11]
            6, 1, 9, 2, 2, 1, 8, 9, 1, 1, 3, 10,    // [0,]     (A_60)
            5, 2, 3, 4, 7, 8, 8, 5, 2, 0, 1, 9,     // [1,]     (A_70)
            7, 8, 4, 7, 4, 1, 6, 3, 3, 6, 7, 8,     // [2,]     (A_80)
            5, 7, 7, 3, 3, 9, 3, 6, 5, 6, 6, 9,     // [3,]     (A_71)
            9, 3, 10, 4, 7, 7, 2, 0, 5, 1, 9, 1,    // [4,]     (A_81)   
            7, 4, 5, 2, 3, 4, 5, 8, 8, 0, 8, 2,     // [5,]     (A_61)
            0, 1, 4, 8, 8, 8, 2, 7, 3, 0, 6, 5,     // [6,]     (A_82)
            2, 7, 9, 9, 5, 0, 7, 4, 9, 10, 1, 5,    // [7,]     (A_62)
            4, 5, 5, 6, 5, 3, 0, 1, 7, 9, 9, 6,     // [8,]     (A_72)
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // [9,]     (A_63)
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // [10,]    (A_73)
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0      // [11,]    (A_83)
        };

        input.length = N * 3;
        input.rows = N;
        input.cols = 3;
        // Allocate memory for a 2D array (array of float pointers)
        input.data = new float*[input.rows];

        float init_values[N][3] = {
            {5, 3, 9},
            {4, 4, 9},
            {6, 1, 3},
            {1, 3, 8},  // (A_63)

            {5, 0, 4},  
            {5, 5, 1},
            {3, 6, 4},
            {7, 2, 10}, // (A_73)

            {2, 4, 10},
            {7, 4, 4},
            {7, 6, 6},
            {1, 1, 3}   // (A_83)
        };

        // Allocate memory for each row and copy values
        for (int i = 0; i < input.rows; i++) {
            input.data[i] = new float[input.cols];
            std::memcpy(input.data[i], init_values[i], input.cols * sizeof(float));
        }

        // Allocate memory for R_k (Row indices)
        R_k.length = C;
        R_k.data = new int[C]{6, 7, 8};  // Mock row indices

        // Allocate memory for Q_i (Selection indices)
        Q_i.length = C + 1;
        Q_i.data = new int[C + 1];

    }

    void TearDown() override {
        delete[] A.data;
        delete[] B.data;
        delete[] R_k.data;
        delete[] Q_i.data;
        delete[] input.data;
    }
};

// ✅ Test: Check that accumulate_B_into_A correctly updates A
TEST_F(AccumulateBIntoATest, ShouldAccumulateCorrectly) {

    // Call function
    accumulate_B_into_A(&config, k, A, B, R_k, Q_i, input);

    // ✅ Expected: A should now have accumulated values from B and input

    // ✅ Expected `A` values
    floatArray expected_A;
    expected_A.length = C * (C + 1) * BLOCK_HEIGHT * BLOCK_LENGTH;
    expected_A.data = new float[N * N]{
    //  (A_60),  (A_61),  (A_62),  (A_63)
        6, 1, 9, 7, 4, 5, 2, 7, 9, 5, 3, 9,  
        2, 2, 1, 2, 3, 4, 9, 5, 0, 4, 4, 9,
        8, 9, 1, 5, 8, 8, 7, 4, 9, 6, 1, 3,
        1, 3, 10, 0, 8, 2, 10, 1, 5, 1, 3, 8,

    //  (A_70),  (A_71),  (A_72),  (A_73)
        5, 2, 3, 5, 7, 7, 4, 5, 5, 5, 0, 4,  
        4, 7, 8, 3, 3, 9, 6, 5, 3, 5, 5, 1,
        8, 5, 2, 3, 6, 5, 0, 1, 7, 3, 6, 4,
        0, 1, 9, 6, 6, 9, 9, 9, 6, 7, 2, 10,

    //  (A_80),  (A_81),  (A_82),  (A_83)
        7, 8, 4, 9, 3, 10, 0, 1, 4, 2, 4, 10,
        7, 4, 1, 4, 7, 7, 8, 8, 8, 7, 4, 4,
        6, 3, 3, 2, 0, 5, 2, 7, 3, 7, 6, 6,
        6, 7, 8, 1, 9, 1, 0, 6, 5, 1, 1, 3
    };

    // ✅ Compare A with expected values
    for (int i = 0; i < A.length; i++) {
        EXPECT_EQ(A.data[i], expected_A.data[i]) << "Mismatch at index " << i;
    }
}

