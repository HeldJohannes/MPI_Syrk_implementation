#include <gtest/gtest.h>
#include <vector>
#include <cstring>  // For memset
#include "two_d_syrk.h"
#include "utils.h"

extern void copy_to_d(double *A_i, float *A, int c, int block_height, int block_length, int index, _Bool trans);
extern void copy_to_f(float *A_i, float *A, int c, int block_height, int block_length, int index, _Bool trans);

extern void copy_to_1D(float *dest, float *source, int r_pos, int w_pos, int block_height, int block_length, int row_length, int w_offset);

// Sample data from the problem statement
static constexpr int C = 3;  // Example value, replace with the actual one
static constexpr int BLOCK_HEIGHT = 4;  // Example value
static constexpr int BLOCK_LENGTH = 3;  // Example value
static constexpr int M = 36;  // Example matrix height
static constexpr int N = 12;  // Example matrix width

class TestArrayCopy : public ::testing::Test {
    protected:
        run_config config;
        doubleArray A;
        floatArray B;
        floatMatrix input;

    void SetUp() override {
        // Configure run_config
        config.m = M;
        config.n = N;
        config.c = C;

        B.row_length = BLOCK_HEIGHT * BLOCK_LENGTH;
        B.length = C * (C + 1) * B.row_length;
        B.data = new float[B.length]{
            6, 1, 9, 2, 2, 1, 8, 9, 1, 1, 3, 10, // (A_60)
            7, 4, 5, 2, 3, 4, 5, 8, 8, 0, 8, 2,  // (A_61)
            2, 7, 9, 9, 5, 0, 7, 4, 9, 10, 1, 5, // (A_62)
            5, 3, 9, 4, 4, 9, 6, 1, 3, 1, 3, 8,  // (A_63)

            5, 2, 3, 4, 7, 8, 8, 5, 2, 0, 1, 9,  // (A_70)
            5, 7, 7, 3, 3, 9, 3, 6, 5, 6, 6, 9,  // (A_71)
            4, 5, 5, 6, 5, 3, 0, 1, 7, 9, 9, 6,  // (A_72)
            5, 0, 4, 5, 5, 1, 3, 6, 4, 7, 2, 10, // (A_73)

            7, 8, 4, 7, 4, 1, 6, 3, 3, 6, 7, 8,  // (A_80)
            9, 3, 10, 4, 7, 7, 2, 0, 5, 1, 9, 1, // (A_81)
            0, 1, 4, 8, 8, 8, 2, 7, 3, 0, 6, 5,  // (A_82)
            2, 4, 10, 7, 4, 4, 7, 6, 6, 1, 1, 3  // (A_83)
        };
    }

    void TearDown() override {
        
    }
};

TEST_F(TestArrayCopy, ShouldCopyToDCorrectly) {
    // Allocate memory for A_i (output matrix)
    A.length = BLOCK_HEIGHT * N;
    A.data = new double[A.length];
    std::memset(A.data, 0, A.length * sizeof(double));

    doubleArray expected_A;
    expected_A.length = BLOCK_HEIGHT * N;
    expected_A.data = new double[expected_A.length]{
    // (A_60)       (A_61)       (A_62)       (A_63)
        6, 1, 9,    7, 4, 5,     2, 7, 9,     5, 3, 9,   
        2, 2, 1,    2, 3, 4,     9, 5, 0,     4, 4, 9,
        8, 9, 1,    5, 8, 8,     7, 4, 9,     6, 1, 3, 
        1, 3, 10,   0, 8, 2,     10, 1, 5,    1, 3, 8
    };

    copy_to_d(A.data, B.data, C, BLOCK_HEIGHT, BLOCK_LENGTH, 0, false);

    for (int i = 0; i < A.length; i++) {
        EXPECT_EQ(A.data[i], expected_A.data[i]) << "Mismatch at index " << i;
    }
}

TEST_F(TestArrayCopy, ShouldCopyToDTransposeCorrectly) {

    // Allocate memory for A_i (output matrix)
    A.length = BLOCK_HEIGHT * N;
    A.data = new double[A.length];
    std::memset(A.data, 0, A.length * sizeof(double));

    doubleArray expected_A;
    expected_A.length = BLOCK_HEIGHT * N;
    expected_A.data = new double[expected_A.length]{
        6, 2, 8, 1,
        1, 2, 9, 3,
        9, 1, 1, 10,    // ((A_60)^T)

        7, 2, 5, 0,     
        4, 3, 8, 8,     
        5, 4, 8, 2,     // ((A_61)^T)

        2, 9, 7, 10,
        7, 5, 4, 1,     
        9, 0, 9, 5,     // ((A_62)^T)

        5, 4, 6, 1,
        3, 4, 1, 3,
        9, 9, 3, 8      // ((A_63)^T)
    };

    copy_to_d(A.data, B.data, C, BLOCK_HEIGHT, BLOCK_LENGTH, 0, true);

    for (int i = 0; i < A.length; i++) {
        EXPECT_EQ(A.data[i], expected_A.data[i]) << "Mismatch at index " << i;
    }

}


TEST_F(TestArrayCopy, shouldCopyCorrectly1D) {

    floatArray result;
    result.length = BLOCK_HEIGHT * BLOCK_HEIGHT;
    result.data = new float[result.length] {
        226, 341, 238, 159,
        311, 362, 360, 236,
        193, 259, 274, 168,
        203, 316, 238, 222
    };

    floatArray result_2;
    result_2.length = BLOCK_HEIGHT * BLOCK_HEIGHT;
    result_2.data = new float[result_2.length] {
        22, 34, 23, 15,
        31, 36, 36, 23,
        19, 25, 27, 16,
        20, 31, 23, 22
    };

    // ✅ Expected values
    floatArray expected_result;
    expected_result.length = M * M;
    expected_result.data = new float[expected_result.length] {
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,   0, 0, 0, 0,      0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,   0, 0, 0, 0,      0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,   0, 0, 0, 0,      0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,   0, 0, 0, 0,      0, 0, 0, 0,     0, 0, 0, 0,
        
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,
        
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    22, 34, 23, 15,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    31, 36, 36, 23,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    19, 25, 27, 16,     0, 0, 0, 0,     0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    20, 31, 23, 22,     0, 0, 0, 0,     0, 0, 0, 0,
        
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    22, 34, 23, 15,    226, 341, 238, 159,   0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    31, 36, 36, 23,    311, 362, 360, 236,   0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    19, 25, 27, 16,    193, 259, 274, 168,   0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,    20, 31, 23, 22,    203, 316, 238, 222,   0, 0, 0, 0
    };

    floatArray destination;
    destination.length = M * M;
    destination.data = new float[destination.length];
    std::memset(destination.data, 0, destination.length * sizeof(float));

    printf("run copy_to_1D:\n");

    int r_pos = 0;
    int w_pos = 8 * C * C;
    int w_offset = 7 * BLOCK_HEIGHT;

    copy_to_1D(destination.data, result.data, r_pos, w_pos, BLOCK_HEIGHT, BLOCK_HEIGHT, M, w_offset);

    w_offset = 6 * BLOCK_HEIGHT;

    copy_to_1D(destination.data, result_2.data, r_pos, w_pos, BLOCK_HEIGHT, BLOCK_HEIGHT, M, w_offset);

    w_pos = 7 * C * C;
    copy_to_1D(destination.data, result_2.data, r_pos, w_pos, BLOCK_HEIGHT, BLOCK_HEIGHT, M, w_offset);

    printArray(destination, M, M, stdout);

    for (int i = 0; i < destination.length; i++) {
        EXPECT_EQ(destination.data[i], expected_result.data[i]) << "Mismatch at index " << i;
    }
}


TEST_F(TestArrayCopy, ShouldCopyCorrectly1Dv2) {

    floatArray input;
    input.length = C * (C + 1) * BLOCK_HEIGHT * BLOCK_LENGTH;
    input.data = new float[input.length]{
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

    floatArray result;
    // Allocate memory for A (output matrix)
    result.length = C * (C + 1) * BLOCK_HEIGHT * BLOCK_LENGTH;
    result.data = new float[result.length];
    std::memset(result.data, 0, result.length * sizeof(float));

    // ✅ Expected values
    floatArray expected;
    expected.length = C * (C + 1) * BLOCK_HEIGHT * BLOCK_LENGTH;
    expected.data = new float[expected.length]{
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // (A_60)
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // (A_61)
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // (A_62)
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // (A_63)

        5, 2, 3, 4, 7, 8, 8, 5, 2, 0, 1, 9,  // (A_70)
        5, 7, 7, 3, 3, 9, 3, 6, 5, 6, 6, 9,  // (A_71)
        4, 5, 5, 6, 5, 3, 0, 1, 7, 9, 9, 6,  // (A_72)
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // (A_73)

        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // (A_80)
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // (A_81)
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // (A_82)
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  // (A_83)
    };

    int read_position = 1;
    int write_position = 0 + ((C + 1) * 1);

    copy_to_1D(result.data, input.data, read_position, write_position, BLOCK_HEIGHT, BLOCK_LENGTH, BLOCK_LENGTH, 0); // A_70


    read_position = 3;
    write_position = 1 + ((C + 1) * 1);

    copy_to_1D(result.data, input.data, read_position, write_position, BLOCK_HEIGHT, BLOCK_LENGTH, BLOCK_LENGTH, 0); // A_71


    read_position = 8;
    write_position = 2 + ((C + 1) * 1);

    copy_to_1D(result.data, input.data, read_position, write_position, BLOCK_HEIGHT, BLOCK_LENGTH, BLOCK_LENGTH, 0); // A_72




    // ✅ Compare A with expected values
    for (int i = 0; i < result.length; i++) {
        EXPECT_EQ(result.data[i], expected.data[i]) << "Mismatch at index " << i;
    }
}
