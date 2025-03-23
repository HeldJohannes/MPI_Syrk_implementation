#include <gtest/gtest.h>
#include <vector>
#include <cstring>  // For memset
#include "three_d_syrk.h"
#include "utils.h"

extern void copy_array_slice(floatArray A_i, floatArray A, int block_height, int block_length, int shift);


static constexpr int M = 3;  // Example value
static constexpr int BLOCK_LENGTH = 2;  // Example value

class TestSliceCopy : public ::testing::Test {
    protected:
        floatArray input;

    void SetUp() override {
        input.row_length = 6;
        input.length = M * input.row_length;
        input.data = new float[input.length] {
            1, 2, 3, 4, 5, 6,
            7, 8, 9, 10, 11, 12,
            13, 14, 15, 16, 17, 18
        };
    }
};

TEST_F(TestSliceCopy, ShouldCopySliceCorrectly) {
    // Allocate memory for A_i (output matrix)
    floatArray A;
    A.length = M * BLOCK_LENGTH;
    A.data = new float[A.length];
    std::memset(A.data, 0, A.length * sizeof(float));


    floatArray expected_A;
    expected_A.length = M * BLOCK_LENGTH;
    expected_A.data = new float[expected_A.length]{
        3, 4,
        9, 10,
        15, 16
    };

    //Call test method:
    copy_array_slice(A, input, M, BLOCK_LENGTH, 1 * BLOCK_LENGTH);

    for (int i = 0; i < expected_A.length; i++) {
        EXPECT_EQ(A.data[i], expected_A.data[i]) << "Mismatch at index " << i;
    }
}