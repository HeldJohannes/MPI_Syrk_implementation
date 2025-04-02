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
        floatArray input, A_i;

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
    A_i.row_length = BLOCK_LENGTH;
    A_i.length = M * A_i.row_length;
    A_i.data = new float[A_i.length];
    std::memset(A_i.data, 0, A_i.length * sizeof(float));


    floatArray expected_A;
    expected_A.length = M * BLOCK_LENGTH;
    expected_A.data = new float[expected_A.length]{
        3, 4,
        9, 10,
        15, 16
    };

    // Call test method:
    EXPECT_NO_FATAL_FAILURE(copy_array_slice(A_i, input, M, BLOCK_LENGTH, 1 * BLOCK_LENGTH));

    for (int i = 0; i < expected_A.length; i++) {
        EXPECT_EQ(A_i.data[i], expected_A.data[i]) << "Mismatch at index " << i;
    }
}

// Test null pointer case for source
TEST_F(TestSliceCopy, SourceNullPointer) {
    input.data = NULL;
    EXPECT_DEATH_IF_SUPPORTED(copy_array_slice(A_i, input, M, BLOCK_LENGTH, 1 * BLOCK_LENGTH), "");
}

// Test null pointer case for destination
TEST_F(TestSliceCopy, DestinationNullPointer) {
    A_i.data = NULL;
    EXPECT_DEATH_IF_SUPPORTED(copy_array_slice(A_i, input, M, BLOCK_LENGTH, 1 * BLOCK_LENGTH), "A_i.data != NULL");
}

// Test invalid row length for source
TEST_F(TestSliceCopy, SourceInvalidRowLength) {
    input.row_length = 0;
    EXPECT_DEATH_IF_SUPPORTED(copy_array_slice(A_i, input, M, BLOCK_LENGTH, 1 * BLOCK_LENGTH), "");
}

// Test invalid row length for destination
TEST_F(TestSliceCopy, DestinationInvalidRowLength) {
    A_i.row_length = 0;
    EXPECT_DEATH_IF_SUPPORTED(copy_array_slice(A_i, input, 2, 2, 0), "");
}

// Test out-of-bounds access
TEST_F(TestSliceCopy, OutOfBoundsIndex) {
    EXPECT_DEATH_IF_SUPPORTED(copy_array_slice(A_i, input, M+1, BLOCK_LENGTH+1, 2 * BLOCK_LENGTH), "");
}