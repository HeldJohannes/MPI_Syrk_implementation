#include <gtest/gtest.h>
#include <vector>
#include <cstring>
#include "two_d_syrk.h"

extern void copy_to_2D(float *B, float **A, int k, int block_height, int block_length, int index);

class CopyTo2DTest : public ::testing::Test {
protected:
    static constexpr int c = 2;
    static constexpr int block_height = 2;
    static constexpr int block_length = 3;

    float** A = nullptr;
    float* B = nullptr;

    void SetUp() override {
        // Allocate 2D source array
        A = new float*[c * block_height];
        for (int i = 0; i < c * block_height; ++i) {
            A[i] = new float[block_length];
            for (int j = 0; j < block_length; ++j) {
                A[i][j] = static_cast<float>(i * block_length + j);
            }
        }

        // Allocate 1D destination array
        B = new float[c * (c + 1) * block_height * block_length]();
    }

    void TearDown() override {
        // Free allocated memory
        for (int i = 0; i < c * block_height; ++i) {
            delete[] A[i];
        }
        delete[] A;
        delete[] B;
    }
};

// Define static constants
constexpr int CopyTo2DTest::c;
constexpr int CopyTo2DTest::block_height;
constexpr int CopyTo2DTest::block_length;

TEST_F(CopyTo2DTest, ShouldCopySingleBlockCorrectly) {
    int k = 1;
    int index = 1;

    // Expected result
    std::vector<float> expected = {
        0, 0, 0, 0, 0, 0,
        6, 7, 8, 9, 10, 11,
        0, 0, 0, 0, 0, 0
    };

    // Perform copy
    copy_to_2D(B, A, k, block_height, block_length, index);

    // Verify result
    ASSERT_TRUE(std::equal(expected.begin(), expected.end(), B)) << "Copy failed!";
}

TEST_F(CopyTo2DTest, ShouldCopyMultipleBlocksCorrectly) {
    int index = 1;

    // Expected result after multiple copies
    std::vector<float> expected = {
        6, 7, 8, 9, 10, 11,
        0, 0, 0, 0, 0, 0,
        6, 7, 8, 9, 10, 11,
        0, 0, 0, 0, 0, 0,
        6, 7, 8, 9, 10, 11,
        0, 0, 0, 0, 0, 0
    };

    // Perform copy multiple times
    for (int i = 0; i < 3; i++) {
        copy_to_2D(B, A, i * 2, block_height, block_length, index);
    }

    // Verify result
    ASSERT_TRUE(std::equal(expected.begin(), expected.end(), B)) << "Multiple block copy failed!";
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
