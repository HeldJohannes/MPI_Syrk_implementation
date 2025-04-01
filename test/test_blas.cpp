#include <gtest/gtest.h>
#include <vector>
#include <cstring>  // For memset
#include <cblas.h>
#include <utils.h>

#ifdef USE_CBLAS_64
    #define CBLAS_SYRK cblas_ssyrk64_
    #define CBLAS_DGEMM cblas_dgemm64_
#else
    #define CBLAS_SYRK cblas_ssyrk
    #define CBLAS_DGEMM cblas_dgemm
#endif

static constexpr int N = 12;  // Number of columns in A, number of rows in B

class TestBlas : public ::testing::Test {
    protected:
        void SetUp() override {
            
        }
};

TEST_F(TestBlas, MatrixMultiplication) {

    // Matrix dimensions
    int block_height = 4;  

    std::vector<double> C, expected_C;

    // Initialize matrices
    doubleArray A;
    A.length = block_height * N;
    A.data = new double[A.length] {
        9, 2, 5, 9, 4, 4, 10, 2, 6, 3, 3, 9, 
        10, 9, 6, 2, 9, 5, 4, 1, 1, 6, 6, 2, 
        10, 8, 8, 5, 3, 6, 3, 2, 8, 8, 0, 2, 
        1, 0, 4, 8, 4, 6, 6, 3, 5, 2, 5, 4
    }; // 4 x 12

    doubleArray B;
    B.length = block_height * N;
    B.data = new double[B.length]{
        5, 6, 8, 2, 
        6, 1, 3, 6, 
        4, 2, 4, 1, 
        2, 2, 1, 5, 
        10, 2, 2, 1, 
        2, 4, 6, 4, 
        6, 4, 6, 2, 
        9, 7, 7, 6, 
        2, 0, 3, 6, 
        5, 2, 10, 4, 
        3, 9, 9, 8, 
        8, 1, 2, 8
    }; // 12 x 4

    // Result matrix C (4x4), initialized to zero
    // Result matrix C (2x2), initialized to zero
    C.assign(block_height * block_height, 0.0);

    // Expected result: C = A * B
    expected_C = {
        329,  204,  306,  276,
        331,  214,  333,  227,
        290,  168,  323,  242,
        219,  160,  221,  214
    }; 

    CBLAS_DGEMM(
        CblasRowMajor, 
        CblasNoTrans, 
        CblasNoTrans,
        block_height,
        block_height, 
        N,
        1.0, 
        A.data, 
        N, 
        B.data, 
        block_height,
        0.0, 
        C.data(), 
        block_height
    );

    // Verify each element in C matches expected_C within a tolerance
    for (int i = 0; i < C.size(); i++) {
        EXPECT_EQ(C[i], expected_C[i]) << "Mismatch at index " << i;
    }
}


TEST_F(TestBlas, MatrixMultiplication_simple) {

    int block_height, k;
    doubleArray A, B;
    std::vector<double> C, expected_C;

    block_height = 2;  
    k = 3; // Number of columns in A, number of rows in B

    // Initialize matrices
    A.data = new double[6]{
        1.0, 2.0, 3.0, 
        4.0, 5.0, 6.0
    }; // 2x3 matrix

    B.data = new double[6]{
        7.0, 8.0, 
        9.0, 10.0, 
        11.0, 12.0
    }; // 3x2 matrix

    // Result matrix C (2x2), initialized to zero
    C.assign(block_height * block_height, 0.0);

    // Expected result: C = A * B
    expected_C = {
        58.0,  64.0,  // [1*7 + 2*9 + 3*11, 1*8 + 2*10 + 3*12]
        139.0, 154.0  // [4*7 + 5*9 + 6*11, 4*8 + 5*10 + 6*12]
    }; 

    CBLAS_DGEMM(
        CblasRowMajor, CblasNoTrans, CblasNoTrans,
        block_height, block_height, k,
        1.0, A.data, k, B.data, block_height,
        0.0, C.data(), block_height
    );

    // Verify each element in C matches expected_C within a tolerance
    for (size_t i = 0; i < C.size(); ++i) {
        ASSERT_NEAR(C[i], expected_C[i], 1e-6) << "Mismatch at index " << i;
    }
}

TEST_F(TestBlas, MatrixSYRK) {

    int block_height, k;
    floatArray A;
    std::vector<float> C, expected_C;

    block_height = 4;  
    k = 12; // Number of columns in A

    // Initialize matrices
    A.data = new float[block_height * k]{
        0,  3,  8,  8,  6,  9,  1,  9,  8,  1,  1,  2,
        4,  9,  1,  9,  9,  1,  8,  9,  7,  9,  9,  8,
        10, 10, 2,  9,  6,  4,  3,  6,  7,  1,  5,  2,
        9,  6,  2,  7,  7,  3,  9,  1,  2,  9,  9,  1
    }; // 4x12 matrix

    // Result matrix C (2x2), initialized to zero
    C.assign(block_height * block_height, 0.0);

    // Expected result: C = A * A^T
    expected_C = {
        406,  0,  0,  0,
        349,  681,  0,  0,
        313,  468,  461,  0,
        213,  486,  374,  477
    }; 

    CBLAS_SYRK(
        CblasRowMajor,          // Row/column order
        CblasLower,             //Upper or Lower triangle of C
        CblasNoTrans,       // How matrix A is to be transposed
        block_height,                      // Number of rows and columns in matrix C
        k,                      // Number of columns of the matrix A if it is not transposed, and number of rows otherwise.
        1.0f,                   // The factor of matrix A
        A.data,                 
        k,                      //
        0.0f,
        C.data(),
        block_height
    );

    // Verify each element in C matches expected_C within a tolerance
    for (size_t i = 0; i < C.size(); ++i) {
        ASSERT_NEAR(C[i], expected_C[i], 1e-6) << "Mismatch at index " << i;
    }
}


TEST_F(TestBlas, MatrixSYRK_simple) {

    int block_height, k;
    floatArray A;
    std::vector<float> C, expected_C;

    block_height = 2;  
    k = 3; // Number of columns in A, number of rows in B

    // Initialize matrices
    A.data = new float[6]{
        1.0, 2.0, 3.0, 
        4.0, 5.0, 6.0
    }; // 2x3 matrix

    // Result matrix C (2x2), initialized to zero
    C.assign(block_height * block_height, 0.0);

    // Expected result: C = A * B
    expected_C = {
        14.0, 0.0,    // [1*1 + 2*2 + 3*3,       -        ]
        32.0, 77.0    // [1*4 + 2*5 + 3*6, 4*4 + 5*5 + 6*6]
    }; 

    CBLAS_SYRK(
        CblasRowMajor,
        CblasLower,
        CblasConjNoTrans,
        block_height,
        k,
        1.0f,
        A.data,
        k,
        0.0f,
        C.data(),
        block_height
    );

    // Verify each element in C matches expected_C within a tolerance
    for (size_t i = 0; i < C.size(); ++i) {
        ASSERT_NEAR(C[i], expected_C[i], 1e-6) << "Mismatch at index " << i;
    }
}

