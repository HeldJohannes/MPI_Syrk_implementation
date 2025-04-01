#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>
#include <bits/stdc++.h>
using namespace std;
#include "MPI_Syrk_implementation.h"

class ParseInputTest : public ::testing::Test {
    protected:
        
        int rank = 0;
};

    TEST_F(ParseInputTest, ValidInput) {

        run_config config;

        const char* argv[] = {"program", "-a", "1", "-m", "10", "-n", "20", "input.txt"};
        int argc = sizeof(argv) / sizeof(argv[0]);

        int ret = parseInput(&config, argc, const_cast<char**>(argv), rank);

        EXPECT_EQ(config.algo, 1);
        EXPECT_EQ(config.m, 10);
        EXPECT_EQ(config.n, 20);
        EXPECT_STREQ(config.fileName, "input.txt");
        
        EXPECT_EQ(ret, 0);
    }

    TEST_F(ParseInputTest, MissingRequiredParameters) {

        run_config config;
        config.m = -1;
        config.n = -1;
        config.algo = -1;
        config.result_File = nullptr;
        config.c = -1;
        config.fileName = nullptr;

        const char* argv[] = {"program", "-a", "1", "-m", "10", "input.txt"};
        int argc = sizeof(argv) / sizeof(argv[0]);
    
        SCOPED_TRACE("Testing missing required parameters");
    
        // Capture the output
        testing::internal::CaptureStderr();
    
        // Run function under test
        parseInput(&config, argc, const_cast<char**>(argv), rank);
    
        // Retrieve captured stderr
        std::string output_stderr = testing::internal::GetCapturedStderr();
    
        // Ensure there is an error message
        EXPECT_FALSE(output_stderr.empty()) << "Expected error message on stderr but got none.";
    
        // Check for the specific missing parameter message
        EXPECT_NE(output_stderr.find("missing parameter n"), std::string::npos)
            << "Expected 'missing parameter [mn]' in stderr output, but got: " << output_stderr;
    }
    

    TEST_F(ParseInputTest, InvalidParameter) {

        run_config config;

        const char* argv[] = {"program", "-a", "1", "-m", "10", "-n", "20", "-x", "input.txt"};
        int argc = sizeof(argv) / sizeof(argv[0]);

        // Capture the output
        testing::internal::CaptureStderr();
    
        // Run function under test
        int ret = parseInput(&config, argc, const_cast<char**>(argv), rank);
    
        // Retrieve captured stderr
        std::string output_stderr = testing::internal::GetCapturedStderr();
    
        // Ensure there is an error message
        EXPECT_FALSE(output_stderr.empty()) << "Expected error message on stderr but got none.";
    
        // Check for the specific missing parameter message
        EXPECT_NE(output_stderr.find("wrong usage: option x doesn't exist"), std::string::npos)
            << "Expected 'wrong usage: option x doesn't exist' in stderr output, but got: " << output_stderr;
    
        EXPECT_EQ(ret, 1);
    }

    TEST_F(ParseInputTest, OptionalParameters) {

        run_config config;
        config.m = -1;
        config.n = -1;
        config.algo = -1;
        config.result_File = nullptr;
        config.c = -1;
        config.fileName = nullptr;
        

        const char* argv[] = {"program", "-a", "1", "-m", "10", "-n", "20", "-o", "result.csv", "-c", "5", "input.txt"};
        int argc = sizeof(argv) / sizeof(argv[0]);

        printf("argc: %d\n", argc);

        parseInput(&config, argc, const_cast<char**>(argv), rank);

        EXPECT_EQ(config.algo, 1);
        EXPECT_EQ(config.m, 10);
        EXPECT_EQ(config.n, 20);
        EXPECT_EQ(config.c, 5);
        EXPECT_STREQ(config.result_File, "result.csv");
        EXPECT_STREQ(config.fileName, "input.txt");
    }

    TEST_F(ParseInputTest, MissingFileName) {

        run_config config;
        config.m = -1;
        config.n = -1;
        config.algo = -1;
        config.result_File = nullptr;
        config.c = -1;
        config.fileName = nullptr;

        const char* argv[] = {"program", "-a", "1", "-m", "10", "-n", "20"};
        int argc = sizeof(argv) / sizeof(argv[0]);

        // Capture the output
        testing::internal::CaptureStdout();
        testing::internal::CaptureStderr();
    
        parseInput(&config, argc, const_cast<char**>(argv), rank);

        std::string output_stdout = testing::internal::GetCapturedStdout();
        std::string output_stderr = testing::internal::GetCapturedStderr();

        EXPECT_FALSE(output_stdout.empty() || output_stderr.empty());

        EXPECT_EQ(config.algo, 1);
        EXPECT_EQ(config.m, 10);
        EXPECT_EQ(config.n, 20);
        EXPECT_EQ(config.fileName, nullptr);

        printf("stdout: %s\n", output_stdout.c_str());
        printf("stderr: %s\n", output_stderr.c_str());

        EXPECT_EQ(output_stderr, "missing input file name --> Generate random input\n");

    }