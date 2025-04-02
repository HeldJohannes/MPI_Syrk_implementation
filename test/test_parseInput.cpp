#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>
#include <bits/stdc++.h>
using namespace std;
using ::testing::HasSubstr;
#include "MPI_Syrk_implementation.h"

class ParseInputTest : public ::testing::Test {
    protected:
        run_config *config = (run_config *)malloc(sizeof(run_config));
        int rank = 0;
        char** argv;

    void SetUpArgv(const std::vector<std::string>& args) {
        argv = (char **) malloc(args.size() * sizeof(char *));
        for (size_t i = 0; i < args.size(); ++i) {
            argv[i] = (char *) malloc((args[i].size() + 1) *sizeof(char));
            strcpy(argv[i], args[i].c_str());
        }
    }

    void TearDownArgv(int argc) {
        for (int i = 0; i < argc; ++i) {
            delete[] argv[i];
        }
        delete[] argv;
    }
};

TEST_F(ParseInputTest, ValidInput) {

    std::vector<std::string> args = {"program", "-a", "1", "-m", "10", "-n", "20", "input.txt"};
    SetUpArgv(args);
    int argc = args.size();

    int ret = parseInput(config, argc, argv, rank);

    EXPECT_EQ(config->algo, 1);
    EXPECT_EQ(config->m, 10);
    EXPECT_EQ(config->n, 20);
    EXPECT_STREQ(config->fileName, "input.txt");
    
    EXPECT_EQ(ret, 0);
    TearDownArgv(argc);
}

TEST_F(ParseInputTest, MissingRequiredParameters) {

    config->m = -1;
    config->n = -1;
    config->algo = -1;
    config->result_File = nullptr;
    config->c = -1;
    config->fileName = nullptr;

    std::vector<std::string> args = {"program", "-a", "1", "-m", "10", "input.txt"};
    SetUpArgv(args);
    int argc = args.size();

    // Capture the output
    testing::internal::CaptureStderr();
    // Run function under test
    int ret = parseInput(config, argc, argv, rank);
    // Retrieve captured stderr
    std::string output_stderr = testing::internal::GetCapturedStderr();

    // Print the string for debugging:
    std::cout << "Captured stderr:\n" << output_stderr << std::endl;

    EXPECT_FALSE(output_stderr.empty()) << "Expected error message on stderr, but got none.";
    EXPECT_THAT(output_stderr, HasSubstr("missing required parameter"))
        << "Expected 'missing required parameter n' in stderr output. Actual output:\n" << output_stderr.c_str();

    // Optional: Check specific config values if applicable
    EXPECT_EQ(config->n, -1) << "Parameter 'n' should remain unchanged.";
    TearDownArgv(argc);
}


TEST_F(ParseInputTest, InvalidParameter) {

    std::vector<std::string> args = {"program", "-a", "1", "-m", "10", "-n", "20", "-x", "input.txt"};
    SetUpArgv(args);
    int argc = args.size();

    testing::internal::CaptureStderr();

    // Run function under test
    int ret = parseInput(config, argc, argv, rank);

    // Retrieve captured stderr
    std::string output_stderr = testing::internal::GetCapturedStderr();

    // Ensure there is an error message
    EXPECT_FALSE(output_stderr.empty()) << "Expected error message on stderr but got none.";

    // Check for the specific missing parameter message
    EXPECT_THAT(output_stderr, HasSubstr("invalid option -- 'x'"))
        << "Expected 'wrong usage: option x doesn't exist' in stderr output, but got: " << output_stderr;

    EXPECT_EQ(ret, 1);
    TearDownArgv(argc);
}

TEST_F(ParseInputTest, OptionalParameters) {

    config->m = -1;
    config->n = -1;
    config->algo = -1;
    config->result_File = nullptr;
    config->c = -1;
    config->fileName = nullptr;
    
    std::vector<std::string> args = {"program", "-a", "1", "-m", "10", "-n", "20", "-o", "result.csv", "-c", "5", "input.txt"};
    SetUpArgv(args);
    int argc = args.size();

    printf("argc: %d\n", argc);

    parseInput(config, argc, argv, rank);

    EXPECT_EQ(config->algo, 1);
    EXPECT_EQ(config->m, 10);
    EXPECT_EQ(config->n, 20);
    EXPECT_EQ(config->c, 5);
    EXPECT_STREQ(config->result_File, "result.csv");
    EXPECT_STREQ(config->fileName, "input.txt");
    TearDownArgv(argc);
}

TEST_F(ParseInputTest, MissingFileName) {

    config->m = -1;
    config->n = -1;
    config->algo = -1;
    config->result_File = nullptr;
    config->c = -1;
    config->fileName = nullptr;

    std::vector<std::string> args = {"program", "-a", "1", "-m", "10", "-n", "20"};
    SetUpArgv(args);
    int argc = args.size();

    // Capture the output
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    parseInput(config, argc, argv, rank);

    std::string output_stdout = testing::internal::GetCapturedStdout();
    std::string output_stderr = testing::internal::GetCapturedStderr();

    EXPECT_FALSE(output_stdout.empty() || output_stderr.empty());

    EXPECT_EQ(config->algo, 1);
    EXPECT_EQ(config->m, 10);
    EXPECT_EQ(config->n, 20);
    EXPECT_EQ(config->fileName, nullptr);

    printf("stdout: %s\n", output_stdout.c_str());
    printf("stderr: %s\n", output_stderr.c_str());

    EXPECT_THAT(output_stderr, HasSubstr("missing input file name --> Generate random input\n"));
    TearDownArgv(argc);
}