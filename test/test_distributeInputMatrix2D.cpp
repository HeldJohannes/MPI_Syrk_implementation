#include <gtest/gtest.h>
#include <mpi.h>
#include "two_d_syrk.h"

void SetUp() override {
    // Initialize MPI
    int provided;
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);
    ASSERT_EQ(provided, MPI_THREAD_MULTIPLE);

    // Initialize run_config
    s.m = 4; // Example values
    s.n = 4; // Example values
    s.c = 2; // Example values
    s.world_size = s.c * (s.c + 1);

    // Allocate input_array and rank_input
    input_array = new float[s.m * s.n];
    rank_input = new float*[s.m / (s.c * s.c)];
    for (int i = 0; i < s.m / (s.c * s.c); ++i) {
        rank_input[i] = new float[s.n / (s.c + 1)];
    }

    // Initialize input_array with example values
    for (int i = 0; i < s.m * s.n; ++i) {
        input_array[i] = static_cast<float>(i);
    }
}

void TearDown() override {
    // Free allocated memory
    delete[] input_array;
    for (int i = 0; i < s.m / (s.c * s.c); ++i) {
        delete[] rank_input[i];
    }
    delete[] rank_input;

    // Finalize MPI
    MPI_Finalize();
}

// Testfall für distribute_input_matrix_2D
TEST(DistributeInputMatrix2DTest, BasicTest) {

    // TODO: Fix this test

    int argc = 0;
    char **argv = nullptr;
    initialize_mpi(&argc, &argv);

    int rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Konfiguration für den Test
    run_config config;
    config.m = 4; // Anzahl der Zeilen
    config.n = 4; // Anzahl der Spalten
    config.c = 2; // Primzahl, sodass P = c(c+1)
    config.world_size = world_size;

    // Eingabematrix initialisieren
    float input_array[16] = {
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16
    };

    // Speicher für rank_input allokieren
    float **rank_input = (float **) calloc(config.c * (config.m / (config.c * config.c)), sizeof(float *));
    for (int i = 0; i < config.c * (config.m / (config.c * config.c)); ++i) {
        rank_input[i] = (float *) calloc(config.n, sizeof(float));
    }

    // Methode aufrufen
    distribute_input_matrix_2D(&config, rank, input_array, rank_input);

    // Überprüfen, ob die Verteilung korrekt ist
    if (rank == 0) {
        EXPECT_FLOAT_EQ(rank_input[0][0], 1);
        EXPECT_FLOAT_EQ(rank_input[0][1], 2);
        EXPECT_FLOAT_EQ(rank_input[0][2], 3);
        EXPECT_FLOAT_EQ(rank_input[0][3], 4);
    } else if (rank == 1) {
        EXPECT_FLOAT_EQ(rank_input[0][0], 5);
        EXPECT_FLOAT_EQ(rank_input[0][1], 6);
        EXPECT_FLOAT_EQ(rank_input[0][2], 7);
        EXPECT_FLOAT_EQ(rank_input[0][3], 8);
    } else if (rank == 2) {
        EXPECT_FLOAT_EQ(rank_input[0][0], 9);
        EXPECT_FLOAT_EQ(rank_input[0][1], 10);
        EXPECT_FLOAT_EQ(rank_input[0][2], 11);
        EXPECT_FLOAT_EQ(rank_input[0][3], 12);
    } else if (rank == 3) {
        EXPECT_FLOAT_EQ(rank_input[0][0], 13);
        EXPECT_FLOAT_EQ(rank_input[0][1], 14);
        EXPECT_FLOAT_EQ(rank_input[0][2], 15);
        EXPECT_FLOAT_EQ(rank_input[0][3], 16);
    }

    // Speicher freigeben
    for (int i = 0; i < config.c * (config.m / (config.c * config.c)); ++i) {
        free(rank_input[i]);
    }
    free(rank_input);

    finalize_mpi();
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}