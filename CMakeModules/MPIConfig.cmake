# MPI finden und konfigurieren
if (HYDRA)
    # Check if mpi/openmpiS module is already loaded
    execute_process(
        COMMAND bash -c "module list 2>&1 | grep mpi/openmpiS"
        OUTPUT_VARIABLE MPI_MODULE_LOADED
        RESULT_VARIABLE MPI_CHECK_RESULT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # If module is not loaded, load it
    if (MPI_CHECK_RESULT)
        message(STATUS "MPI module not found, loading mpi/openmpiS...")
        execute_process(
            COMMAND bash -c "module load mpi/openmpiS"
        )
    else()
        message(STATUS "MPI module already loaded: mpi/openmpiS")
    endif()


    list(APPEND CMAKE_PREFIX_PATH "/opt/mpi/openmpi-4.1.4")
endif()

find_package(MPI REQUIRED)

if (MPI_FOUND)
    message(STATUS "MPI gefunden: ${MPI_C_LIBRARIES}")
    message(STATUS "MPI Include-Pfade: ${MPI_C_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR "MPI nicht gefunden")
endif()

# Includes setzen
include_directories(${MPI_C_INCLUDE_DIRS})
