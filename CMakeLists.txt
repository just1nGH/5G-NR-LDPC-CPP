cmake_minimum_required(VERSION 3.15)

project(5G-NR-LDPC-CPP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



# Library composed of two implementation files
add_library(nr_ldpc STATIC
    "Source files/nrLDPCTables.cpp"
)

file(GLOB MAIN_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/Source files/example_run.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/Source files/LDPC.cpp"
)

# Single executable target using the example runner
add_executable(example_run ${MAIN_SOURCES})

# Private include directory for the example (quotes needed because the folder name contains a space)
target_include_directories(example_run PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Source files")

# Link the example executable against the nr_ldpc library
target_link_libraries(example_run PRIVATE nr_ldpc)