# Bioinformatics 1
University of Zagreb
Faculty of Electrical Engineering and Computing
www.fer.unizg.hr/predmet/bio1

# Bamboo Filter

This project contains a C++ implementation of a Bamboo-style Cuckoo Filter.

## Prerequisites

* A C++ compiler that supports C++17 (e.g., g++, clang++)
* CMake (version 3.10 or newer)
* Make (or your chosen build system generator)

## Building the Project

To build the project, follow these steps from the root directory of the project:

1.  Create a build directory:
    ```bash
    mkdir build
    ```

2.  Navigate into the build directory:
    ```bash
    cd build
    ```

3.  Run CMake to configure the project:
    ```bash
    cmake ..
    ```
    If you want a release build (optimized), you can specify it:
    ```bash
    cmake .. -DCMAKE_BUILD_TYPE=Release
    ```
    For a debug build:
    ```bash
    cmake .. -DCMAKE_BUILD_TYPE=Debug
    ```

4.  Compile the project:
    ```bash
    make
    ```

## Running

After a successful build, the executable will be located in the `build` directory. You can run it from the `build` directory:

```bash
./BambooFilterTest