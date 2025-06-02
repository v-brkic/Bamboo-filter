 # Bioinformatics 1
University of Zagreb
Faculty of Electrical Engineering and Computing
www.fer.unizg.hr/predmet/bio1

# Bamboo Filter Implementation

This project contains a C++ implementation of a Cuckoo-style filter with dynamic resizing capabilities, inspired by the principles of Bamboo Filters.

## Reference and Inspiration

This implementation draws inspiration from and aims to address challenges discussed in the following academic paper:

* Wang, H., Geng, Z., Li, Q., Zhang, Y., Liu, X. S., Jin, R., Jiang, H., Wu, X., & Cong, J. (2022). **Bamboo Filters: Make Resizing Smooth and Adaptive**. In *2022 IEEE 38th International Conference on Data Engineering (ICDE)*  IEEE.
    * DOI: [10.1109/ICDE53745.2022.00241](https://doi.org/10.1109/TNET.2024.3403997)
    * IEEE Xplore: [Link to paper](https://ieeexplore.ieee.org/abstract/document/9835626)

The reference C++ implementation by the original authors, which was used for comparison in this project, can be found on GitHub:
* **Original BambooFilters GitHub Repository:** [https://github.com/wanghanchengchn/bamboofilters](https://github.com/wanghanchengchn/bamboofilters)

While our implementation utilizes a global table rebuild for dynamic resizing (doubling capacity when a load threshold is met) rather than the incremental, segment-by-segment splitting detailed in the original Bamboo Filter paper, the core objectives of achieving dynamic capacity, maintaining low false positive rates, and ensuring efficient set membership testing are shared. Our approach focuses on correctness with full hash-based rebuilding and leverages 16-bit fingerprints for enhanced accuracy.

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

After a successful build, the executable will be located in the build directory. You can run it from the build directory:

```bash
./MyBambooFilterTest
