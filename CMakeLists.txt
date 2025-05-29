cmake_minimum_required(VERSION 3.10)
project(Bamboo-filter LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

set(MY_SOURCES
        src/bamboo_filter.cpp
        main.cpp
)

include_directories(src)


add_executable(BambooFilterTest ${MY_SOURCES})

message(STATUS "Konfiguracija za Bamboo-filter je završena.")
message(STATUS "Za build, koristite 'make' unutar build direktorija.")
message(STATUS "Izvršna datoteka će biti: build/BambooFilterTest")