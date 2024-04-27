cmake_minimum_required(VERSION 3.10)
project(falcon-lazy2)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_C_FLAGS "-Wall -Wextra -Wshadow -Wundef -Werror")
set (CMAKE_CXX_FLAGS "-Wall -Wextra -Wshadow -Wundef -Werror")
set (CMAKE_LD_FLAGS "-Wl-no-undefined")

set (SRCS
        codec.c
        common.c
        config.h
        falcon.c
        falcon.h
        fft.c
        fpr.c
        fpr.h
        inner.h
        keygen.c
        rng.c
        sign.c
        #speed.c
        vrfy.c
        shake.c
)

set(TEST_SRCS
        tests/unittest.cpp
)

add_library(falcon STATIC ${SRCS})
target_link_libraries(falcon m)

add_executable(unittest ${TEST_SRCS})
target_link_libraries(unittest falcon gtest gtest_main)
