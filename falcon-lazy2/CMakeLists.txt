cmake_minimum_required(VERSION 3.10)
project(falcon-lazy2)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -Wshadow -Wundef -Werror -g3 -O0")
set (CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -Wshadow -Wundef -Werror -g3 -O0")
set (CMAKE_LD_FLAGS_DEBUG "-Wl-no-undefined")
set (CMAKE_C_FLAGS_RELEASE "-Wall -Wextra -Wshadow -Wundef -Werror -g3 -O9")
set (CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -Wshadow -Wundef -Werror -g3 -O9")
set (CMAKE_LD_FLAGS_RELEASE "-Wl-no-undefined")

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

add_executable(speed speed.c)
target_link_libraries(speed falcon m)

add_executable(test_falcon test_falcon.c)
target_link_libraries(test_falcon falcon m)
target_compile_options(test_falcon PRIVATE -Wno-unused)

add_executable(unittest ${TEST_SRCS})
target_link_libraries(unittest falcon gtest gtest_main)

