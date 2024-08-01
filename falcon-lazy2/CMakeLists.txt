cmake_minimum_required(VERSION 3.10)
project(falcon-lazy2)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -Wshadow -Wundef -Werror -g3 -O0")
set (CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -Wshadow -Wundef -Werror -g3 -O0")
set (CMAKE_LD_FLAGS_DEBUG "-Wl-no-undefined")
set (CMAKE_C_FLAGS_RELEASE "-Wall -Wextra -Wshadow -Wundef -Werror -g3 -O9")
set (CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -Wshadow -Wundef -Werror -g3 -O9")
set (CMAKE_LD_FLAGS_RELEASE "-Wl-no-undefined")

set(TEST_INCS)
set(BENCHMARK_LIBS)
set(UTESTS_LIBS)

find_library(gtest NAMES gtest)
find_library(gtest_main NAMES gtest_main)
find_path(gtest_incs gtest/gtest.h)
find_library(benchmark NAMES benchmark)
find_library(benchmark_main NAMES benchmark_main)
find_path(benchmark_incs benchmark/benchmark.h)
if (gtest AND gtest_main AND gtest_incs)
    message(STATUS "Found gtest: L=${gtest},${gtest_main} I=${gtest_incs}")
    set(UTESTS_LIBS ${UTESTS_LIBS} ${gtest} ${gtest_main})
    set(TEST_INCS ${TEST_INCS} ${gtest_incs})
else()
    message(FATAL_ERROR "ERROR: Cannot Find gtest: L=${gtest},${gtest_main} I=${gtest_incs}")
endif()
if (benchmark AND benchmark_main AND benchmark_incs)
    message(STATUS "Found benchmark: L=${benchmark},${benchmark_main} I=${benchmark_incs}")
    set(BENCHMARK_LIBS ${BENCHMARK_LIBS} ${benchmark} ${benchmark_main})
    set(TEST_INCS ${TEST_INCS} ${benchmark_incs})
else()
    message(FATAL_ERROR "ERROR: Cannot Find benchmark: L=${benchmark},${benchmark_main} I=${benchmark_incs}")
endif()


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

set (TESTLIB_SRCS
        tests/testlib.cpp
        tests/testlib.h
)

set(TEST_SRCS

)

add_library(falcon STATIC ${SRCS})
target_link_libraries(falcon m)

add_library(falcon_testlib STATIC ${TESTLIB_SRCS})
target_link_libraries(falcon_testlib falcon)

add_executable(speed speed.c)
target_link_libraries(speed falcon m)

add_executable(test_falcon test_falcon.c)
target_link_libraries(test_falcon falcon m)
target_compile_options(test_falcon PRIVATE -Wno-unused)

add_executable(unittest tests/unittest.cpp)
target_link_libraries(unittest falcon falcon_testlib ${UTESTS_LIBS})
target_include_directories(unittest PRIVATE ${TEST_INCS})

add_executable(falcon_bench tests/bench_lazy_falcon.cpp)
target_link_libraries(falcon_bench falcon falcon_testlib ${BENCHMARK_LIBS})
target_include_directories(falcon_bench PRIVATE ${TEST_INCS})

