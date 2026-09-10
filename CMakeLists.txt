cmake_minimum_required(VERSION 3.20)
project(klpt_engine CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

find_path(GMP_INCLUDE_DIR NAMES gmp.h gmpxx.h)
find_library(GMP_LIBRARY NAMES gmp)
find_library(GMPXX_LIBRARY NAMES gmpxx)

include_directories(${GMP_INCLUDE_DIR})
add_executable(klpt_test src/main.cpp)
target_link_libraries(klpt_test PRIVATE ${GMPXX_LIBRARY} ${GMP_LIBRARY})
