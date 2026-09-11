cmake_minimum_required(VERSION 3.20)
project(orthrus VERSION 0.3.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# "Control your environment, or the defaults will control you."
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -fno-omit-frame-pointer -g -Wall -Wextra -Wpedantic -Wshadow -Wconversion")

include(GNUInstallDirs)

find_path(GMP_INCLUDE_DIR NAMES gmp.h gmpxx.h)
find_library(GMP_LIBRARY NAMES gmp)
find_library(GMPXX_LIBRARY NAMES gmpxx)

if (NOT GMP_LIBRARY OR NOT GMPXX_LIBRARY)
    message(FATAL_ERROR "GMP/GMPXX not found! Run: sudo pacman -S gmp")
endif()

# Библиотеки
add_library(orthrus SHARED src/c_api.cpp)
target_include_directories(orthrus PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    ${GMP_INCLUDE_DIR}
)
target_link_libraries(orthrus PRIVATE ${GMPXX_LIBRARY} ${GMP_LIBRARY})

add_library(orthrus_static STATIC src/c_api.cpp)
target_include_directories(orthrus_static PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    ${GMP_INCLUDE_DIR}
)
target_link_libraries(orthrus_static PRIVATE ${GMPXX_LIBRARY} ${GMP_LIBRARY})

# Исполняемые файлы
add_executable(orthrus-cli src/cli.cpp)
target_link_libraries(orthrus-cli PRIVATE orthrus_static ${GMPXX_LIBRARY} ${GMP_LIBRARY})

add_executable(c_demo tests/c_demo.c)
target_link_libraries(c_demo PRIVATE orthrus)

add_executable(orthrus_fuzz tests/fuzz_cornacchia.cpp)
target_include_directories(orthrus_fuzz PRIVATE include ${GMP_INCLUDE_DIR})
target_link_libraries(orthrus_fuzz PRIVATE ${GMPXX_LIBRARY} ${GMP_LIBRARY})

add_executable(orthrus_bench src/bench.cpp)
target_include_directories(orthrus_bench PRIVATE include ${GMP_INCLUDE_DIR})
target_link_libraries(orthrus_bench PRIVATE ${GMPXX_LIBRARY} ${GMP_LIBRARY})

add_executable(test_uint256 tests/test_uint256.cpp)
target_include_directories(test_uint256 PRIVATE include ${GMP_INCLUDE_DIR})
target_link_libraries(test_uint256 PRIVATE ${GMPXX_LIBRARY} ${GMP_LIBRARY})

# pkg-config и инсталляция
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/orthrus.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/orthrus.pc
    @ONLY
)

install(TARGETS orthrus orthrus_static orthrus-cli
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/orthrus DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/orthrus.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
