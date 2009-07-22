project(mhd)

cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Hermes2D REQUIRED)
find_package(UMFPACK REQUIRED)
find_package(BLAS REQUIRED)


link_directories(${HERMES2D_LIB})
add_executable(${PROJECT_NAME} main.cpp)

include_directories(${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${HERMES2D_INCLUDE})

target_link_libraries(${PROJECT_NAME} hermes2d-real ${UMFPACK_LIBRARY}
        ${AMD_LIBRARY} ${BLAS_LIBRARIES})
