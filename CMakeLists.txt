project(mhd)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Hermes2D REQUIRED)
find_package(PythonLibs REQUIRED)


link_directories(${HERMES2D_LIB})
add_executable(${PROJECT_NAME} main.cpp)

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-g -Wno-write-strings")

include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${HERMES2D_INCLUDE} ${HERMES2D_PXD_INCLUDE})

target_link_libraries(${PROJECT_NAME} hermes2d-real ${PYTHON_LIBRARIES})
