cmake_minimum_required(VERSION 3.15)
project(
    "LLMDSL"
    VERSION 1.0
    LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
message(STATUS "Started CMake for ${PROJECT_NAME} v${PROJECT_VERSION}...\n")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# set vcpkg path
set(CMAKE_TOOLCHAIN_FILE "${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake")

# find packages
set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/vcpkg_installed/x64-windows")
find_package(boost_graph REQUIRED CONFIG PATHS CMAKE_PREFIX_PATH)
find_package(nlohmann_json CONFIG REQUIRED PATHS CMAKE_PREFIX_PATH)

# set include
include_directories(${CMAKE_SOURCE_DIR}/include)

# set Gurobi
include(FindGUROBI.cmake)
include_directories(${GUROBI_INCLUDE_DIRS})

# add executable
file(GLOB_RECURSE SOLVER_SOURCES src/*.cpp)
add_executable(LLMDSL ${GUI_SOURCES} ${SOLVER_SOURCES} ${VIEW_SOURCES} ${SHADER_SOURCES})

# link libraries
target_link_libraries(LLMDSL PRIVATE ${catkin_LIBRARIES} ${GUROBI_LIBRARIES})
target_link_libraries(LLMDSL PRIVATE Boost::graph nlohmann_json::nlohmann_json)