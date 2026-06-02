cmake_minimum_required(VERSION 3.12)
project(OpticalMC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Threads REQUIRED)

add_executable(OpticalMC src/OpticalMC.cpp)
target_link_libraries(OpticalMC PRIVATE Threads::Threads)
if (MSVC)
  target_compile_options(OpticalMC PRIVATE /W4)
else()
  target_compile_options(OpticalMC PRIVATE -Wall -Wextra -pedantic)
endif()
