set(CMAKE_VERBOSE_MAKEFILE 1)
cmake_minimum_required(VERSION 2.8)
project(SSP CXX)

set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )
set( CMAKE_CXX_FLAGS "-std=c++11 -O2 -W -Wall")

# サブディレクトリを登録
add_subdirectory(src)
add_subdirectory(common)
add_subdirectory(test)
