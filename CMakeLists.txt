cmake_minimum_required(VERSION 3.0)
project(SSP CXX)

# Need C++11 due to JsonCpp nowadays
set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )

# サブディレクトリを登録
add_subdirectory(src)
add_subdirectory(common)
add_subdirectory(test)
