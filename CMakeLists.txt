cmake_minimum_required(VERSION 2.8)
project(SSP CXX)

if(COMMAND cmake_policy)
      cmake_policy(SET CMP0015 NEW)
endif(COMMAND cmake_policy)

set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )
set( CMAKE_CXX_FLAGS "-std=c++11 -O2 -W -Wall")

include_directories("/usr/local/include")

# サブディレクトリを登録
add_subdirectory(src)
add_subdirectory(common)
add_subdirectory(test)
