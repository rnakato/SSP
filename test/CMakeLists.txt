set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

if(APPLE)
  set(BOOST_LIBS "-lpthread -lboost_thread-mt -lboost_program_options -lboost_system -lboost_filesystem")
else()
  set(BOOST_LIBS "-lpthread -lboost_thread -lboost_program_options -lboost_system -lboost_filesystem")
endif()

#find_package(Boost 1.53.0 REQUIRED thread program_options system filesystem)
#include_directories(${Boost_INCLUDE_DIRS})

add_executable(ssp ssp_main.cpp)

link_directories(ssp "/usr/local/lib")

target_link_libraries(ssp
                      ssp_func
                      ssp_common
                      ${PROJECT_SOURCE_DIR}/src/htslib-1.10.2/libhts.a
                      ${BOOST_LIBS}
                      -lcurl -llzma -lbz2 -lz
)
