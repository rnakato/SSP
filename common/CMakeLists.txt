add_library(ssp_common
  STATIC
 BoostOptions.cpp statistics.cpp util.cpp gzstream.C
  )

target_include_directories(ssp_func PUBLIC ${PROJECT_SOURCE_DIR}/common )
