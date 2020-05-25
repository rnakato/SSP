add_library(ssp_func
  STATIC
    Mapfile.cpp
    ParseMapfile.cpp
    LibraryComplexity.cpp
    ShiftProfile.cpp
    FragmentClusterScore.cpp
  )

target_include_directories(ssp_func
	 PUBLIC ${PROJECT_SOURCE_DIR}/src/
)
target_include_directories(ssp_func
	 PRIVATE ${PROJECT_SOURCE_DIR}/src/htslib-1.10.2/
)
