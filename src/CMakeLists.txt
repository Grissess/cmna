set(HEADER_LIST "${cmna_SOURCE_DIR}/include/cmna/cmna.h")

add_library(cmna cmna.c ${HEADER_LIST})
target_include_directories(cmna PUBLIC ../include)
target_link_libraries(cmna PRIVATE LAPACK::LAPACK)
