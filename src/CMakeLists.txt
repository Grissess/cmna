set(HEADER_LIST "${cmna_SOURCE_DIR}/include/cmna/cmna.h")

add_library(libcmna cmna.c ${HEADER_LIST})
target_include_directories(libcmna PUBLIC ../include)
target_link_libraries(libcmna PRIVATE LAPACK::LAPACK)
