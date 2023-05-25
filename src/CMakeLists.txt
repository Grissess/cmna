set(HEADER_LIST "${cmna_SOURCE_DIR}/include/cmna/cmna.h")
set(ALL_SOURCES cmna.c ${HEADER_LIST})

add_library(cmna-static STATIC ${ALL_SOURCES})
target_include_directories(cmna-static PUBLIC $<BUILD_INTERFACE:${cmna_SOURCE_DIR}/include/>)
target_link_libraries(cmna-static PRIVATE LAPACK::LAPACK)

add_library(cmna-shared SHARED ${ALL_SOURCES})
target_include_directories(cmna-shared PUBLIC $<BUILD_INTERFACE:${cmna_SOURCE_DIR}/include/>)
target_link_libraries(cmna-shared PRIVATE LAPACK::LAPACK)


set_target_properties(cmna-static cmna-shared
	PROPERTIES OUTPUT_NAME cmna)

install(TARGETS cmna-static cmna-shared
	EXPORT cmna-config)
install(FILES ${HEADER_LIST}
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cmna/")
export(TARGETS cmna-static cmna-shared
	NAMESPACE cmna::
	FILE "${CMAKE_CURRENT_BINARY_DIR}/cmna-config.cmake")
install(EXPORT cmna-config
	NAMESPACE cmna::
	DESTINATION "${CMAKE_INSTALL_DATADIR}/cmna/cmake")
