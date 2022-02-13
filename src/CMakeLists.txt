set(HEADER_LIST "${cmna_SOURCE_DIR}/include/cmna/cmna.h")

add_library(cmna cmna.c ${HEADER_LIST})
target_include_directories(cmna PUBLIC $<BUILD_INTERFACE:${cmna_SOURCE_DIR}/include/>)
target_link_libraries(cmna PRIVATE LAPACK::LAPACK)

install(TARGETS cmna
	EXPORT cmna-config)
install(FILES ${HEADER_LIST}
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cmna/")
export(TARGETS cmna
	NAMESPACE cmna::
	FILE "${CMAKE_CURRENT_BINARY_DIR}/cmna-config.cmake")
install(EXPORT cmna-config
	NAMESPACE cmna::
	DESTINATION "${CMAKE_INSTALL_DATADIR}/cmna/cmake")
