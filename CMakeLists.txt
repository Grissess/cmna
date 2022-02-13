cmake_minimum_required(VERSION 3.4...3.22)

project(cmna
	VERSION 0.1
	DESCRIPTION "A lightweight library for Modified Nodal Analysis"
	LANGUAGES C)

include(GNUInstallDirs)

find_package(LAPACK REQUIRED)

add_subdirectory(src)
add_subdirectory(app)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	${CMAKE_CURRENT_BINARY_DIR}/cmna-config-version.cmake
	COMPATIBILITY SameMajorVersion)
