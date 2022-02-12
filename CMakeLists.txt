cmake_minimum_required(VERSION 3.4...3.22)

project(cmna
	VERSION 0.1
	DESCRIPTION "A lightweight library for Modified Nodal Analysis"
	LANGUAGES C)

find_package(LAPACK REQUIRED)

add_subdirectory(src)
add_subdirectory(app)
