add_executable(cmna-bin main.c)
target_link_libraries(cmna-bin PRIVATE cmna-static)
set_target_properties(cmna-bin 
	PROPERTIES OUTPUT_NAME cmna)

install(TARGETS cmna-bin)
