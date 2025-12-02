# Define the test executable
add_executable(iapws_test iapws_test.F90)

# Link against the internal alias
# This automatically propagates the module directory (${CMAKE_BINARY_DIR}/.mod) defined in src/
target_link_libraries(iapws_test PRIVATE IAPWS::IAPWS)

# Register the test
add_test(NAME run_iapws_test COMMAND iapws_test)