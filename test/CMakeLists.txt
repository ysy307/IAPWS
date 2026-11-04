# Define the test executable
add_executable(iapws_test iapws_test.F90)

# Link against the internal alias
# This automatically propagates the module directory (${CMAKE_BINARY_DIR}/.mod) defined in src/
target_link_libraries(iapws_test PRIVATE IAPWS::IAPWS)

# Register the test
add_test(NAME run_iapws_test COMMAND iapws_test "${PROJECT_SOURCE_DIR}/validations")
foreach(test_name gibbs_derivatives region3_derivatives ice_derivatives saturation_derivatives ordinary_derivative_edges saturation_reference)
    add_executable(${test_name} ${test_name}.F90)
    target_link_libraries(${test_name} PRIVATE IAPWS::IAPWS)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME validation_harness COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/validation_harness.py"
        "$<TARGET_FILE:iapws_test>" "${PROJECT_SOURCE_DIR}/validations")
endif()
