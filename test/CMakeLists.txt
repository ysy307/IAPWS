# Define the test executable
add_executable(iapws_test iapws_test.F90)

# Link against the internal alias
# This automatically propagates the module directory (${CMAKE_BINARY_DIR}/.mod) defined in src/
target_link_libraries(iapws_test PRIVATE IAPWS::IAPWS)

# Register the test
add_test(NAME run_iapws_test COMMAND iapws_test "${PROJECT_SOURCE_DIR}/validations")
foreach(test_name saturation_batch batch_properties ordinary_properties_regression bulk_properties gibbs_derivatives region3_derivatives ice_derivatives saturation_derivatives ordinary_derivative_edges saturation_reference)
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

# Manual timing only: no unstable wall-clock assertions in CTest.
if(EXISTS "${PROJECT_SOURCE_DIR}/benchmarks/bulk_benchmark.F90")
    add_executable(bulk_benchmark EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/benchmarks/bulk_benchmark.F90)
    target_link_libraries(bulk_benchmark PRIVATE IAPWS::IAPWS)
endif()
if(EXISTS "${PROJECT_SOURCE_DIR}/benchmarks/kernel_benchmark.F90")
    add_executable(kernel_benchmark EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/benchmarks/kernel_benchmark.F90)
    target_link_libraries(kernel_benchmark PRIVATE IAPWS::IAPWS)
endif()

find_package(OpenMP QUIET COMPONENTS Fortran)

if(EXISTS "${PROJECT_SOURCE_DIR}/benchmarks/batch_benchmark.F90")
    add_executable(batch_benchmark EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/benchmarks/batch_benchmark.F90)
    target_link_libraries(batch_benchmark PRIVATE IAPWS::IAPWS)
endif()

add_executable(batch_shapes batch_shapes.F90)
target_link_libraries(batch_shapes PRIVATE IAPWS::IAPWS)
if(Python3_Interpreter_FOUND)
    add_test(NAME batch_shapes COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/batch_shapes.py" "$<TARGET_FILE:batch_shapes>")
endif()

add_executable(batch_extension batch_extension.F90)
set_target_properties(batch_extension PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/batch_extension_mod")
target_link_libraries(batch_extension PRIVATE IAPWS::IAPWS)
add_test(NAME batch_extension COMMAND batch_extension)

if(OpenMP_Fortran_FOUND)
    foreach(test_name bulk_properties batch_properties)
        target_link_libraries(${test_name} PRIVATE OpenMP::OpenMP_Fortran)
        set_tests_properties(${test_name} PROPERTIES
            ENVIRONMENT "OMP_NUM_THREADS=4;OMP_DYNAMIC=FALSE;OMP_THREAD_LIMIT=4")
    endforeach()
endif()
