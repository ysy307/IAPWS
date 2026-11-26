# Define the test executable
add_executable(iapws_test iapws_test.F90)

# Link against the internal alias
# This automatically propagates the module directory (${CMAKE_BINARY_DIR}/.mod) defined in src/
target_link_libraries(iapws_test PRIVATE IAPWS::IAPWS)

# Register the test
add_test(NAME run_iapws_test COMMAND iapws_test "${PROJECT_SOURCE_DIR}/validations")
foreach(test_name saturation_batch batch_properties ordinary_properties_regression bulk_properties gibbs_derivatives region3_derivatives ice_derivatives saturation_derivatives ordinary_derivative_edges saturation_reference state_contract water_state region3_volume solver_jacobian)
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
    foreach(test_name bulk_properties batch_properties solver_jacobian)
        target_link_libraries(${test_name} PRIVATE OpenMP::OpenMP_Fortran)
        set_tests_properties(${test_name} PROPERTIES
            ENVIRONMENT "OMP_NUM_THREADS=4;OMP_DYNAMIC=FALSE;OMP_THREAD_LIMIT=4")
    endforeach()
endif()

add_executable(water_performance EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/benchmarks/water_performance.F90)
target_link_libraries(water_performance PRIVATE IAPWS::IAPWS)
if(OpenMP_Fortran_FOUND)
    target_link_libraries(water_performance PRIVATE OpenMP::OpenMP_Fortran)
endif()

if(Python3_Interpreter_FOUND)
    add_executable(edge_states ${PROJECT_SOURCE_DIR}/validations/edge_states.F90)
    target_link_libraries(edge_states PRIVATE IAPWS::IAPWS)
    add_test(NAME coolprop_reference COMMAND "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/validations/compare_water.py" "$<TARGET_FILE:edge_states>")
    set_tests_properties(coolprop_reference PROPERTIES SKIP_RETURN_CODE 77)
    add_executable(python_values ${PROJECT_SOURCE_DIR}/validations/python_values.F90)
    target_link_libraries(python_values PRIVATE IAPWS::IAPWS)
    add_test(NAME python_reference COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/python_reference.py" "$<TARGET_FILE:python_values>")
    set_tests_properties(python_reference PROPERTIES SKIP_RETURN_CODE 77)
endif()
