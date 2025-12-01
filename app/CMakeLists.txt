add_executable(main main.F90)

target_link_libraries(main PRIVATE iapws)

target_include_directories(main PRIVATE
    ${CMAKE_BINARY_DIR}/.mod
)

add_dependencies(main iapws)

add_test(NAME run_main COMMAND main)
