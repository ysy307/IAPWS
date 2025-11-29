cmake_minimum_required(VERSION 3.25)
project(IAPWS LANGUAGES Fortran)

add_executable(main src/main.f90)

enable_testing()

add_test(NAME run_main COMMAND main)