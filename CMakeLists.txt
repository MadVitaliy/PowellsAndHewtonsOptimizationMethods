cmake_minimum_required(VERSION 3.23)
set(CMAKE_CXX_STANDARD 17)

project(PowellsMethod)

function(out_module_data I_TargetName I_HEADER_FILES I_SRC_FILES )
    message(TargetName)
    message(${TargetName})
    message(HEADER_FILES)
    foreach (_header_file ${HEADER_FILES})
        message(STATUS "${_header_file}")
    endforeach()
    message(SRC_FILES)
    foreach (_src_file ${SRC_FILES})
        message(STATUS "${_src_file}")
    endforeach()
endfunction()

include(GenerateExportHeader)
include_directories(${CMAKE_BINARY_DIR}/exports)
set(EXPORT_HEADER_DIR ${CMAKE_BINARY_DIR}/exports)
set(RESOURCES_DIR ${CMAKE_BINARY_DIR}/resources)

set(CMAKE_CXX_STANDARD 17)
include(FetchContent)
FetchContent_Declare(
        googletest
        googlemock
        URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_subdirectory(src)

#add_executable(PowellsMethod main.cpp point.h #point.cpp Optimization.cpp Optimization.h)
