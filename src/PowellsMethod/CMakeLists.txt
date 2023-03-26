set(TargetName "PowellsMethod")

file(GLOB_RECURSE HEADER_FILES include/*)
file(GLOB_RECURSE SRC_FILES src/*)

out_module_data(TargetName SRC_FILES HEADER_FILES)
add_executable(${TargetName}  ${HEADER_FILES} ${SRC_FILES})

target_include_directories(${TargetName}
        PUBLIC
        include
        )
target_link_libraries(${TargetName}
        PUBLIC
        PowellsMethod.Math
        )