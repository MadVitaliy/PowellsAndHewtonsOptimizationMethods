set(TargetName "PowellsMethod.Math.UnitTests")

file(GLOB_RECURSE SRC_FILES *)
message(tests)
out_module_data(TargetName "" SRC_FILES)
add_executable(${TargetName} ${SRC_FILES})

target_link_libraries(${TargetName}
        PUBLIC
        gtest
        gmock
        gtest_main
        PowellsMethod.Math
        )

add_test(NAME ${TargetName} COMMAND ${TargetName})