set(TargetName "PowellsMethod.Math")

file(GLOB_RECURSE HEADER_FILES include/*)
file(GLOB_RECURSE SRC_FILES src/*)

#out_module_data(TargetName SRC_FILES HEADER_FILES)
add_library(${TargetName} ${HEADER_FILES} ${SRC_FILES})

string(TOUPPER ${TargetName} TargetNameUpper)
string(REPLACE "." "" TargetNameUpper ${TargetNameUpper})

#message(EXPORT_HEADER_DIR)
#message(${EXPORT_HEADER_DIR})
generate_export_header(${TargetName} EXPORT_MACRO_NAME ${TargetNameUpper}_API EXPORT_FILE_NAME "${EXPORT_HEADER_DIR}/${TargetName}/${TargetName}.API.h")

target_include_directories(${TargetName}
        PUBLIC
        include
        )

add_subdirectory(tests)