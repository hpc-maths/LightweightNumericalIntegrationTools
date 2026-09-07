add_executable(test_StridedView test_StridedView.cpp)
target_include_directories(test_StridedView PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_compile_options(test_StridedView PRIVATE
    $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:GNU>>:${LNIT_GNU_CXX_COMPILE_OPTIONS}>
    $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:Clang>>:${LNIT_Clang_CXX_COMPILE_OPTIONS}>
)

add_test(NAME StridedView COMMAND test_StridedView)
