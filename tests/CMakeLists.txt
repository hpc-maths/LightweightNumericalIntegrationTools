# Each test is a small framework-free executable registered with ctest.
# LNIT_TESTING exposes the quadrature tables (nodes and weights) to the tests.
function(lnit_add_test name)
    add_executable(${name} ${name}.cpp)
    target_include_directories(${name} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_compile_definitions(${name} PRIVATE LNIT_TESTING)
    target_compile_options(${name} PRIVATE
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:GNU>>:${LNIT_GNU_CXX_COMPILE_OPTIONS}>
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:Clang>>:${LNIT_Clang_CXX_COMPILE_OPTIONS}>
    )
    if(ARGN)
        target_link_libraries(${name} PRIVATE ${ARGN})
    endif()
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES TIMEOUT 30)
endfunction()

lnit_add_test(test_StridedView)
lnit_add_test(test_ClenshawCurtisTables   LNIT)
lnit_add_test(test_AdaptiveQuadratureBase LNIT)
