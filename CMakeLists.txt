cmake_minimum_required(VERSION 3.14)
project(LightweightNumericalIntegrationTools VERSION 1.0.0 LANGUAGES CXX)

# List of compiler flags
set(LNIT_GNU_CXX_COMPILE_OPTIONS
	-fopenmp-simd
	-Wundef
	-Wvarargs
	-Wall
	-Wextra
	-Winit-self
	-Wpedantic
	-Werror
	-Wconversion
	-Wuninitialized
	-Wmissing-declarations
	-Wsign-conversion
	-Wshadow
	-Wcast-align
	-Wnull-dereference
	-Wformat=2
	-flax-vector-conversions
	-pedantic
	-pedantic-errors
	-Wno-error=array-bounds
	-Wno-c99-extensions
	-Wdouble-promotion
	-Wswitch-enum
	-Wrange-loop-construct
	-Wnon-virtual-dtor
	-Wold-style-cast
	-Woverloaded-virtual
	-Wvexing-parse)

set(LNIT_Clang_CXX_COMPILE_OPTIONS
	-fopenmp-simd
	-Wundef
	-Wvarargs
	-Wall
	-Wextra
	-Winit-self
	-Wpedantic
	-Werror
	-Wconversion
	-Wuninitialized
	-Wmissing-declarations
	-Wsign-conversion
	-Wshadow
	-Wcast-align
	-Wnull-dereference
	-Wformat=2
	-flax-vector-conversions
	-pedantic
	-pedantic-errors
	-Wno-error=array-bounds
	-Wno-c99-extensions
	-Wdouble-promotion
	-Wswitch-enum
	-Wrange-loop-construct
	-Wnon-virtual-dtor
	-Wold-style-cast
	-Woverloaded-virtual
	-Wvexing-parse
	-Wno-error=switch-enum)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)

    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

# === Options ===
option(LNIT_CXX23           "Use the C++23 standard" OFF)
option(LNIT_BUILD_DEMO      "Build demo executable" OFF)
option(LNIT_NATIVE_OPT      "Enable CPU-specific optimizations (-march=native -mtune=native)" OFF)
option(LNIT_BUILD_BENCHMARK "Build benchmark executables" OFF)
option(LNIT_BUILD_DOC       "Build Doxygen documentation" OFF)
option(LNIT_BUILD_TESTS     "Build unit tests (run with ctest)" OFF)

if(LNIT_CXX23)
	set(CMAKE_CXX_STANDARD 23)
else()
	set(CMAKE_CXX_STANDARD 20)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# === Dependencies ===
find_package(fmt REQUIRED)

# === Doxygen ===
if(LNIT_BUILD_DOC)
    find_package(Doxygen REQUIRED dot)

    # Set default values for Doxyfile
    set(DOXYGEN_PROJECT_NAME "Lightweight Numerical Integration Tools (LNIT)")
    set(DOXYGEN_INPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
    set(DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/doc")

    set(DOXYGEN_FILE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${DOXYGEN_FILE} @ONLY)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_FILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
endif()

# === Demo subdir ===
if(LNIT_BUILD_DEMO)
    add_subdirectory(demo)
endif()

# === Tests subdir ===
if(LNIT_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# === Benchmark subdir ===
if(LNIT_BUILD_BENCHMARK)
    add_subdirectory(benchmark)
endif()

# === Library Target ===
set(LNIT_SRC
	src/LNIT/GaussLaguerreQuadrature.cpp
	src/LNIT/GaussHermiteQuadrature.cpp
	#
	src/LNIT/AdaptiveQuadratures/ClenshawCurtisAdaptiveQuadrature.cpp
	src/LNIT/AdaptiveQuadratures/GaussLegendreAdaptiveQuadrature.cpp
	src/LNIT/AdaptiveQuadratures/ClenshawCurtisHybridAdaptiveQuadrature.cpp
	src/LNIT/AdaptiveQuadratures/GLCCAdaptiveQuadrature.cpp
)

add_library(LNIT STATIC ${LNIT_SRC})

target_include_directories(LNIT
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(LNIT PUBLIC fmt::fmt)

target_compile_options(LNIT PRIVATE
    $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:GNU>>:${LNIT_GNU_CXX_COMPILE_OPTIONS}>
    $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:Clang>>:${LNIT_Clang_CXX_COMPILE_OPTIONS}>
)

if(LNIT_NATIVE_OPT)
    message(STATUS "Enabling native CPU optimizations (-march=native -mtune=native)")
    target_compile_options(LNIT PRIVATE -march=native -mtune=native)
endif()

# === Installation ===
include(GNUInstallDirs)

install(DIRECTORY include/LNIT DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS LNIT
    EXPORT LNITTargets
)

install(EXPORT LNITTargets
    FILE LNITTargets.cmake
    NAMESPACE LNIT::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LNIT
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/LNITConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/LNITConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/LNITConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LNIT
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/LNITConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/LNITConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LNIT
)

# === pkg-config ===
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/lnit.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/lnit.pc
    @ONLY
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lnit.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)
