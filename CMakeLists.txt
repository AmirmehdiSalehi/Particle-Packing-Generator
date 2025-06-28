cmake_minimum_required(VERSION 3.14)
project(RandomPackingGenerator VERSION 1.0)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add position independent code flag
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Options
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
option(BUILD_EXAMPLES "Build example programs" ON)

# Platform-specific configurations
if(WIN32)
    # Windows-specific settings
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Find required packages with better cross-platform support
find_package(PkgConfig QUIET)

# Try to find libspatialindex using multiple methods
set(SPATIALINDEX_FOUND FALSE)

# Method 1: Try standard CMake find_package
find_package(libspatialindex QUIET)
if(libspatialindex_FOUND)
    set(SPATIALINDEX_FOUND TRUE)
    set(SPATIALINDEX_LIBRARIES ${libspatialindex_LIBRARIES})
    set(SPATIALINDEX_INCLUDE_DIRS ${libspatialindex_INCLUDE_DIRS})
    message(STATUS "Found libspatialindex using find_package")
endif()

# Method 2: Try pkg-config
if(NOT SPATIALINDEX_FOUND AND PkgConfig_FOUND)
    pkg_check_modules(SPATIALINDEX QUIET libspatialindex)
    if(SPATIALINDEX_FOUND)
        message(STATUS "Found libspatialindex using pkg-config")
    endif()
endif()

# Method 3: Manual search in common locations
if(NOT SPATIALINDEX_FOUND)
    # Common installation prefixes
    set(SEARCH_PREFIXES
        /usr
        /usr/local
        /opt/local  # MacPorts
        /opt/homebrew  # Homebrew on Apple Silicon
        $ENV{CONDA_PREFIX}  # Conda environment
        $ENV{VCPKG_ROOT}/installed/${VCPKG_TARGET_TRIPLET}  # vcpkg on Windows
    )
    
    # Find include directory
    find_path(SPATIALINDEX_INCLUDE_DIR
        NAMES spatialindex/SpatialIndex.h
        PATHS ${SEARCH_PREFIXES}
        PATH_SUFFIXES include
        DOC "SpatialIndex include directory"
    )
    
    # Find library with platform-specific names
    find_library(SPATIALINDEX_LIBRARY
        NAMES 
            spatialindex          # Standard name
            libspatialindex       # With lib prefix
            spatialindex_64       # 64-bit version
        PATHS ${SEARCH_PREFIXES}
        PATH_SUFFIXES 
            lib 
            lib64 
            lib/x86_64-linux-gnu  # Ubuntu/Debian multiarch
            lib/aarch64-linux-gnu # ARM64 Linux
        DOC "SpatialIndex library"
    )
    
    if(SPATIALINDEX_INCLUDE_DIR AND SPATIALINDEX_LIBRARY)
        set(SPATIALINDEX_FOUND TRUE)
        set(SPATIALINDEX_INCLUDE_DIRS ${SPATIALINDEX_INCLUDE_DIR})
        set(SPATIALINDEX_LIBRARIES ${SPATIALINDEX_LIBRARY})
        message(STATUS "Found libspatialindex manually:")
        message(STATUS "  Include: ${SPATIALINDEX_INCLUDE_DIRS}")
        message(STATUS "  Library: ${SPATIALINDEX_LIBRARIES}")
    endif()
endif()

# Error if not found
if(NOT SPATIALINDEX_FOUND)
    message(FATAL_ERROR "libspatialindex not found. Please install it using one of these methods:
    
    Ubuntu/Debian: sudo apt-get install libspatialindex-dev
    CentOS/RHEL:   sudo yum install spatialindex-devel
    macOS:         brew install spatialindex
    Conda:         conda install -c conda-forge libspatialindex
    Windows:       vcpkg install spatialindex
    
    Or build from source: https://github.com/libspatialindex/libspatialindex")
endif()

# Find TIFF library
find_package(TIFF REQUIRED)
if(NOT TIFF_FOUND)
    message(FATAL_ERROR "TIFF library not found. Please install it:
    
    Ubuntu/Debian: sudo apt-get install libtiff-dev
    CentOS/RHEL:   sudo yum install libtiff-devel
    macOS:         brew install libtiff
    Conda:         conda install -c conda-forge libtiff
    Windows:       vcpkg install tiff")
endif()

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${SPATIALINDEX_INCLUDE_DIRS}
    ${TIFF_INCLUDE_DIRS}
)

# Add main library
add_library(PackingGenerator SHARED
    src/PackingGenerator.cpp
    src/CApi.cpp
    src/Interface.cpp
    src/Particle.cpp
    src/Sphere.cpp
    src/VoxelGrid.cpp
)

# Link dependencies to main library
target_link_libraries(PackingGenerator
    ${SPATIALINDEX_LIBRARIES}
    ${TIFF_LIBRARIES}
)

# Platform-specific linking
if(WIN32)
    # Windows may need additional libraries
    target_link_libraries(PackingGenerator ws2_32)
endif()

# Set output name and directories
set_target_properties(PackingGenerator PROPERTIES
    OUTPUT_NAME "packing_generator"
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

# Set up examples if enabled
if(BUILD_EXAMPLES)
    add_executable(example tests/Example.cpp)
    target_link_libraries(example PackingGenerator)
    
    # Set output directory for example
    set_target_properties(example PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
endif()

# Set up Python bindings if enabled
if(BUILD_PYTHON_BINDINGS)
    # Find Python
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    
    # Find pybind11 - try different methods
    find_package(pybind11 CONFIG QUIET)
    if(NOT pybind11_FOUND)
        # Try to find via pip install
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "import pybind11; print(pybind11.get_cmake_dir())"
            OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(PYBIND11_CMAKE_DIR)
            find_package(pybind11 PATHS ${PYBIND11_CMAKE_DIR})
        endif()
    endif()
    
    if(NOT pybind11_FOUND)
        message(FATAL_ERROR "pybind11 not found. Install it using:
        pip install pybind11
        or
        conda install -c conda-forge pybind11")
    endif()
    
    # Add Python module
    pybind11_add_module(particle_packing bindings/binding.cpp)
    
    # Link with our library
    target_link_libraries(particle_packing PRIVATE
        PackingGenerator
        ${SPATIALINDEX_LIBRARIES}
        ${TIFF_LIBRARIES}
    )
    
    # Set properties for Python module
    set_target_properties(particle_packing PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
        # Ensure proper Python extension naming across platforms
        PREFIX ""
    )
    
    # Platform-specific Python module settings
    if(WIN32)
        set_target_properties(particle_packing PROPERTIES SUFFIX ".pyd")
    endif()
    
    # Install Python module
    install(TARGETS particle_packing
        LIBRARY DESTINATION ${Python3_SITELIB}
        RUNTIME DESTINATION ${Python3_SITELIB}  # For Windows .pyd files
    )
    
    # Create __init__.py for package
    file(WRITE ${CMAKE_BINARY_DIR}/python/__init__.py 
        "# Random Packing Generator Python Package\n")
    install(FILES ${CMAKE_BINARY_DIR}/python/__init__.py
        DESTINATION ${Python3_SITELIB})
endif()

# Install library and headers
install(TARGETS PackingGenerator
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# Install headers
install(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)

# Testing support
enable_testing()
if(BUILD_EXAMPLES)
    add_test(NAME BasicTest COMMAND example)
endif()

# Output configuration information
message(STATUS "")
message(STATUS "Configuration Summary:")
message(STATUS "======================")
message(STATUS "  Project:              ${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "  CMake version:        ${CMAKE_VERSION}")
message(STATUS "  System:               ${CMAKE_SYSTEM_NAME}")
message(STATUS "  C++ compiler:         ${CMAKE_CXX_COMPILER}")
message(STATUS "  C++ standard:         ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build type:           ${CMAKE_BUILD_TYPE}")
message(STATUS "  Install prefix:       ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
message(STATUS "Options:")
message(STATUS "  Build Python bindings: ${BUILD_PYTHON_BINDINGS}")
message(STATUS "  Build examples:         ${BUILD_EXAMPLES}")
message(STATUS "")
message(STATUS "Dependencies:")
message(STATUS "  SpatialIndex found:     ${SPATIALINDEX_FOUND}")
if(SPATIALINDEX_FOUND)
    message(STATUS "    Include dir:          ${SPATIALINDEX_INCLUDE_DIRS}")
    message(STATUS "    Libraries:            ${SPATIALINDEX_LIBRARIES}")
endif()
message(STATUS "  TIFF found:             ${TIFF_FOUND}")
if(TIFF_FOUND)
    message(STATUS "    Include dir:          ${TIFF_INCLUDE_DIRS}")
    message(STATUS "    Libraries:            ${TIFF_LIBRARIES}")
endif()
if(BUILD_PYTHON_BINDINGS)
    message(STATUS "  Python:                 ${Python3_VERSION}")
    message(STATUS "    Executable:           ${Python3_EXECUTABLE}")
    message(STATUS "    Site packages:        ${Python3_SITELIB}")
    message(STATUS "  pybind11 found:         ${pybind11_FOUND}")
endif()
message(STATUS "")