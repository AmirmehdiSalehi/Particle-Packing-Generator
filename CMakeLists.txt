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

# Method 3: Search for local source builds (sibling directories)
if(NOT SPATIALINDEX_FOUND)
    # Common local build locations relative to project root
    set(LOCAL_BUILD_PREFIXES
        ${CMAKE_CURRENT_SOURCE_DIR}/../libspatialindex
        ${CMAKE_CURRENT_SOURCE_DIR}/../spatialindex
        ${CMAKE_CURRENT_SOURCE_DIR}/../deps/libspatialindex
        ${CMAKE_CURRENT_SOURCE_DIR}/../deps/spatialindex
        ${CMAKE_CURRENT_SOURCE_DIR}/../third-party/libspatialindex
        ${CMAKE_CURRENT_SOURCE_DIR}/../third-party/spatialindex
        ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/libspatialindex
        ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/spatialindex
        ${CMAKE_CURRENT_SOURCE_DIR}/../external/libspatialindex
        ${CMAKE_CURRENT_SOURCE_DIR}/../external/spatialindex
    )
    
    # Find include directory in local builds
    find_path(SPATIALINDEX_INCLUDE_DIR
        NAMES spatialindex/SpatialIndex.h
        PATHS ${LOCAL_BUILD_PREFIXES}
        PATH_SUFFIXES 
            include
            src
            # CMake build directories
            build/include
            build/install/include
            # Make install directories  
            install/include
            # Direct source directory
            .
        NO_DEFAULT_PATH  # Only search specified paths
        DOC "SpatialIndex include directory (local build)"
    )
    
    # Find library in local builds
    find_library(SPATIALINDEX_LIBRARY
        NAMES 
            spatialindex
            libspatialindex
            spatialindex_64
        PATHS ${LOCAL_BUILD_PREFIXES}
        PATH_SUFFIXES 
            lib
            lib64
            # CMake build directories
            build/lib
            build/lib64
            build/src
            build/install/lib
            build/install/lib64
            # Make build directories
            src/.libs
            lib/.libs
            # Install directories
            install/lib
            install/lib64
        NO_DEFAULT_PATH  # Only search specified paths
        DOC "SpatialIndex library (local build)"
    )
    
    if(SPATIALINDEX_INCLUDE_DIR AND SPATIALINDEX_LIBRARY)
        set(SPATIALINDEX_FOUND TRUE)
        set(SPATIALINDEX_INCLUDE_DIRS ${SPATIALINDEX_INCLUDE_DIR})
        set(SPATIALINDEX_LIBRARIES ${SPATIALINDEX_LIBRARY})
        message(STATUS "Found libspatialindex in local build:")
        message(STATUS "  Include: ${SPATIALINDEX_INCLUDE_DIRS}")
        message(STATUS "  Library: ${SPATIALINDEX_LIBRARIES}")
    endif()
endif()

# Method 4: Manual search in system locations
if(NOT SPATIALINDEX_FOUND)
    # Common system installation prefixes
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
        message(STATUS "Found libspatialindex in system installation:")
        message(STATUS "  Include: ${SPATIALINDEX_INCLUDE_DIRS}")
        message(STATUS "  Library: ${SPATIALINDEX_LIBRARIES}")
    endif()
endif()

# Error if not found
if(NOT SPATIALINDEX_FOUND)
    message(FATAL_ERROR "libspatialindex not found. Please install it using one of these methods:

    SYSTEM PACKAGES:
    ================
    Ubuntu/Debian: sudo apt-get install libspatialindex-dev
    CentOS/RHEL:   sudo yum install spatialindex-devel
    macOS:         brew install spatialindex
    Conda:         conda install -c conda-forge libspatialindex
    Windows:       vcpkg install spatialindex

    BUILD FROM SOURCE (will be auto-detected):
    ==========================================
    # Clone and build in parent directory:
    cd ..
    git clone https://github.com/libspatialindex/libspatialindex.git
    cd libspatialindex
    mkdir build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install
    make -j$(nproc)
    make install
    cd ../../your-project  # Back to your project
    
    # Alternative locations supported:
    # ../spatialindex/, ../deps/libspatialindex/, ../third-party/libspatialindex/
    # ../dependencies/libspatialindex/, ../external/libspatialindex/
    
    BUILD FROM SOURCE (manual):
    ===========================
    Download from: https://github.com/libspatialindex/libspatialindex
    Configure with: cmake -DCMAKE_INSTALL_PREFIX=/usr/local
    Build with: make && sudo make install")
endif()

# Find TIFF library
find_package(TIFF QUIET)
if(NOT TIFF_FOUND)
    # Try to find TIFF in local builds
    set(TIFF_LOCAL_PREFIXES
        ${CMAKE_CURRENT_SOURCE_DIR}/../libtiff
        ${CMAKE_CURRENT_SOURCE_DIR}/../tiff
        ${CMAKE_CURRENT_SOURCE_DIR}/../deps/libtiff
        ${CMAKE_CURRENT_SOURCE_DIR}/../deps/tiff
        ${CMAKE_CURRENT_SOURCE_DIR}/../third-party/libtiff
        ${CMAKE_CURRENT_SOURCE_DIR}/../third-party/tiff
        ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/libtiff
        ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/tiff
        ${CMAKE_CURRENT_SOURCE_DIR}/../external/libtiff
        ${CMAKE_CURRENT_SOURCE_DIR}/../external/tiff
    )
    
    find_path(TIFF_INCLUDE_DIR
        NAMES tiff.h
        PATHS ${TIFF_LOCAL_PREFIXES}
        PATH_SUFFIXES 
            include
            libtiff
            build/libtiff
            build/install/include
            install/include
        NO_DEFAULT_PATH
    )
    
    find_library(TIFF_LIBRARY
        NAMES tiff libtiff tiff_d libtiff_d
        PATHS ${TIFF_LOCAL_PREFIXES}
        PATH_SUFFIXES 
            lib
            lib64
            build/libtiff
            build/lib
            build/install/lib
            install/lib
        NO_DEFAULT_PATH
    )
    
    if(TIFF_INCLUDE_DIR AND TIFF_LIBRARY)
        set(TIFF_FOUND TRUE)
        set(TIFF_INCLUDE_DIRS ${TIFF_INCLUDE_DIR})
        set(TIFF_LIBRARIES ${TIFF_LIBRARY})
        message(STATUS "Found TIFF in local build:")
        message(STATUS "  Include: ${TIFF_INCLUDE_DIRS}")
        message(STATUS "  Library: ${TIFF_LIBRARIES}")
    endif()
endif()

if(NOT TIFF_FOUND)
    message(FATAL_ERROR "TIFF library not found. Please install it:
    
    SYSTEM PACKAGES:
    ================
    Ubuntu/Debian: sudo apt-get install libtiff-dev
    CentOS/RHEL:   sudo yum install libtiff-devel
    macOS:         brew install libtiff
    Conda:         conda install -c conda-forge libtiff
    Windows:       vcpkg install tiff
    
    BUILD FROM SOURCE (will be auto-detected):
    ==========================================
    cd ..
    git clone https://gitlab.com/libtiff/libtiff.git
    cd libtiff
    mkdir build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install
    make -j$(nproc)
    make install
    cd ../../your-project")
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
message(STATUS "TIP: For local source builds, use this directory structure:")
message(STATUS "  parent-dir/")
message(STATUS "    ├── libspatialindex/     (or spatialindex/)")
message(STATUS "    │   └── build/")
message(STATUS "    │       ├── include/")
message(STATUS "    │       └── lib/")
message(STATUS "    ├── libtiff/             (or tiff/)")
message(STATUS "    │   └── build/")
message(STATUS "    │       ├── include/")
message(STATUS "    │       └── lib/")
message(STATUS "    └── your-project/        (this project)")
message(STATUS "        └── build/")
message(STATUS "")