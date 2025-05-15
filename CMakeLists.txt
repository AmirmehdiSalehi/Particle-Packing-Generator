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

# Find required packages
find_package(SpatialIndex REQUIRED)
find_package(TIFF REQUIRED)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${SPATIALINDEX_INCLUDE_DIRS}
    ${TIFF_INCLUDE_DIRS}
)

# Add main library
add_library(PackingGenerator SHARED
    PackingGenerator.cpp
    PackingGenerator.h
)

# Link dependencies to main library
target_link_libraries(PackingGenerator
    ${SPATIALINDEX_LIBRARIES}
    ${TIFF_LIBRARIES}
)

# Set output name
set_target_properties(PackingGenerator PROPERTIES
    OUTPUT_NAME "packing_generator"
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

# Set up examples if enabled
if(BUILD_EXAMPLES)
    add_executable(example example.cpp)
    target_link_libraries(example PackingGenerator)
endif()

# Set up Python bindings if enabled
if(BUILD_PYTHON_BINDINGS)
    # Find Python
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    
    # Find pybind11
    find_package(pybind11 CONFIG REQUIRED)
    
    # Add Python module
    pybind11_add_module(particle_packing binding.cpp)
    
    # Link with our library
    target_link_libraries(particle_packing PRIVATE
        PackingGenerator
        ${SPATIALINDEX_LIBRARIES}
        ${TIFF_LIBRARIES}
    )
    
    # Set output directory
    set_target_properties(particle_packing PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
    )
    
    # Install Python module
    install(TARGETS particle_packing
        LIBRARY DESTINATION ${Python3_SITELIB}
    )
endif()

# Install library and headers
install(TARGETS PackingGenerator
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES PackingGenerator.h
    DESTINATION include
)

# Add uninstall target
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY
    )
    
    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
    )
endif()

# Output configuration information
message(STATUS "Configuration Summary:")
message(STATUS "  CMake version: ${CMAKE_VERSION}")
message(STATUS "  C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  Build Python bindings: ${BUILD_PYTHON_BINDINGS}")
message(STATUS "  Build examples: ${BUILD_EXAMPLES}")
message(STATUS "  Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  SpatialIndex include dir: ${SPATIALINDEX_INCLUDE_DIRS}")
message(STATUS "  TIFF include dir: ${TIFF_INCLUDE_DIRS}")