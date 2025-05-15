/**
 * @file CApi.h
 * @brief C API for the Random Packing Generator library
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file provides a C-compatible interface to the Random Packing
 * Generator library, allowing it to be used from C programs and
 * other languages that can interface with C libraries.
 */

#ifndef PACKING_C_API_H
#define PACKING_C_API_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
    #ifdef PACKING_EXPORTS
        #define PACKING_API __declspec(dllexport)
    #else
        #define PACKING_API __declspec(dllimport)
    #endif
#else
    #define PACKING_API
#endif

/**
 * @brief Creates a new packing generator instance
 * @param size Size of the cubic domain in voxels
 * @param coreRadiusMin Minimum radius for core spheres
 * @param coreRadiusMax Maximum radius for core spheres
 * @param secondaryRadiusMin Minimum radius for secondary spheres
 * @param secondaryRadiusMax Maximum radius for secondary spheres
 * @param tertiaryRadiusMin Minimum radius for tertiary spheres
 * @param tertiaryRadiusMax Maximum radius for tertiary spheres
 * @param tertiaryVolumeFraction Volume fraction of tertiary spheres (0-1)
 * @param targetDensity Target packing density
 * @param compactnessFactor Controls sphere overlap (0-1)
 * @return Opaque pointer to the generator instance
 */
PACKING_API void* CreatePackingGenerator(
    uint32_t size,
    int coreRadiusMin, int coreRadiusMax,
    int secondaryRadiusMin, int secondaryRadiusMax,
    int tertiaryRadiusMin, int tertiaryRadiusMax,
    double tertiaryVolumeFraction,
    double targetDensity,
    double compactnessFactor
);

/**
 * @brief Frees a packing generator instance
 * @param generator Pointer to the generator to free
 */
PACKING_API void FreePackingGenerator(void* generator);

/**
 * @brief Generates the particle packing
 * @param generator Pointer to the generator
 * @return 1 if successful, 0 otherwise
 */
PACKING_API int Generate(void* generator);

/**
 * @brief Gets the current packing density
 * @param generator Pointer to the generator
 * @return Density value (0-1)
 */
PACKING_API double GetDensity(void* generator);

/**
 * @brief Gets the total number of particles
 * @param generator Pointer to the generator
 * @return Number of particles
 */
PACKING_API uint32_t GetParticleCount(void* generator);

/**
 * @brief Gets the total number of inter-particle contacts
 * @param generator Pointer to the generator
 * @return Number of contacts
 */
PACKING_API uint32_t GetContactCount(void* generator);

/**
 * @brief Gets the average coordination number
 * @param generator Pointer to the generator
 * @return Average coordination number
 */
PACKING_API double GetAverageCoordinationNumber(void* generator);

/**
 * @brief Gets coordination numbers for all particles
 * @param generator Pointer to the generator
 * @param result Pre-allocated array to store results
 * @return Pointer to the result array
 * 
 * The caller must allocate an array of size GetParticleCount()
 * before calling this function.
 */
PACKING_API uint32_t* GetCoordinationNumbers(void* generator, uint32_t* result);

/**
 * @brief Gets the average sphericity of particles
 * @param generator Pointer to the generator
 * @return Average sphericity (0-1)
 */
PACKING_API double GetAverageSphericity(void* generator);

/**
 * @brief Saves the packing as a 3D TIFF file
 * @param generator Pointer to the generator
 * @param filename Output filename
 * @param binary 1 for binary output, 0 for particle IDs
 * @return 1 if successful, 0 otherwise
 */
PACKING_API int SaveTIFF(void* generator, const char* filename, int binary);

#ifdef __cplusplus
}
#endif

#endif // PACKING_C_API_H