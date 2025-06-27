/**
 * @file PackingGenerator.h
 * @brief Main PackingGenerator class definition
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file defines the PackingGenerator class, which implements the
 * Random Sequential Addition (RSA) algorithm for generating packings
 * of non-spherical particles with tunable topology.
 */

#ifndef PACKING_GENERATOR_H
#define PACKING_GENERATOR_H

#include "Core.h"
#include "Particle.h"
#include "VoxelGrid.h"
#include <memory>
#include <vector>
#include <string>
#include <cstdint>
#include <utility>

// Forward declaration for spatial indexing
namespace SpatialIndex {
    namespace RTree { class RTree; }
}

namespace Packing {

/**
 * @class PackingGenerator
 * @brief Main class for generating random packings of non-spherical particles
 * 
 * This class implements a three-stage Random Sequential Addition (RSA) algorithm
 * to create packings of complex particles. The algorithm proceeds as follows:
 * 1. Core sphere insertion to establish particle centers
 * 2. Secondary sphere addition to shape particle forms
 * 3. Tertiary sphere placement to create surface roughness and contacts
 * 
 * The generator allows control over particle shape, size distribution,
 * packing density, and inter-particle connectivity.
 */
class PackingGenerator {
public:
    /**
     * @brief Constructs a packing generator with specified parameters
     * @param size Cubic domain dimension (size × size × size voxels)
     * @param coreRadiusMin Minimum radius for core spheres
     * @param coreRadiusMax Maximum radius for core spheres
     * @param secondaryRadiusMin Minimum radius for secondary spheres
     * @param secondaryRadiusMax Maximum radius for secondary spheres
     * @param tertiaryRadiusMin Minimum radius for tertiary spheres
     * @param tertiaryRadiusMax Maximum radius for tertiary spheres
     * @param tertiaryVolumeFraction Volume fraction of tertiary spheres (0-1)
     * @param targetDensity Target packing density to achieve
     * @param compactnessFactor Controls sphere overlap (0-1, default 0.5)
     */
    PackingGenerator(
        uint32_t size,
        int coreRadiusMin = 10, int coreRadiusMax = 20,
        int secondaryRadiusMin = 7, int secondaryRadiusMax = 12,
        int tertiaryRadiusMin = 2, int tertiaryRadiusMax = 7,
        double targetDensity = 0.6,
        double compactnessFactor = 0.5,
        uint32_t randomSeed = 0
    );
    
    /**
     * @brief Destructor
     */
    ~PackingGenerator();
    
    /**
     * @brief Generates the particle packing
     * @return true if generation was successful
     * 
     * This method executes the three-stage algorithm to create
     * the packing. Progress is reported to standard output.
     */
    bool generate();
    
    /**
     * @brief Gets the current packing density
     * @return Density value (0-1) as filled volume fraction
     */
    double getCurrentDensity() const;
    
    /**
     * @brief Gets a particle by index
     * @param index Particle index (0-based)
     * @return Pointer to particle or nullptr if index is invalid
     */
    const Particle* getParticle(uint32_t index) const;
    
    /**
     * @brief Gets a sphere by index
     * @param index Particle index (0-based)
     * @return Pointer to particle or nullptr if index is invalid
     */
    const std::shared_ptr<Sphere> getSphere(uint32_t index) const;
    
    /**
     * @brief Gets the total number of particles
     * @return Number of particles in the packing
     */
    uint32_t getParticleCount() const { return particles.size(); }
    
    /**
     * @brief Gets the total number of particles
     * @return Number of particles in the packing
     */
    uint32_t getSphereCount() const { return spheres.size(); }
    
    /**
     * @brief Gets all particles in the packing
     * @return Const reference to the particle vector
     */
    const std::vector<Particle>& getParticles() const { return particles; }
    
    /**
     * @brief Gets the total number of inter-particle contacts
     * @return Number of distinct contact regions
     */
    uint32_t getContactCount() const;
    
    /**
     * @brief Calculates the average coordination number
     * @return Mean number of contacts per particle
     */
    double getAverageCoordinationNumber() const;
    
    /**
     * @brief Gets coordination numbers for all particles
     * @return Vector of coordination numbers
     */
    std::vector<uint32_t> getCoordinationNumbers() const;
    
    /**
     * @brief Calculates the average sphericity index
     * @return Mean sphericity value (0-1)
     */
    double getAverageSphericity() const;
    
    /**
     * @brief Calculates the average particle radius 
     * @return Mean particle radius using the equivalent sphere volume
     */
    double getAverageParticleRadius() const;
    
    /**
     * @brief Calculates the total number of filled voxels
     * @return The total volume occupied in the voxel grid by the particles
     */
    uint32_t getTotalVolume() const;

    /**
     * @brief Saves the packing as a 3D TIFF image
     * @param filename Output filename
     * @param binary If true, save as binary (filled/empty),
     *               otherwise save particle IDs
     * @return true if save was successful
     */
    bool saveTIFF(const std::string& filename, bool binary = true) const;

    /**
     * @brief Gets all contact pairs between particles 
     * @return Vector of pairs (particleID1, particleID2) representing contacts
     * 
     * Used by the RL environemnt to build edges among the core sphere nodes in the graph representation
     * Each contact is represented only once (no duplicates).
     * Uses the contact information stored in each particle.
     */
    std::set<std::pair<uint16_t, uint16_t>> getContactPairs() const;

    std::vector<uint64_t> getSphereNeighbors(int x, int y, int z, int radius) const;

    /**
     * @brief Gets the volume of a particle given its index
     * @return particle volume
     * 
     * Used by the RL environemnt to build the core sphere nodes in the graph representation
     * Uses the number of 'bulk' voxels in each particle and takes particle overlaps into account
     */
    uint32_t getParticleVolume(uint16_t particleId) const;

    /**
     * @brief Inserts a supplementary sphere into a specific particle
     * @param center Center position of the sphere
     * @param radius Radius of the sphere
     * @param particleID ID of the particle to add the sphere to
     * @return true if sphere was successfully added, false otherwise
     * 
     * This method adds a sphere to an existing particle following the same
     * compactness criteria as secondary spheres, but without overlap checking
     * with other particles. If the sphere doesn't meet compactness criteria
     * with existing spheres in the target particle, the method returns false.
     */
    bool insertSuppSpheres(int x, int y, int z, int radius, uint16_t particleID);

    /**
     * @brief Inserts core spheres into the domain
     * @return Number of core spheres successfully placed
     * 
     * Core spheres are placed randomly without overlap to establish
     * particle centers. Placement continues until the target density
     * fraction is reached or no more spheres can be placed.
     * Made public for RL integration.
     */
    uint32_t insertCoreSpheres();
 
private:
    // Configuration parameters
    uint32_t size;                  ///< Domain size in voxels
    uint32_t numSpheres;            ///< Total spheres added
    int coreRadiusMin;              ///< Min core sphere radius
    int coreRadiusMax;              ///< Max core sphere radius
    int secondaryRadiusMin;         ///< Min secondary radius
    int secondaryRadiusMax;         ///< Max secondary radius
    int tertiaryRadiusMin;          ///< Min tertiary radius
    int tertiaryRadiusMax;          ///< Max tertiary radius
    double tertiaryVolumeFraction;  ///< Tertiary volume fraction
    double targetDensity;           ///< Target packing density
    double compactnessFactor;       ///< Sphere overlap control
    
    // Core data structures
    std::vector<Particle> particles;          ///< All particles
    std::vector<std::shared_ptr<Sphere>> spheres; 
    VoxelGrid voxelGrid;                      ///< 3D voxel representation
    std::unique_ptr<SpatialIndex::ISpatialIndex> spatialIndex; ///< Sphere R-tree
    

    
    /**
     * @brief Stage 2: Adds secondary spheres to existing particles
     * @param maxAttempts Maximum placement attempts before stopping
     * @return Number of secondary spheres successfully placed
     * 
     * Secondary spheres are added around core spheres to define the
     * overall particle shape. They must satisfy the compactness
     * criteria with existing spheres in the same particle.
     */
    uint32_t addSecondarySpheres(uint32_t maxAttempts = 5000);
    
    /**
     * @brief Stage 3: Adds tertiary spheres for roughness and contacts
     * @param maxAttempts Maximum placement attempts before stopping
     * @return Number of tertiary spheres successfully placed
     * 
     * Tertiary spheres are added to create surface roughness and
     * inter-particle contacts. They are smaller and placed at the
     * particle boundaries.
     */
    uint32_t addTertiarySpheres(uint32_t maxAttempts = 5000);
    
    /**
     * @brief Generates a random point within the domain
     * @return Random 3D point
     */
    Point3D getRandomPoint() const;
    
    /**
     * @brief Generates a random radius within specified bounds
     * @param minRadius Minimum radius value
     * @param maxRadius Maximum radius value
     * @return Random radius value
     */
    int getRandomRadius(int minRadius, int maxRadius) const;

};

} // namespace Packing

#endif // PACKING_GENERATOR_H