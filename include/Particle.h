/**
 * @file Particle.h
 * @brief Particle class definition for the Random Packing Generator
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file defines the Particle class, which represents non-spherical
 * particles composed of multiple spheres arranged hierarchically.
 */

#ifndef PACKING_PARTICLE_H
#define PACKING_PARTICLE_H

#include "Core.h"
#include "Sphere.h"
#include <vector>
#include <set>
#include <cstdint>

namespace Packing {

/**
 * @class Particle
 * @brief Represents a non-spherical particle composed of multiple spheres
 * 
 * Particles are created by aggregating spheres in three stages:
 * 1. Core sphere(s) defining the center
 * 2. Secondary spheres shaping the overall form
 * 3. Tertiary spheres adding surface roughness and inter-particle contacts
 * 
 * The class tracks particle properties including volume, surface area,
 * sphericity, and coordination number (number of contacts).
 */
class Particle {
public:
    /**
     * @brief Constructs a particle with a unique identifier
     * @param id Unique identifier for this particle (1-based)
     */
    explicit Particle(uint16_t id);
    
    /**
     * @brief Gets the unique identifier of this particle
     * @return Particle ID
     */
    uint16_t getId() const { return id; }
    
    /**
     * @brief Adds a sphere to this particle
     * @param center Center coordinates of the sphere
     * @param radius Radius of the sphere
     * @param type Hierarchy level (core, secondary, tertiary)
     * @param sphereID Unique ID for spatial indexing
     * @return Const reference to the newly added sphere
     */
    const std::shared_ptr<Sphere>& addSphere(const Point3D& center, int radius, 
                           SphereType type, uint32_t sphereID);

    /**
     * @brief Gets all spheres composing this particle
     * @return Const reference to the sphere vector
     */
    const std::vector<std::shared_ptr<Sphere>>& getSpheres() const { return spheres; }
    
    /**
     * @brief Finds and returns the core sphere of this particle
     * @return Pointer to the core sphere, or nullptr if none found
     */
    const std::shared_ptr<Sphere> getCoreSphere() const;
    
    /**
     * @brief Gets the volume of this particle in voxels
     * @return Number of bulk voxels
     */
    uint32_t getVolume() const { return bulkCount; }
    
    /**
     * @brief Gets the surface area of this particle in voxels
     * @return Number of surface voxels
     */
    uint32_t getArea() const { return surfaceCount; }

    /**
     * @brief Increments the bulk voxel count
     */
    void incrementVolume() { bulkCount++; }

    /**
     * @brief Decrements the bulk voxel count
     */
    void decrementVolume() { bulkCount--; }

    /**
     * @brief Increments the surface voxel count
     */
    void incrementArea() { surfaceCount++; }
 
    /**
     * @brief Decrements the surface voxel count
     */
    void decrementArea() { surfaceCount--; }
    
    /**
     * @brief Calculates the sphericity index of this particle
     * @return Sphericity value (0-1, where 1 is a perfect sphere)
     * 
     * Uses the formula: Sphericity = (36π × V²) / S³
     * where V is volume and S is surface area
     */
    double calculateSphericity() const;
    
    /**
     * @brief Gets the set of particles in contact with this one
     * @return Set of contacting particle IDs
     */
    const std::set<uint16_t>& getContacts() const { return contacts; }
    
    /**
     * @brief Records a contact with another particle
     * @param particleId ID of the contacted particle
     */
    void addContact(uint16_t particleId);
    
    /**
     * @brief Gets the coordination number (number of contacts)
     * @return Number of particles in contact with this one
     */
    uint32_t getCoordinationNumber() const { return contacts.size(); }

    // Contact set made public to match original implementation
    std::set<uint16_t> contacts; ///< IDs of contacting particles
 
private:
    uint16_t id;                ///< Unique particle identifier
    uint32_t bulkCount;         ///< Number of interior voxels (volume)
    uint32_t surfaceCount;      ///< Number of surface voxels (area)
    std::vector<std::shared_ptr<Sphere>> spheres; ///< Component spheres of this particle
};

} // namespace Packing

#endif // PACKING_PARTICLE_H