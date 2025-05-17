/**
 * @file Sphere.h
 * @brief Sphere class definition for the Random Packing Generator
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file defines the Sphere class, which represents individual spherical
 * components that combine to form non-spherical particles in the packing.
 */

#ifndef PACKING_SPHERE_H
#define PACKING_SPHERE_H

#include "Core.h"
#include <vector>
#include <cstdint>
#include <spatialindex/SpatialIndex.h>

namespace Packing {

/**
 * @class Sphere
 * @brief Represents a single sphere in the particle system
 * 
 * Each sphere is a building block of a particle, with specific radius,
 * type, and association to a parent particle. Spheres are added in three
 * stages (core, secondary, tertiary) to build complex particle shapes.
 */
class Sphere {
public:
    /**
     * @brief Constructs a sphere with specified properties
     * @param center Center coordinates in voxel space
     * @param radius Sphere radius in voxels
     * @param type Hierarchy level (core, secondary, or tertiary)
     * @param particleId ID of the parent particle
     * @param sphereId Unique identifier for spatial indexing
     */
    Sphere(const Point3D& center, int radius, SphereType type, 
           uint16_t particleId, uint32_t sphereId);
    
    /**
     * @brief Gets the center coordinates of the sphere
     * @return Const reference to the center point
     */
    const Point3D& getCenter() const { return center; }
    
    /**
     * @brief Gets the sphere radius
     * @return Radius value in voxels
     */
    int getRadius() const { return radius; }
    
    /**
     * @brief Gets the hierarchy type of the sphere
     * @return SphereType enum value
     */
    SphereType getType() const { return type; }
    
    /**
     * @brief Gets the ID of the parent particle
     * @return Parent particle's unique identifier
     */
    uint16_t getParticleId() const { return particleId; }

    /**
     * @brief Gets the unique identifier of this sphere
     * @return Sphere's ID for spatial indexing
     */
    uint32_t getSphereId() const { return sphereId; }
    
    /**
     * @brief Checks if this sphere satisfies the compactness criteria with another
     * @param other The sphere to check against
     * @param beta Compactness factor (0-1) controlling allowed overlap
     * @return true if spheres meet the compactness criteria, false otherwise
     * 
     * This method implements the compactness criteria from the paper:
     * For two spheres to be properly packed, they must satisfy:
     * r2 - r1 < d ≤ r2 - (β × r1), where d is the center distance
     */
    bool intersectsWith(const Sphere& other, double beta = 0.5) const;
      
private:
    Point3D center;      ///< Center coordinates in voxel space
    int radius;          ///< Sphere radius in voxels
    SphereType type;     ///< Hierarchy level of this sphere
    uint16_t particleId; ///< ID of the parent particle
    uint32_t sphereId;   ///< Unique ID for spatial indexing
};



/**
 * @class SphereData
 * @brief Data wrapper for storing Sphere in spatial index
 * 
 * This class implements the SpatialIndex::IData interface to allow
 * spheres to be stored and queried in an R-tree.
 */
class SphereData : public SpatialIndex::IData {
public:
    /**
     * @brief Constructs a sphere data wrapper
     * @param id Unique identifier for the sphere
     * @param sphere Pointer to the Sphere
     */
    SphereData(uint32_t id, const Sphere* sphere);
    
    /**
     * @brief Destructor
     */
    ~SphereData() override = default;
    
    /**
     * @brief Gets the identifier for this data item
     * @return Sphere ID
     */
    SpatialIndex::id_type getIdentifier() const override;
    
    /**
     * @brief Gets the data payload
     * @param length Output parameter for data length
     * @param data Output parameter for data pointer
     */

    void getData(uint32_t& length, uint8_t** data) const override;

    // const std::shared_ptr<Sphere> getPointer() const;


    SphereData* clone() override;
    
    /**
     * @brief Creates a shape representing the sphere
     * @return Dynamically allocated shape object
     */
    // SpatialIndex::IShape* getShape() const;
    
    /**
     * @brief Creates a shape via output parameter
     * @param shape Output parameter for the shape
     */
    void getShape(SpatialIndex::IShape** shape) const override;

private:
    uint32_t m_id;        ///< Sphere identifier
    const Sphere* m_sphere; ///< Pointer to the sphere
};

/**
 * @class SphereVisitor
 * @brief Visitor for collecting spheres from spatial queries
 * 
 * This visitor collects all spheres found during a spatial query operation.
 */
class SphereVisitor : public SpatialIndex::IVisitor {
public:
    /**
     * @brief Default constructor
     */
    SphereVisitor() = default;
    
    /**
     * @brief Destructor
     */
    ~SphereVisitor() override = default;
    
    /**
     * @brief Called when visiting a node (not used)
     * @param node The node being visited
     */
    void visitNode(const SpatialIndex::INode& node) override;
    
    /**
     * @brief Called when visiting data
     * @param data The data item found
     */
    void visitData(const SpatialIndex::IData& data) override;
    
    /**
     * @brief Called when visiting multiple data items
     * @param v Vector of data items
     */
    void visitData(std::vector<const SpatialIndex::IData*>& v) override;
    
    /// Collection of found spheres
    std::vector<const Sphere*> foundSpheres;
};

} // namespace Packing

#endif // PACKING_SPHERE_H

