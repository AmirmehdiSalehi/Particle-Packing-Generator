/**
 * @file SpatialIndex.h
 * @brief Spatial indexing helper classes for the Random Packing Generator
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file contains helper classes for interfacing with the libspatialindex
 * library. These classes provide data wrappers and visitors for efficient
 * spatial queries of spheres and interface segments.
 */

#ifndef PACKING_SPATIAL_INDEX_H
#define PACKING_SPATIAL_INDEX_H

#include "Sphere.h"
#include "Interface.h"
#include <vector>
#include <cstdint>
#include <spatialindex/SpatialIndex.h>

namespace Packing {

/**
 * @class SegmentData
 * @brief Data wrapper for storing Interface in spatial index
 * 
 * This class implements the SpatialIndex::IData interface to allow
 * interface segments to be stored and queried in an R-tree.
 */
class SegmentData : public SpatialIndex::IData {
public:
    /**
     * @brief Constructs a segment data wrapper
     * @param id Unique identifier for the segment
     * @param segment Pointer to the Interface
     */
    SegmentData(uint32_t id, const Interface* segment);
    
    /**
     * @brief Destructor
     */
    ~SegmentData() override = default;
    
    /**
     * @brief Gets the identifier for this data item
     * @return Segment ID
     */
    SpatialIndex::id_type getIdentifier() const override;
    
    /**
     * @brief Gets the data payload
     * @param length Output parameter for data length
     * @param data Output parameter for data pointer
     */
    void getData(uint32_t& length, uint8_t** data) const override;
    
    /**
     * @brief Creates a shape representing the segment's bounding box
     * @return Dynamically allocated shape object
     */
    SpatialIndex::IShape* getShape() const override;
    
    /**
     * @brief Creates a shape via output parameter
     * @param shape Output parameter for the shape
     */
    void getShape(SpatialIndex::IShape** shape) const override;

private:
    uint32_t m_id;                    ///< Segment identifier
    const Interface* m_segment; ///< Pointer to the segment
};

/**
 * @class SegmentVisitor
 * @brief Visitor for collecting interface segments from spatial queries
 * 
 * This visitor collects all interface segments found during a spatial
 * query operation.
 */
class SegmentVisitor : public SpatialIndex::IVisitor {
public:
    /**
     * @brief Default constructor
     */
    SegmentVisitor() = default;
    
    /**
     * @brief Destructor
     */
    ~SegmentVisitor() override = default;
    
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
    
    /// Collection of found segments
    std::vector<const Interface*> foundSegments;
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
    
    /**
     * @brief Creates a shape representing the sphere
     * @return Dynamically allocated shape object
     */
    SpatialIndex::IShape* getShape() const override;
    
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

#endif // PACKING_SPATIAL_INDEX_H