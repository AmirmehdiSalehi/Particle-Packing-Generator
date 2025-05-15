/**
 * @file Interface.h
 * @brief Interface class definition for the Random Packing Generator
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file defines the Interface class, which represents regions
 * where two particles are in contact. These segments are crucial for
 * tracking inter-particle connections and calculating coordination numbers.
 */

#ifndef PACKING_INTERFACE_H
#define PACKING_INTERFACE_H

#include "Core.h"
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <vector>
#include <cstdint>

namespace Packing {

/**
 * @class Interface
 * @brief Represents a contact region between two particles
 * 
 * Interface segments are collections of voxels at the boundary between
 * particles. They are used to track particle contacts and calculate
 * coordination numbers. The class maintains a set of voxels and a
 * bounding box for efficient spatial queries.
 */
class Interface {
public:
    /**
     * @brief Constructs an interface segment with a unique identifier
     * @param id Unique identifier for this segment
     */
    explicit Interface(int id);
    
    /**
     * @brief Destructor
     */
    ~Interface();
    
    /**
     * @brief Adds a point to this interface segment
     * @param point The voxel coordinate to add
     * 
     * When adding a point, the segment also includes its 18-connected
     * neighbors (face and edge connections) to ensure proper coverage
     * of the interface region.
     */
    void addPoint(Point3D& point);

    /**
     * @brief Gets the axis-aligned bounding box of this segment
     * @return Array containing [minX, minY, minZ, maxX, maxY, maxZ]
     */
    const std::array<uint32_t, 6>& getBoundingBox() const { return bbox; }

    /**
     * @brief Checks if a point is contained within this segment
     * @param point The point to check
     * @return true if the point is in this segment, false otherwise
     */
    bool contains(Point3D& point) const;

    /**
     * @brief Checks if a point is within the bounding box
     * @param point The point to check
     * @return true if the point is inside the bounding box
     */
    bool withinBbox(Point3D& point) const;

    /**
     * @brief Merges another bounding box with this segment's bbox
     * @param otherBbox The bounding box to merge
     */
    void mergeBbox(std::array<uint32_t, 6>& otherBbox);

    /**
     * @brief Merges multiple segments into this one
     * @param segments The segments to absorb
     * @param interfacialSegments Global segment map for cleanup
     * @return true if merge was successful
     * 
     * This method combines multiple adjacent segments into a single
     * larger segment, updating the spatial index accordingly.
     */
    bool mergeSegments(std::vector<Interface*>& segments, 
                      std::unordered_map<int, Interface>& interfacialSegments);

    /**
     * @brief Gets the unique identifier of this segment
     * @return Segment ID
     */
    int getId() const { return id; }

private:
    int id;                                            ///< Unique segment identifier
    std::unordered_set<Point3D, Point3D::Hash> segment; ///< Set of voxels in this segment
    std::array<uint32_t, 6> bbox;                      ///< Bounding box [minX,minY,minZ,maxX,maxY,maxZ]
};

} // namespace Packing

#endif // PACKING_INTERFACE_H