/**
 * @file Core.h
 * @brief Core data structures and enumerations for the Random Packing Generator
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file contains fundamental data structures and enumerations used throughout
 * the Random Packing Generator library, including 3D point representation and
 * type definitions for spheres and voxels.
 */

#ifndef PACKING_CORE_H
#define PACKING_CORE_H

#include <cstdint>
#include <cmath>
#include <functional>

namespace Packing {

/**
 * @enum SphereType
 * @brief Defines the hierarchy level of spheres in the particle clustering process
 * 
 * The algorithm constructs particles through three stages, each adding a different
 * type of sphere to build the final non-spherical particle shape.
 */
enum class SphereType {
    CORE,       ///< Core spheres forming the particle centers (Stage 1)
    SECONDARY,  ///< Secondary spheres defining the overall shape (Stage 2)
    TERTIARY    ///< Tertiary spheres for surface roughness and connectivity (Stage 3)
};

/**
 * @enum VoxelType
 * @brief Categorizes voxels in the 3D grid representation
 * 
 * Each voxel in the grid is classified based on whether it's inside a particle,
 * on the surface, or in the surrounding space.
 */
enum class VoxelType {
    SURFACE,  ///< Surface voxels marking particle boundaries
    BULK,     ///< Interior voxels within particles
    AIR,      ///< Empty space outside particles
    GAP       ///< Transition voxels between particles (used for contact detection)
};

/**
 * @struct Point3D
 * @brief Represents a discrete 3D coordinate in the voxel grid
 * 
 * This structure provides basic 3D point functionality including
 * distance calculations and hash operations for use in containers.
 */
struct Point3D {
    int x;  ///< X coordinate in voxel space
    int y;  ///< Y coordinate in voxel space
    int z;  ///< Z coordinate in voxel space
    
    /**
     * @brief Constructs a Point3D with specified coordinates
     * @param x X coordinate (default: 0)
     * @param y Y coordinate (default: 0)
     * @param z Z coordinate (default: 0)
     */
    Point3D(int x = 0, int y = 0, int z = 0) : x(x), y(y), z(z) {}
    
    /**
     * @brief Checks if two points have the same coordinates
     * @param other The point to compare with
     * @return true if coordinates are identical, false otherwise
     */
    bool operator==(const Point3D& other) const {
        return x == other.x && y == other.y && z == other.z;
    }

    /**
     * @brief Calculates the Euclidean distance to another point
     * @param other The target point
     * @return The distance as a double value
     */
    double distanceTo(const Point3D& other) const {
        double dx = static_cast<double>(x - other.x);
        double dy = static_cast<double>(y - other.y);
        double dz = static_cast<double>(z - other.z);
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    /**
     * @struct Hash
     * @brief Hash function for using Point3D in unordered containers
     * 
     * Combines the three coordinates into a single hash value
     * using bitwise operations for efficient map/set operations.
     */
    struct Hash {
        size_t operator()(const Point3D& v) const {
            return std::hash<int>()(v.x) ^ 
                  (std::hash<int>()(v.y) << 1) ^ 
                  (std::hash<int>()(v.z) << 2);
        }
    };
};

} // namespace Packing

#endif // PACKING_CORE_H