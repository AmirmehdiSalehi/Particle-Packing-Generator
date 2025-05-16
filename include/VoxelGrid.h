/**
 * @file VoxelGrid.h
 * @brief VoxelGrid class definition for the Random Packing Generator
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file defines the VoxelGrid class, which provides a 3D voxel
 * representation of the particle packing with efficient chunk-based
 * memory management and spatial indexing for interface segments.
 */

#ifndef PACKING_VOXEL_GRID_H
#define PACKING_VOXEL_GRID_H

#include "Core.h"
#include "Particle.h"
#include "Interface.h"
#include <memory>
#include <unordered_map>
#include <string>
#include <cstdint>

namespace SpatialIndex {
    namespace RTree { class RTree; }
}

namespace Packing {

/**
 * @class VoxelGrid
 * @brief 3D voxel grid representation of the particle packing
 * 
 * The VoxelGrid class provides an efficient representation of the 3D
 * particle packing using a chunk-based storage system. It tracks
 * particle volumes, surface areas, and inter-particle contacts.
 * The grid uses 16-bit voxel values encoding both particle ID
 * and voxel type (surface/bulk).
 */
class VoxelGrid {
public:
    /**
     * @brief Constructs a voxel grid of specified size
     * @param size Cubic grid dimension (size × size × size)
     */
    explicit VoxelGrid(uint32_t size);
    
    /**
     * @brief Destructor - cleans up allocated memory
     */
    ~VoxelGrid();
    
    /**
     * @brief Gets the size of the grid
     * @return Grid dimension
     */
    uint32_t getSize() const { return size; }
           
    /**
     * @brief Adds a sphere to the voxel grid
     * @param particles Vector of all particles for tracking volumes
     * @param center Center coordinates of the sphere
     * @param radius Radius of the sphere
     * @param particleID ID of the particle this sphere belongs to
     * 
     * This method uses pre-computed sphere masks for efficiency and
     * automatically detects surface/bulk voxels and particle contacts.
     */
    void addSphere(std::vector<Particle>& particles, 
                   const Point3D& center, int radius, uint16_t particleID);

    /**
     * @brief Sets a voxel value and updates particle statistics
     * @param point Voxel coordinates
     * @param particles Vector of particles for updating counts
     * @param currentValue Reference to the current voxel value
     * @param newVoxelID New particle ID to set
     * @param newVoxelType New voxel type to set
     * @return true if value was changed, false otherwise
     */
    bool setVoxel(const Point3D& point, std::vector<Particle>& particles, 
                  uint16_t& currentValue, uint16_t newVoxelID, VoxelType newVoxelType);

    /**
     * @brief Records an interface segment at a contact point
     * @param point The voxel where contact occurs
     * 
     * This method manages interface segments, merging adjacent ones
     * and maintaining the spatial index for efficient queries.
     */
    void addInterfaceSegment(const Point3D& point);

    /**
     * @brief Gets the number of interface segments
     * @return Count of distinct contact regions
     */
    int getInterfaceSegmentCount() const { return interfacialSegments.size(); }

    /**
     * @brief Gets the value of a voxel
     * @param point Voxel coordinates
     * @return 16-bit voxel value encoding ID and type
     */
    uint16_t getVoxel(const Point3D& point) const;

    /**
     * @brief Extracts the voxel type from a voxel value
     * @param value The 16-bit voxel value
     * @return VoxelType (SURFACE or BULK)
     */
    VoxelType getVoxelType(const uint16_t& value) const;

    /**
     * @brief Sets the type portion of a voxel value
     * @param currentValue Reference to the voxel value
     * @param newVoxelType New type to set
     * @return true if successful
     */
    bool setVoxelType(uint16_t& currentValue, VoxelType newVoxelType);

    /**
     * @brief Extracts the particle ID from a voxel value
     * @param value The 16-bit voxel value
     * @return Particle ID (lower 15 bits)
     */
    uint16_t getVoxelID(const uint16_t& value) const;

    /**
     * @brief Sets the ID portion of a voxel value
     * @param currentValue Reference to the voxel value
     * @param newIndex New particle ID to set
     * @return true if successful
     */
    bool setVoxelID(uint16_t& currentValue, uint16_t newIndex);

    /**
     * @brief Gets the total number of filled voxels
     * @return Count of non-empty voxels
     */
    uint32_t getFilledVoxelCount() const { return filledVoxelCount; }

    /**
     * @brief Creates or retrieves a cached sphere mask
     * @param radius Radius of the sphere
     * @return Iterator to the mask and boolean for creation status
     * 
     * Sphere masks are pre-computed templates that speed up the
     * sphere addition process by avoiding redundant calculations.
     */
    std::pair<std::unordered_map<int, std::vector<VoxelType>>::iterator, bool> 
        createMask(int radius);

    /**
     * @brief Checks if coordinates are within grid bounds
     * @param point Point to check
     * @return true if within bounds
     */
    bool isInBounds(const Point3D& point) const;

    /**
     * @brief Saves the voxel grid as a multi-page TIFF file
     * @param filename Output filename
     * @param binary If true, save as binary (filled/empty), 
     *               otherwise save particle IDs
     * @return true if successful
     */
    bool saveToTIFF(const std::string& filename, bool binary = true) const;
 
private:
    // Grid dimensions and chunk management
    uint32_t size;                ///< Grid size in each dimension
    uint32_t chunkSize;           ///< Size of each memory chunk
    uint32_t numChunks;           ///< Number of chunks per dimension
    uint16_t** grid;              ///< 2D array of chunk pointers
    
    // Statistics and counters
    uint32_t filledVoxelCount;    ///< Number of non-empty voxels
    uint32_t interfaceSegsCount;  ///< Running counter for segment IDs

    // Optimization caches and spatial indices
    std::unordered_map<int, std::vector<VoxelType>> sphericalMasks; ///< Cached sphere templates
    std::unordered_map<int, Interface> interfacialSegments;  ///< Contact regions
    std::unique_ptr<SpatialIndex::ISpatialIndex> spatialIndexSegs;   ///< R-tree for segments

    /**
     * @brief Calculates chunk coordinates for a point
     * @param center Point coordinates
     * @return Chunk coordinates
     */
    Point3D getChunkCoords(const Point3D& center) const;

    /**
     * @brief Calculates local coordinates within a chunk
     * @param center Point coordinates
     * @return Local coordinates within the chunk
     */
    Point3D getLocalCoords(const Point3D& center) const;

    /**
     * @brief Gets the linear index of a chunk
     * @param center Point coordinates
     * @return Chunk index in the grid array
     */
    uint32_t getChunkIndex(const Point3D& center) const;

    /**
     * @brief Gets the local index within a chunk
     * @param center Point coordinates
     * @return Index within the chunk array
     */
    uint32_t getLocalIndex(const Point3D& center) const;
};

} // namespace Packing

#endif // PACKING_VOXEL_GRID_H