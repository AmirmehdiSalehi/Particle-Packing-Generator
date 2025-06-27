/**
 * @file VoxelGrid.cpp
 * @brief Implementation of the VoxelGrid class (Part 1)
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file contains the implementation of the VoxelGrid class methods
 * for managing the 3D voxel representation of particle packings.
 */

#include "VoxelGrid.h"
#include <tiffio.h>
#include <algorithm>
#include <iostream>
#include <list>

namespace Packing {

//=============================================================================
// VoxelGrid Implementation - Construction and Basic Operations
//=============================================================================

VoxelGrid::VoxelGrid(uint32_t size) 
    : size(size), 
      chunkSize(512), 
      filledVoxelCount(0), 
      interfaceSegsCount(0) {
    
    // Calculate number of chunks needed in each dimension
    numChunks = (size + chunkSize - 1) / chunkSize;
    
    // Allocate grid as array of chunk pointers
    // Using chunked storage for better memory efficiency
    grid = new uint16_t*[numChunks * numChunks * numChunks];
    
    // Initialize each chunk
    for (uint32_t i = 0; i < numChunks * numChunks * numChunks; ++i) {
        grid[i] = new uint16_t[chunkSize * chunkSize * chunkSize]();
    }

    // Initialize spatial index for interface segments
    SpatialIndex::IStorageManager* memoryManager = 
        SpatialIndex::StorageManager::createNewMemoryStorageManager();
    
    // R-tree configuration
    SpatialIndex::id_type indexIdentifier = 1; 
    SpatialIndex::RTree::RTreeVariant variant = SpatialIndex::RTree::RV_RSTAR;
    double fillFactor = 0.7;
    uint32_t indexCapacity = 10;
    uint32_t leafCapacity = 10;
    uint32_t dimension = 3;
    
    // Create R-tree for efficient segment queries
    spatialIndexSegs.reset(
        SpatialIndex::RTree::createNewRTree(
            *memoryManager, 
            fillFactor, 
            indexCapacity, 
            leafCapacity, 
            dimension, 
            variant, 
            indexIdentifier
        )
    );
}

VoxelGrid::~VoxelGrid() {
    // Clean up each chunk
    for (uint32_t i = 0; i < numChunks * numChunks * numChunks; ++i) {
        delete[] grid[i];
    }
    
    // Clean up the grid array
    delete[] grid;
}

//=============================================================================
// Sphere Addition Operations
//=============================================================================

void VoxelGrid::addSphere(std::vector<Particle>& particles, const Point3D& center, 
                          int radius, uint16_t particleID) {
    // Look up or create the sphere mask for this radius
    auto it = sphericalMasks.find(radius);
    if (it == sphericalMasks.end()) {
        auto result = createMask(radius);
        it = result.first;
    }
    
    // Get the pre-computed spherical mask
    const auto& mask = it->second;
    size_t maskSize = 2 * (radius + 1) + 1;

    // Calculate bounding box for the sphere
    Point3D minCorner(center.x - (radius + 1), 
                     center.y - (radius + 1), 
                     center.z - (radius + 1)); 
    Point3D maxCorner(center.x + (radius + 1), 
                     center.y + (radius + 1), 
                     center.z + (radius + 1));
    
    // Clamp to grid boundaries
    minCorner.x = std::max(minCorner.x, 0);
    minCorner.y = std::max(minCorner.y, 0);
    minCorner.z = std::max(minCorner.z, 0);
    
    maxCorner.x = std::min(maxCorner.x, static_cast<int>(size - 1));
    maxCorner.y = std::min(maxCorner.y, static_cast<int>(size - 1));
    maxCorner.z = std::min(maxCorner.z, static_cast<int>(size - 1));
    
    // Determine affected chunks
    const Point3D minChunk = getChunkCoords(minCorner);
    const Point3D maxChunk = getChunkCoords(maxCorner);

    // Process each affected chunk
    for (int chunkX = minChunk.x; chunkX <= maxChunk.x; chunkX++) {
        for (int chunkY = minChunk.y; chunkY <= maxChunk.y; chunkY++) {
            for (int chunkZ = minChunk.z; chunkZ <= maxChunk.z; chunkZ++) {
                // Calculate intersection with this chunk
                int startX = std::max(minCorner.x, static_cast<int>(chunkX * chunkSize));
                int startY = std::max(minCorner.y, static_cast<int>(chunkY * chunkSize));
                int startZ = std::max(minCorner.z, static_cast<int>(chunkZ * chunkSize));
                
                int endX = std::min(maxCorner.x, static_cast<int>((chunkX + 1) * chunkSize - 1));
                int endY = std::min(maxCorner.y, static_cast<int>((chunkY + 1) * chunkSize - 1));
                int endZ = std::min(maxCorner.z, static_cast<int>((chunkZ + 1) * chunkSize - 1));

                uint32_t chunkIndex = chunkZ + chunkY * numChunks + 
                                     chunkX * numChunks * numChunks;

                // Apply mask to voxels in this chunk
                for (int x = startX; x <= endX; x++) {
                    for (int y = startY; y <= endY; y++) {
                        for (int z = startZ; z <= endZ; z++) { 
                            Point3D voxelPoint(x, y, z);
                            uint32_t localIndex = getLocalIndex(voxelPoint);

                            // Calculate position in the mask
                            int maskX = x - (center.x - (radius + 1));
                            int maskY = y - (center.y - (radius + 1));
                            int maskZ = z - (center.z - (radius + 1));
                            
                            // Check mask bounds
                            if (maskX >= 0 && maskX < static_cast<int>(maskSize) && 
                                maskY >= 0 && maskY < static_cast<int>(maskSize) && 
                                maskZ >= 0 && maskZ < static_cast<int>(maskSize)) {
                                
                                uint32_t maskIndex = maskZ + maskY * maskSize + 
                                                   maskX * maskSize * maskSize; 
                                
                                // Apply mask value to voxel
                                setVoxel(voxelPoint, particles, 
                                        grid[chunkIndex][localIndex], 
                                        particleID, mask[maskIndex]);
                            }
                        }
                    }
                }
            }
        }
    }
}

//=============================================================================
// Voxel Manipulation Operations
//=============================================================================

bool VoxelGrid::setVoxel(const Point3D& point, std::vector<Particle>& particles, 
                         uint16_t& currentValue, uint16_t newVoxelID, 
                         VoxelType newVoxelType) {
    // Empty voxels don't change
    if (newVoxelType == VoxelType::AIR) {
        return false;
    }

    // Extract current voxel properties
    VoxelType currentVoxelType = getVoxelType(currentValue);
    uint16_t currentVoxelID = getVoxelID(currentValue);

    // Find particles involved
    Particle* newParticle = nullptr;
    Particle* currentParticle = nullptr;
    
    for (auto& particle : particles) {
        if (particle.getId() == newVoxelID) {
            newParticle = &particle;
        }
        if (particle.getId() == currentVoxelID) {
            currentParticle = &particle;
        }
        if (newParticle && currentParticle) {
            break;
        }
    }

    // Handle empty voxel
    if (currentValue == 0) {
        setVoxelID(currentValue, newVoxelID);
        setVoxelType(currentValue, newVoxelType);
        
        // Update particle statistics
        if (newParticle) {
            if (newVoxelType == VoxelType::SURFACE) { 
                newParticle->incrementArea(); 
            }
            if (newVoxelType == VoxelType::BULK) { 
                newParticle->incrementVolume(); 
            }
        }
        
        filledVoxelCount++;
        return true;
    } 
    
    // Handle already filled voxel
    if (currentVoxelID == newVoxelID) {
        // Same particle - possibly convert surface to bulk
        if (currentVoxelType == VoxelType::SURFACE && 
            newVoxelType == VoxelType::BULK) {
            setVoxelType(currentValue, VoxelType::BULK);
            
            if (currentParticle) {
                currentParticle->decrementArea();
                currentParticle->incrementVolume();
            }
            return true;
        }
        return false;
    } 
    else {
        // Different particles - potential contact
        if (currentVoxelType == VoxelType::SURFACE) {
            // Record contact between particles
            if (currentParticle && newParticle) {
                currentParticle->addContact(newVoxelID);
                newParticle->addContact(currentVoxelID);
            }

            // Create interface segment
            addInterfaceSegment(point);
            return true;
        }
        return false;
    }
}

//=============================================================================
// Interface Segment Management
//=============================================================================

void VoxelGrid::addInterfaceSegment(const Point3D& point) {
    // Create query point for spatial search
    double coords[3] = {
        static_cast<double>(point.x), 
        static_cast<double>(point.y), 
        static_cast<double>(point.z)
    };
    SpatialIndex::Point queryPoint(coords, 3);
    
    // Find segments that might contain this point
    IdVisitor visitor;
    spatialIndexSegs->intersectsWithQuery(queryPoint, visitor);
    
    std::vector<const uint64_t> containingSegments;
    
    // Check which segments actually contain the point
    for (const auto foundSegmentId : visitor.GetResults()) {
        auto it = interfacialSegments.find(foundSegmentId);
        if (it != interfacialSegments.end() && (it->second)->contains(point)) {
            containingSegments.push_back(foundSegmentId);
        }
    }
    
    // Case 1: No existing segment - create new
    if (containingSegments.empty()) {
        int segmentId = interfaceSegsCount++;
        auto result = interfacialSegments.emplace(segmentId, std::make_shared<Interface>(segmentId));
        std::shared_ptr<Interface>& segment = result.first->second;
        
        segment->addPoint(point);
        
        // Add to spatial index
        const auto& bbox = segment->getBoundingBox();
        double pLow[3] = {
            static_cast<double>(bbox[0]), 
            static_cast<double>(bbox[1]), 
            static_cast<double>(bbox[2])
        };
        double pHigh[3] = {
            static_cast<double>(bbox[3]), 
            static_cast<double>(bbox[4]), 
            static_cast<double>(bbox[5])
        };

        SpatialIndex::Region queryRegion(pLow, pHigh, 3);
        spatialIndexSegs->insertData(0, nullptr, queryRegion, segmentId);

    }
    // Case 2: One existing segment - add point to it
    else if (containingSegments.size() == 1) {
        const uint64_t& segmentId = containingSegments[0];
        
        auto it = interfacialSegments.find(segmentId);
        if (it != interfacialSegments.end()) {
            std::shared_ptr<Interface>& segment = it->second;
            
            // Only update if not already within bounding box
            if (!segment->withinBbox(point)) {
                segment->addPoint(point);
                
                // Update spatial index
                const auto& bbox = segment->getBoundingBox();
                double pLow[3] = {
                    static_cast<double>(bbox[0]), 
                    static_cast<double>(bbox[1]), 
                    static_cast<double>(bbox[2])
                };
                double pHigh[3] = {
                    static_cast<double>(bbox[3]), 
                    static_cast<double>(bbox[4]), 
                    static_cast<double>(bbox[5])
                };
                SpatialIndex::Region queryRegion(pLow, pHigh, 3);
                
                spatialIndexSegs->deleteData(queryRegion, segmentId);
                spatialIndexSegs->insertData(0, nullptr, queryRegion, segmentId);
            }
        }
    }
    // Case 3: Multiple segments - merge them
    else if (containingSegments.size() > 1) {
        std::vector<std::shared_ptr<Interface>> segmentsToMerge;
        for (const auto segmentId : containingSegments) {
            auto it = interfacialSegments.find(segmentId);
            if (it != interfacialSegments.end()) {
                segmentsToMerge.push_back(it->second);
            }
        }
        
        if (!segmentsToMerge.empty()) {
            // Use first segment as merge target
            std::shared_ptr<Interface>& targetSegment = segmentsToMerge[0];
            uint64_t targetId = targetSegment->getId();
            
            targetSegment->addPoint(point);
            targetSegment->mergeSegments(segmentsToMerge, interfacialSegments);
            
            // Update spatial index
            const auto& bbox = targetSegment->getBoundingBox();
            double pLow[3] = {
                static_cast<double>(bbox[0]), 
                static_cast<double>(bbox[1]), 
                static_cast<double>(bbox[2])
            };
            double pHigh[3] = {
                static_cast<double>(bbox[3]), 
                static_cast<double>(bbox[4]), 
                static_cast<double>(bbox[5])
            };
            SpatialIndex::Region queryRegion(pLow, pHigh, 3);
            
            // Remove old entries
            for (const auto segmentId : containingSegments) {
                spatialIndexSegs->deleteData(queryRegion, segmentId);
            }
            
            // Insert merged segment
            spatialIndexSegs->insertData(0, nullptr, queryRegion, targetId);

        }
    }
}

//=============================================================================
// Voxel Property Operations
//=============================================================================

uint16_t VoxelGrid::getVoxel(const Point3D& point) const {
    if (!isInBounds(point)) {
        return 0;
    }

    uint32_t chunkIndex = getChunkIndex(point);
    uint32_t localIndex = getLocalIndex(point);
    
    return grid[chunkIndex][localIndex];
}

VoxelType VoxelGrid::getVoxelType(const uint16_t& value) const {
    if (value == 0) {
        return VoxelType::AIR;
    }
    
    // MSB indicates BULK (1) or SURFACE (0)
    bool isBulk = (value >> 15) & 0x1;
    
    return isBulk ? VoxelType::BULK : VoxelType::SURFACE;
}

bool VoxelGrid::setVoxelType(uint16_t& currentValue, VoxelType newVoxelType) {
    if (newVoxelType == VoxelType::BULK) {
        currentValue |= (1 << 15);  // Set MSB for BULK
        return true;
    } else if (newVoxelType == VoxelType::SURFACE) {
        currentValue &= ~(1 << 15);  // Clear MSB for SURFACE
        return true;
    }
    
    return false;
}

uint16_t VoxelGrid::getVoxelID(const uint16_t& value) const {
    // Lower 15 bits contain the particle ID
    return value & 0x7FFF;
}

bool VoxelGrid::setVoxelID(uint16_t& currentValue, uint16_t newIndex) {
    // Check if ID fits in 15 bits
    if (newIndex > 0x7FFF) {
        return false;
    }
    
    // Preserve the type bit (MSB) and set new ID
    uint16_t typeBit = currentValue & 0x8000;
    currentValue = (typeBit | (newIndex & 0x7FFF));
    return true;
}

//=============================================================================
// Chunk Management Operations
//=============================================================================

Point3D VoxelGrid::getChunkCoords(const Point3D& center) const {
    return Point3D(center.x / chunkSize,
                   center.y / chunkSize,
                   center.z / chunkSize);
}

Point3D VoxelGrid::getLocalCoords(const Point3D& center) const {
    return Point3D(center.x % chunkSize,
                   center.y % chunkSize,
                   center.z % chunkSize);
}

uint32_t VoxelGrid::getChunkIndex(const Point3D& center) const {
    Point3D chunk_coords = getChunkCoords(center);
    
    return chunk_coords.z + 
           chunk_coords.y * numChunks + 
           chunk_coords.x * numChunks * numChunks;
}

uint32_t VoxelGrid::getLocalIndex(const Point3D& center) const {
    Point3D local_coords = getLocalCoords(center);
    
    return local_coords.z + 
           local_coords.y * chunkSize + 
           local_coords.x * chunkSize * chunkSize;
}

bool VoxelGrid::isInBounds(const Point3D& point) const {
    return (point.x >= 0 && point.x < static_cast<int>(size) &&
            point.y >= 0 && point.y < static_cast<int>(size) &&
            point.z >= 0 && point.z < static_cast<int>(size));
}

//=============================================================================
// Sphere Mask Creation
//=============================================================================

std::pair<std::unordered_map<int, std::vector<VoxelType>>::iterator, bool> 
VoxelGrid::createMask(int radius) {
    // List to track potential surface voxels
    std::list<Point3D> potentialSurfaceVoxels;
    
    // Create mask grid
    size_t maskSize = 2 * (radius + 1) + 1;
    std::vector<VoxelType> maskGrid(maskSize * maskSize * maskSize, VoxelType::AIR);
    
    // Fill sphere with bulk voxels
    for (int dx = -(radius + 1); dx <= (radius + 1); ++dx) {
        for (int dy = -(radius + 1); dy <= (radius + 1); ++dy) {
            for (int dz = -(radius + 1); dz <= (radius + 1); ++dz) {
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                // Convert to positive indices
                int x = dx + (radius + 1);
                int y = dy + (radius + 1);
                int z = dz + (radius + 1);
                
                if (dist <= radius) {
                    maskGrid[z + y * maskSize + x * maskSize * maskSize] = VoxelType::BULK;

                    // Mark near-surface voxels for later processing
                    if (radius - dist <= 0.8) {
                        potentialSurfaceVoxels.emplace_back(x, y, z);
                    }
                }
            }
        }
    }

    // Identify surface voxels
    const int dx[] = {-1, 1, 0, 0, 0, 0};
    const int dy[] = {0, 0, -1, 1, 0, 0};
    const int dz[] = {0, 0, 0, 0, -1, 1};
    
    for (const auto& voxel : potentialSurfaceVoxels) {
        for (int i = 0; i < 6; i++) {
            int nx = voxel.x + dx[i];
            int ny = voxel.y + dy[i];
            int nz = voxel.z + dz[i];
            
            if (nx >= 0 && nx < static_cast<int>(maskSize) && 
                ny >= 0 && ny < static_cast<int>(maskSize) && 
                nz >= 0 && nz < static_cast<int>(maskSize)) {
                
                size_t neighborIndex = nz + ny * maskSize + nx * maskSize * maskSize;
                
                // If neighbor is empty, mark current as surface
                if (maskGrid[neighborIndex] == VoxelType::AIR) {
                    size_t voxelIndex = voxel.z + voxel.y * maskSize + 
                                       voxel.x * maskSize * maskSize;
                    maskGrid[voxelIndex] = VoxelType::SURFACE;
                    break;
                }
            }
        }
    }

    // Cache the mask
    return sphericalMasks.emplace(radius, maskGrid);
}

//=============================================================================
// File Output Operations
//=============================================================================

bool VoxelGrid::saveToTIFF(const std::string& filename, bool binary) const {
    // Open TIFF file for writing
    TIFF* tif = TIFFOpen(filename.c_str(), "w");
    if (!tif) {
        std::cerr << "Failed to open TIFF file: " << filename << std::endl;
        return false;
    }
    
    // Buffer for one Z-slice
    uint16_t* buffer = new uint16_t[size * size];
    
    // Write each Z-slice as a TIFF page
    for (uint32_t z = 0; z < size; ++z) {
        // Configure TIFF tags for this page
        TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, size);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, size);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
        TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, size));
        TIFFSetField(tif, TIFFTAG_PAGENUMBER, z, size);
        
        // Fill buffer with slice data
        for (uint32_t y = 0; y < size; ++y) {
            for (uint32_t x = 0; x < size; ++x) {
                Point3D point(x, y, z);
                uint16_t val = getVoxel(point);
                
                // Convert to binary if requested
                if (binary) {
                    buffer[y * size + x] = val > 0 ? 65535 : 0;
                } else {
                    buffer[y * size + x] = val;
                }
            }
        }
        
        // Write scanlines
        for (uint32_t row = 0; row < size; ++row) {
            if (TIFFWriteScanline(tif, &buffer[row * size], row, 0) < 0) {
                std::cerr << "Failed to write scanline " << row << std::endl;
                delete[] buffer;
                TIFFClose(tif);
                return false;
            }
        }
        
        // Write directory for this page
        if (!TIFFWriteDirectory(tif)) {
            std::cerr << "Failed to write directory for slice " << z << std::endl;
            delete[] buffer;
            TIFFClose(tif);
            return false;
        }
    }
    
    // Cleanup
    delete[] buffer;
    TIFFClose(tif);
    
    return true;
}

} // namespace Packing