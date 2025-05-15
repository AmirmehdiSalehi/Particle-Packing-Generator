/**
 * @file PackingGenerator.cpp
 * @brief Implementation of the PackingGenerator class (Part 1)
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file contains the implementation of the main PackingGenerator
 * class that orchestrates the three-stage particle packing algorithm.
 */

#include "PackingGenerator.h"
#include "SpatialIndex.h"
#include <random>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <ctime>

namespace Packing {

//=============================================================================
// Constructor and Initialization
//=============================================================================

PackingGenerator::PackingGenerator(
    uint32_t size,
    int coreRadiusMin, int coreRadiusMax,
    int secondaryRadiusMin, int secondaryRadiusMax,
    int tertiaryRadiusMin, int tertiaryRadiusMax,
    double tertiaryVolumeFraction,
    double targetDensity,
    double compactnessFactor
) : size(size),
    numSpheres(0),
    coreRadiusMin(coreRadiusMin), coreRadiusMax(coreRadiusMax),
    secondaryRadiusMin(secondaryRadiusMin), secondaryRadiusMax(secondaryRadiusMax),
    tertiaryRadiusMin(tertiaryRadiusMin), tertiaryRadiusMax(tertiaryRadiusMax),
    tertiaryVolumeFraction(tertiaryVolumeFraction),
    targetDensity(targetDensity),
    compactnessFactor(compactnessFactor),
    voxelGrid(size) {
    
    // Initialize spatial index for sphere overlap detection
    SpatialIndex::IStorageManager* memoryManager = 
        SpatialIndex::StorageManager::createNewMemoryStorageManager();
    
    // R-tree configuration for optimal performance
    SpatialIndex::id_type indexIdentifier = 1; 
    SpatialIndex::RTree::RTreeVariant variant = SpatialIndex::RTree::RV_RSTAR;
    double fillFactor = 0.7;
    uint32_t indexCapacity = 10;
    uint32_t leafCapacity = 10;
    uint32_t dimension = 3;
    
    // Create R-tree for sphere queries
    spatialIndex.reset(
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
    
    // Initialize random seed for reproducible results
    srand(static_cast<unsigned int>(time(nullptr)));
}

PackingGenerator::~PackingGenerator() {
    // Smart pointers automatically clean up resources
}

//=============================================================================
// Main Generation Method
//=============================================================================

bool PackingGenerator::generate() {
    std::cout << "Starting packing generation..." << std::endl;
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Domain size: " << size << "×" << size << "×" << size << std::endl;
    std::cout << "  Target density: " << targetDensity << std::endl;
    std::cout << "  Core radius range: [" << coreRadiusMin << ", " << coreRadiusMax << "]" << std::endl;
    std::cout << "  Secondary radius range: [" << secondaryRadiusMin << ", " << secondaryRadiusMax << "]" << std::endl;
    std::cout << "  Tertiary radius range: [" << tertiaryRadiusMin << ", " << tertiaryRadiusMax << "]" << std::endl;
    std::cout << "  Compactness factor: " << compactnessFactor << std::endl;
    std::cout << std::endl;
    
    // Stage 1: Insert core spheres
    std::cout << "Stage 1: Inserting core spheres..." << std::endl;
    auto startTime = std::chrono::steady_clock::now();
    uint32_t coreCount = insertCoreSpheres();
    auto endTime = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    std::cout << "  Inserted " << coreCount << " core spheres in " 
              << duration.count() << " ms" << std::endl;
    std::cout << "  Density after stage 1: " << getCurrentDensity() << std::endl;
    std::cout << std::endl;
    
    // Stage 2: Add secondary spheres
    std::cout << "Stage 2: Adding secondary spheres..." << std::endl;
    startTime = std::chrono::steady_clock::now();
    uint32_t secondaryCount = addSecondarySpheres();
    endTime = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    std::cout << "  Added " << secondaryCount << " secondary spheres in " 
              << duration.count() << " ms" << std::endl;
    std::cout << "  Density after stage 2: " << getCurrentDensity() << std::endl;
    std::cout << std::endl;
    
    // Stage 3: Add tertiary spheres
    std::cout << "Stage 3: Adding tertiary spheres..." << std::endl;
    startTime = std::chrono::steady_clock::now();
    uint32_t tertiaryCount = addTertiarySpheres();
    endTime = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    std::cout << "  Added " << tertiaryCount << " tertiary spheres in " 
              << duration.count() << " ms" << std::endl;
    std::cout << "  Final density: " << getCurrentDensity() << std::endl;
    std::cout << std::endl;
    
    // Update particle contacts
    updateParticleContacts();
    
    // Report final statistics
    std::cout << "Packing generation complete!" << std::endl;
    std::cout << "Final statistics:" << std::endl;
    std::cout << "  Total particles: " << getParticleCount() << std::endl;
    std::cout << "  Total spheres: " << numSpheres << std::endl;
    std::cout << "  Packing density: " << getCurrentDensity() << std::endl;
    std::cout << "  Inter-particle contacts: " << getContactCount() << std::endl;
    std::cout << "  Average coordination number: " << getAverageCoordinationNumber() << std::endl;
    std::cout << "  Average sphericity: " << getAverageSphericity() << std::endl;
    
    return true;
}

//=============================================================================
// Stage 1: Core Sphere Insertion
//=============================================================================

uint32_t PackingGenerator::insertCoreSpheres() {
    uint32_t totalCoreSpheres = 0;
    uint32_t attempts = 0;
    uint32_t maxAttempts = 10000;  // Safety limit to prevent infinite loops
    uint32_t consecutiveFailures = 0;
    uint32_t maxConsecutiveFailures = 1000;  // Stop if too many failures in a row
    
    while (attempts < maxAttempts && consecutiveFailures < maxConsecutiveFailures) {
        // Generate random position
        int x = rand() % size;
        int y = rand() % size;
        int z = rand() % size;
        Point3D center(x, y, z);
        
        // Generate random radius
        int radius = getRandomRadius(coreRadiusMin, coreRadiusMax);
        
        // Check if sphere would be within bounds
        if (x - radius < 0 || x + radius >= static_cast<int>(size) ||
            y - radius < 0 || y + radius >= static_cast<int>(size) ||
            z - radius < 0 || z + radius >= static_cast<int>(size)) {
            attempts++;
            consecutiveFailures++;
            continue;
        }
        
        // Query spatial index for potential overlaps
        double centerCoords[3] = {
            static_cast<double>(x), 
            static_cast<double>(y), 
            static_cast<double>(z)
        };
        SpatialIndex::Ball queryBall(radius + coreRadiusMax, centerCoords, 3);
        
        SphereVisitor visitor;
        spatialIndex->intersectsWithQuery(queryBall, visitor);
        
        bool hasExcessiveOverlap = false;
        
        // Check each found sphere for actual overlap
        for (const auto* foundSphere : visitor.foundSpheres) {
            double distance = center.distanceTo(foundSphere->getCenter());
            double minDistance = radius + foundSphere->getRadius() - 5; // Small tolerance
            
            if (distance < minDistance) {
                hasExcessiveOverlap = true;
                break;
            }
        }
        
        if (hasExcessiveOverlap) {
            attempts++;
            consecutiveFailures++;
            continue;
        }
        
        // Create new particle
        uint16_t particleId = particles.size() + 1;
        particles.emplace_back(particleId);
        
        // Add core sphere to the particle
        uint32_t sphereId = numSpheres++;
        const Sphere& sphere = particles.back().addSphere(center, radius, 
                                                         SphereType::CORE, sphereId);
        
        // Add to voxel grid
        voxelGrid.addSphere(particles, center, radius, particleId);
        
        // Add to spatial index
        SphereData* data = new SphereData(sphereId, &sphere);
        SpatialIndex::Ball sphereBall(radius, centerCoords, 3);
        spatialIndex->insertData(sizeof(data), reinterpret_cast<uint8_t*>(&data), 
                               sphereBall, sphereId);
        
        totalCoreSpheres++;
        consecutiveFailures = 0;  // Reset consecutive failures on success
        attempts++;
        
        // Check if we've reached target density for core spheres
        double density = getCurrentDensity();
        if (density > 0.7 * targetDensity) {  // Core spheres should achieve ~70% of target
            break;
        }
        
        // Report progress periodically
        if (totalCoreSpheres % 50 == 0) {
            std::cout << "    Progress: " << totalCoreSpheres << " core spheres, "
                     << "density: " << density << std::endl;
        }
    }
    
    if (consecutiveFailures >= maxConsecutiveFailures) {
        std::cout << "    Core sphere insertion stopped due to consecutive failures" << std::endl;
    }
    
    return totalCoreSpheres;
}

/**
 * @file PackingGenerator.cpp (Part 2)
 * @brief Implementation of the PackingGenerator class - Continued
 */

//=============================================================================
// Stage 2: Secondary Sphere Addition
//=============================================================================

uint32_t PackingGenerator::addSecondarySpheres(uint32_t maxAttempts) {
    uint32_t attempts = 0;
    uint32_t successfulInsertions = 0;
    uint32_t consecutiveFailures = 0;
    uint32_t maxConsecutiveFailures = 500;
    
    while (attempts < maxAttempts && consecutiveFailures < maxConsecutiveFailures) {
        // Check if we have particles to work with
        if (particles.empty()) {
            break;
        }
        
        // Choose random particle
        uint32_t particleIndex = rand() % particles.size();
        Particle& particle = particles[particleIndex];
        
        // Get the core sphere of this particle
        const Sphere* coreSphere = particle.getCoreSphere();
        if (!coreSphere) {
            attempts++;
            consecutiveFailures++;
            continue;
        }
        
        // Generate position for secondary sphere
        double angle1 = (rand() % 360) * M_PI / 180.0;
        double angle2 = (rand() % 360) * M_PI / 180.0;
        
        int secondaryRadius = getRandomRadius(secondaryRadiusMin, secondaryRadiusMax);
        
        // Distance from core should allow for proper clustering
        int maxDistance = static_cast<int>(2 * coreSphere->getRadius());
        int distanceFromCore = rand() % maxDistance;
        
        // Convert spherical to Cartesian coordinates
        int dx = static_cast<int>(distanceFromCore * sin(angle1) * cos(angle2));
        int dy = static_cast<int>(distanceFromCore * sin(angle1) * sin(angle2));
        int dz = static_cast<int>(distanceFromCore * cos(angle1));
        
        Point3D center(
            coreSphere->getCenter().x + dx,
            coreSphere->getCenter().y + dy,
            coreSphere->getCenter().z + dz
        );
        
        // Check bounds
        if (center.x - secondaryRadius < 0 || center.x + secondaryRadius >= static_cast<int>(size) ||
            center.y - secondaryRadius < 0 || center.y + secondaryRadius >= static_cast<int>(size) ||
            center.z - secondaryRadius < 0 || center.z + secondaryRadius >= static_cast<int>(size)) {
            attempts++;
            consecutiveFailures++;
            continue;
        }
        
        // Create temporary sphere for testing
        uint32_t sphereId = numSpheres;
        Sphere newSphere(center, secondaryRadius, SphereType::SECONDARY, particle.getId(), sphereId);
        
        // Verify compactness criteria with existing spheres in the particle
        bool meetsCriteria = false;
        
        // Should be overlapping the core sphere compactly
        if (newSphere.intersectsWith(*coreSphere, compactnessFactor)) {
            meetsCriteria = true;
        } 
        // Otherwise, should be at least connected to core and compactly to any other secondary spheres
        else if (newSphere.intersectsWith(*coreSphere, 1)){
            for (const auto& existingSphere : particle.getSpheres()) {
                if (existingSphere.getType() == SphereType::SECONDARY && 
                    newSphere.intersectsWith(existingSphere, compactnessFactor)) {
                    meetsCriteria = true;
                    break;
                }
            }
        }
        
        if (!meetsCriteria) {
            attempts++;
            consecutiveFailures++;
            continue;
        }
        
        // Check for overlap with other particles
        double sphereCenterCoords[3] = {
            static_cast<double>(center.x), 
            static_cast<double>(center.y), 
            static_cast<double>(center.z)
        };
        SpatialIndex::Ball queryBall(secondaryRadius + secondaryRadiusMax, sphereCenterCoords, 3);
        
        SphereVisitor visitor;
        spatialIndex->intersectsWithQuery(queryBall, visitor);
        
        bool hasOverlap = false;
        
        for (const auto* foundSphere : visitor.foundSpheres) {
            if (foundSphere->getParticleId() != particle.getId()) {
                double distance = center.distanceTo(foundSphere->getCenter());
                double minDistance = secondaryRadius + foundSphere->getRadius() -5;
                
                if (distance < minDistance) {
                    hasOverlap = true;
                    break;
                }
            }
        }
        
        if (hasOverlap) {
            attempts++;
            consecutiveFailures++;
            continue;
        }
        
        // Add the sphere to the particle
        const Sphere& sphere = particle.addSphere(center, secondaryRadius, 
                                                 SphereType::SECONDARY, sphereId);
        
        // Add to voxel grid
        voxelGrid.addSphere(particles, center, secondaryRadius, particle.getId());
        
        // Add to spatial index
        SphereData* data = new SphereData(sphereId, &sphere);
        SpatialIndex::Ball sphereBall(secondaryRadius, sphereCenterCoords, 3);
        spatialIndex->insertData(sizeof(data), reinterpret_cast<uint8_t*>(&data), 
                               sphereBall, sphereId);
        
        numSpheres++;
        successfulInsertions++;
        consecutiveFailures = 0;
        
        // Check density progress
        double density = getCurrentDensity();
        if (density > 0.85 * targetDensity) {  // Secondary spheres should reach ~85% of target
            break;
        }
        
        // Report progress
        if (successfulInsertions % 100 == 0) {
            std::cout << "    Progress: " << successfulInsertions << " secondary spheres, "
                     << "density: " << density << std::endl;
        }
        
        attempts++;
    }
    
    return successfulInsertions;
}

//=============================================================================
// Stage 3: Tertiary Sphere Addition
//=============================================================================

uint32_t PackingGenerator::addTertiarySpheres(uint32_t maxAttempts) {
    uint32_t attempts = 0;
    uint32_t successfulInsertions = 0;
    uint32_t consecutiveFailures = 0;
    uint32_t maxConsecutiveFailures = 500;
    
    while (attempts < maxAttempts && consecutiveFailures < maxConsecutiveFailures) {
        if (particles.empty()) {
            break;
        }
        
        // Choose random particle
        uint32_t particleIndex = rand() % particles.size();
        Particle& particle = particles[particleIndex];
        
        // Get spheres from the particle
        const std::vector<Sphere>& particleSpheres = particle.getSpheres();
        if (particleSpheres.empty()) {
            attempts++;
            consecutiveFailures++;
            continue;
        }
        
        // Choose a random sphere as base
        uint32_t sphereIndex = rand() % particleSpheres.size();
        const Sphere& baseSphere = particleSpheres[sphereIndex];
        
        // Generate position for tertiary sphere on the surface
        double angle1 = (rand() % 360) * M_PI / 180.0;
        double angle2 = (rand() % 360) * M_PI / 180.0;
        
        int tertiaryRadius = getRandomRadius(tertiaryRadiusMin, tertiaryRadiusMax);
        
        // Place on the surface of the base sphere
        int distanceFromCenter = baseSphere.getRadius() - 
                               static_cast<int>(compactnessFactor * tertiaryRadius);
        
        // Convert to Cartesian
        int dx = static_cast<int>(distanceFromCenter * sin(angle1) * cos(angle2));
        int dy = static_cast<int>(distanceFromCenter * sin(angle1) * sin(angle2));
        int dz = static_cast<int>(distanceFromCenter * cos(angle1));
        
        Point3D center(
            baseSphere.getCenter().x + dx,
            baseSphere.getCenter().y + dy,
            baseSphere.getCenter().z + dz
        );
        
        // Check bounds
        if (center.x - tertiaryRadius < 0 || center.x + tertiaryRadius >= static_cast<int>(size) ||
            center.y - tertiaryRadius < 0 || center.y + tertiaryRadius >= static_cast<int>(size) ||
            center.z - tertiaryRadius < 0 || center.z + tertiaryRadius >= static_cast<int>(size)) {
            attempts++;
            consecutiveFailures++;
            continue;
        }
        
        // Create temporary sphere
        uint32_t sphereId = numSpheres;
        Sphere newSphere(center, tertiaryRadius, SphereType::TERTIARY, particle.getId(), sphereId);
        
        // Verify connection with at least one sphere in the particle
        bool isConnected = false;
        for (const auto& existingSphere : particleSpheres) {
            if (newSphere.intersectsWith(existingSphere, compactnessFactor)) {
                isConnected = true;
                break;
            }
        }
        
        if (!isConnected) {
            attempts++;
            consecutiveFailures++;
            continue;
        }
        
        // Check for opportunities to create inter-particle contacts
        double sphereCenterCoords[3] = {
            static_cast<double>(center.x), 
            static_cast<double>(center.y), 
            static_cast<double>(center.z)
        };
        SpatialIndex::Ball queryBall(tertiaryRadius + tertiaryRadiusMax, sphereCenterCoords, 3);
        
        SphereVisitor visitor;
        spatialIndex->intersectsWithQuery(queryBall, visitor);
        
        bool hasExcessiveOverlap = false;
        bool createsContact = false;
                
        // Add the sphere
        const Sphere& sphere = particle.addSphere(center, tertiaryRadius, 
                                                 SphereType::TERTIARY, sphereId);
        
        // Add to voxel grid (will automatically detect contacts)
        voxelGrid.addSphere(particles, center, tertiaryRadius, particle.getId());
        
        // Add to spatial index
        SphereData* data = new SphereData(sphereId, &sphere);
        SpatialIndex::Ball sphereBall(tertiaryRadius, sphereCenterCoords, 3);
        spatialIndex->insertData(sizeof(data), reinterpret_cast<uint8_t*>(&data), 
                               sphereBall, sphereId);
        
        numSpheres++;
        successfulInsertions++;
        consecutiveFailures = 0;
        
        // Check density
        double density = getCurrentDensity();
        if (density > targetDensity) {
            break;
        }
        
        // Report progress
        if (successfulInsertions % 200 == 0) {
            std::cout << "    Progress: " << successfulInsertions << " tertiary spheres, "
                     << "density: " << density << std::endl;
        }
        
        attempts++;
    }
    
    return successfulInsertions;
}

//=============================================================================
// Property Calculation Methods
//=============================================================================

double PackingGenerator::getCurrentDensity() const {
    double domainVolume = static_cast<double>(size * size * size);
    double filledVolume = static_cast<double>(voxelGrid.getFilledVoxelCount());
    return filledVolume / domainVolume;
}

const Particle* PackingGenerator::getParticle(uint32_t index) const {
    if (index < particles.size()) {
        return &particles[index];
    }
    return nullptr;
}

uint32_t PackingGenerator::getContactCount() const {
    return voxelGrid.getInterfaceSegmentCount();
}

double PackingGenerator::getAverageCoordinationNumber() const {
    if (particles.empty()) {
        return 0.0;
    }
    
    uint32_t totalContacts = 0;
    for (const auto& particle : particles) {
        totalContacts += particle.getCoordinationNumber();
    }
    
    // Each contact is counted twice (once per particle)
    return static_cast<double>(totalContacts) / particles.size();
}

std::vector<uint32_t> PackingGenerator::getCoordinationNumbers() const {
    std::vector<uint32_t> result;
    result.reserve(particles.size());
    
    for (const auto& particle : particles) {
        result.push_back(particle.getCoordinationNumber());
    }
    
    return result;
}

double PackingGenerator::getAverageSphericity() const {
    if (particles.empty()) {
        return 0.0;
    }
    
    double totalSphericity = 0.0;
    for (const auto& particle : particles) {
        totalSphericity += particle.calculateSphericity();
    }
    
    return totalSphericity / particles.size();
}

//=============================================================================
// File Output Method
//=============================================================================

bool PackingGenerator::saveTIFF(const std::string& filename, bool binary) const {
    return voxelGrid.saveToTIFF(filename, binary);
}

//=============================================================================
// Utility Methods
//=============================================================================

Point3D PackingGenerator::getRandomPoint() const {
    return Point3D(
        rand() % size,
        rand() % size,
        rand() % size
    );
}

int PackingGenerator::getRandomRadius(int minRadius, int maxRadius) const {
    return minRadius + rand() % (maxRadius - minRadius + 1);
}

void PackingGenerator::updateParticleContacts() {
    // Contact tracking is now handled automatically by the voxel grid
    // during sphere addition, so this method is essentially a no-op.
    // It's kept for API compatibility and potential future enhancements.
}

} // namespace Packing
