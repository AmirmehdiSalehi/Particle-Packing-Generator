/**
 * @file Particle.cpp
 * @brief Implementation of the Particle class
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file contains the implementation of the Particle class methods.
 */

#include "Particle.h"
#include <cmath>
#include <algorithm>

namespace Packing {

//=============================================================================
// Particle Implementation
//=============================================================================

Particle::Particle(uint16_t id) 
    : id(id), bulkCount(0), surfaceCount(0) {
}

const std::shared_ptr<Sphere>& Particle::addSphere(const Point3D& center, int radius, 
                                 SphereType type, uint32_t sphereID) {
    // Create and add the sphere to this particle
    spheres.emplace_back(std::make_shared<Sphere>(center, radius, type, id, sphereID)); 
    return spheres.back();                          
}

const std::shared_ptr<Sphere>& Particle::getCoreSphere() const {
    // Search for the core sphere in the collection
    for (const auto& sphere : spheres) {
        if (sphere->getType() == SphereType::CORE) {
            return sphere;
        }
    }
    return nullptr;
}

double Particle::calculateSphericity() const {
    // Get volume and surface area in voxels
    double volume = static_cast<double>(getVolume());
    double area = static_cast<double>(getArea());
    
    // Avoid division by zero
    if (area == 0) {
        return 0.0;
    }
              
    // Apply sphericity formula from the paper
    // Sphericity = (36π × V²) / S³
    double sphericity = (36.0 * M_PI * std::pow(volume, 2)) / std::pow(area, 3);
    
    // Clamp result to valid range [0, 1]
    return std::min(1.0, std::max(0.0, sphericity));
}

void Particle::addContact(uint16_t particleId) {
    // Prevent self-contacts
    if (particleId == id) {
        return;
    }
    
    // Insert into set (automatically handles duplicates)
    contacts.insert(particleId);
}

} // namespace Packing