/**
 * @file Sphere.cpp
 * @brief Implementation of the Sphere class
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file contains the implementation of the Sphere class methods.
 */

#include "Sphere.h"
#include <algorithm>

namespace Packing {

//=============================================================================
// Sphere Implementation
//=============================================================================

Sphere::Sphere(const Point3D& center, int radius, SphereType type, 
               uint16_t particleId, uint32_t sphereId)
    : center(center), radius(radius), type(type), 
      particleId(particleId), sphereId(sphereId) {
}

bool Sphere::intersectsWith(const Sphere& other, double beta) const {
    // Calculate center-to-center distance
    double distance = center.distanceTo(other.center);
    
    // Ensure the larger radius is r2 and smaller is r1
    double r1 = std::min(static_cast<double>(radius), 
                        static_cast<double>(other.radius));
    double r2 = std::max(static_cast<double>(radius), 
                        static_cast<double>(other.radius));
    
    // Apply compactness criteria from the paper:
    // Condition 1: Spheres must be tightly packed (not inside each other)
    // Condition 2: A newly added sphere should lead to particle growth
    // Formula: r2 - r1 < d ≤ r2 - (β × r1)
    bool meetsCriteria = (r2 - r1 < distance && distance <= r2 - (beta * r1));
    
    return meetsCriteria;
}

bool Sphere::intersectsWith(const double x, const double y, const double z, const double radiusOther, double beta) const {
    
    Point3D centerOther(x, y, z);

    // Calculate center-to-center distance
    double distance = center.distanceTo(centerOther);
    
    // Ensure the larger radius is r2 and smaller is r1
    double r1 = std::min(radiusOther, 
                        static_cast<double>(radius));
    double r2 = std::max(radiusOther, 
                        static_cast<double>(radius));
    
    // Apply compactness criteria from the paper:
    // Condition 1: Spheres must be tightly packed (not inside each other)
    // Condition 2: A newly added sphere should lead to particle growth
    // Formula: r2 - r1 < d ≤ r2 - (β × r1)
    bool meetsCriteria = (r2 - r1 < distance && distance <= r2 - (beta * r1));
    
    return meetsCriteria;
}

} // namespace Packing