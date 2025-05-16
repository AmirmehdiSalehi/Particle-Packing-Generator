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

//=============================================================================
// SphereData Implementation
//=============================================================================

SphereData::SphereData(uint32_t id, const Sphere* sphere) 
    : m_id(id), m_sphere(sphere) {
}

SpatialIndex::id_type SphereData::getIdentifier() const { 
    return m_id; 
}

void SphereData::getData(uint32_t& length, uint8_t** data) const {
    length = sizeof(Sphere*);
    *data = reinterpret_cast<uint8_t*>(const_cast<Sphere**>(&m_sphere));
}

// SpatialIndex::IShape* SphereData::getShape() const {
//     const Point3D& center = m_sphere->getCenter();
//     double centerCoords[3] = {
//         static_cast<double>(center.x),
//         static_cast<double>(center.y),
//         static_cast<double>(center.z)
//     };
//     return new SpatialIndex::Ball(
//         static_cast<double>(m_sphere->getRadius()),
//         centerCoords,
//         3
//     );
// }

SphereData* SphereData::clone() {
  return new SphereData(m_id, m_sphere);
}

void SphereData::getShape(SpatialIndex::IShape** shape) const {
    const Point3D& center = m_sphere->getCenter();
    double centerCoords[3] = {
        static_cast<double>(center.x),
        static_cast<double>(center.y),
        static_cast<double>(center.z)
    };
    *shape = new SpatialIndex::Ball(
        static_cast<double>(m_sphere->getRadius()),
        centerCoords,
        3
    );
}


//=============================================================================
// SphereVisitor Implementation
//=============================================================================

void SphereVisitor::visitNode(const SpatialIndex::INode& node) {
    // Not needed for data collection
}

void SphereVisitor::visitData(const SpatialIndex::IData& data) {
    uint32_t len;
    uint8_t* ptr;
    
    data.getData(len, &ptr);
    const Sphere* sphere = *reinterpret_cast<const Sphere**>(ptr);
    
    foundSpheres.push_back(sphere);
}

void SphereVisitor::visitData(std::vector<const SpatialIndex::IData*>& v) {
    for (const auto* data : v) {
        visitData(*data);
    }
}

} // namespace Packing