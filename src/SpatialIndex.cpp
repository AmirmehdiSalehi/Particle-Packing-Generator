/**
 * @file SpatialIndex.cpp
 * @brief Implementation of spatial indexing helper classes
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file contains the implementation of helper classes for spatial
 * indexing using the libspatialindex library.
 */

#include "SpatialIndex.h"

namespace Packing {

//=============================================================================
// SegmentData Implementation
//=============================================================================

SegmentData::SegmentData(uint32_t id, const Interface* segment) 
    : m_id(id), m_segment(segment) {
}

SpatialIndex::id_type SegmentData::getIdentifier() const { 
    return m_id; 
}

void SegmentData::getData(uint32_t& length, uint8_t** data) const {
    length = sizeof(Interface*);
    *data = reinterpret_cast<uint8_t*>(const_cast<Interface**>(&m_segment));
}

SpatialIndex::IShape* SegmentData::getShape() const {
    const auto& bbox = m_segment->getBoundingBox();
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
    return new SpatialIndex::Region(pLow, pHigh, 3);
}

void SegmentData::getShape(SpatialIndex::IShape** shape) const {
    const auto& bbox = m_segment->getBoundingBox();
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
    *shape = new SpatialIndex::Region(pLow, pHigh, 3);
}

//=============================================================================
// SegmentVisitor Implementation
//=============================================================================

void SegmentVisitor::visitNode(const SpatialIndex::INode& node) {
    // Not needed for data collection
}

void SegmentVisitor::visitData(const SpatialIndex::IData& data) {
    uint32_t len;
    uint8_t* ptr;
    
    data.getData(len, &ptr);
    const Interface* segment = *reinterpret_cast<const Interface**>(ptr);
    
    foundSegments.push_back(segment);
}

void SegmentVisitor::visitData(std::vector<const SpatialIndex::IData*>& v) {
    for (const auto* data : v) {
        visitData(*data);
    }
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

SpatialIndex::IShape* SphereData::getShape() const {
    const Point3D& center = m_sphere->getCenter();
    double centerCoords[3] = {
        static_cast<double>(center.x),
        static_cast<double>(center.y),
        static_cast<double>(center.z)
    };
    return new SpatialIndex::Ball(
        static_cast<double>(m_sphere->getRadius()),
        centerCoords,
        3
    );
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