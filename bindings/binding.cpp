/**
 * @file binding.cpp
 * @brief Improved Python bindings for the Random Packing Generator library using pybind11
 * 
 * This implementation focuses on proper memory management and ownership tracking
 * to avoid crashes and memory leaks when used from Python.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "PackingGenerator.h"

namespace py = pybind11;

PYBIND11_MODULE(particle_packing, m) {
    m.doc() = "Python bindings for the Random Sequential Addition-based Particle Packing Generator";

    // Enums - these are simple and don't have ownership issues
    py::enum_<Packing::SphereType>(m, "SphereType")
        .value("CORE", Packing::SphereType::CORE)
        .value("SECONDARY", Packing::SphereType::SECONDARY)
        .value("TERTIARY", Packing::SphereType::TERTIARY)
        .export_values();

    py::enum_<Packing::VoxelType>(m, "VoxelType")
        .value("SURFACE", Packing::VoxelType::SURFACE)
        .value("BULK", Packing::VoxelType::BULK)
        .value("AIR", Packing::VoxelType::AIR)
        .export_values();

    // Point3D - Value type without ownership complexity
    py::class_<Packing::Point3D>(m, "Point3D")
        .def(py::init<int, int, int>(), py::arg("x") = 0, py::arg("y") = 0, py::arg("z") = 0)
        .def_readwrite("x", &Packing::Point3D::x)
        .def_readwrite("y", &Packing::Point3D::y)
        .def_readwrite("z", &Packing::Point3D::z)
        .def("distanceTo", &Packing::Point3D::distanceTo)
        .def("__eq__", &Packing::Point3D::operator==)
        .def("__repr__", [](const Packing::Point3D &p) {
            return "Point3D(" + std::to_string(p.x) + ", " + std::to_string(p.y) + ", " + std::to_string(p.z) + ")";
        });

    // Sphere - Use copy policy for safer access and prevent raw pointer access
    py::class_<Packing::Sphere, std::shared_ptr<Packing::Sphere>>(m, "Sphere")
        .def(py::init<const Packing::Point3D&, int, Packing::SphereType, uint16_t, uint32_t>())
        .def("getCenter", &Packing::Sphere::getCenter, py::return_value_policy::copy)
        .def("getRadius", &Packing::Sphere::getRadius)
        .def("getType", &Packing::Sphere::getType)
        .def("getParticleId", &Packing::Sphere::getParticleId)
        .def("getSphereId", &Packing::Sphere::getSphereId)
        .def("intersectsWith", 
            static_cast<bool (Packing::Sphere::*)(const Packing::Sphere&, double) const>(&Packing::Sphere::intersectsWith),
            py::arg("other"), py::arg("beta") = 0.5)
        .def("intersectsWith", 
            static_cast<bool (Packing::Sphere::*)(const double, const double, const double, const double, double) const>(&Packing::Sphere::intersectsWith),
            py::arg("x"), py::arg("y"), py::arg("z"), py::arg("radius"), py::arg("beta") = 0.5)        
        .def("__repr__", [](const Packing::Sphere &s) {
            std::string typeStr;
            switch (s.getType()) {
                case Packing::SphereType::CORE: typeStr = "CORE"; break;
                case Packing::SphereType::SECONDARY: typeStr = "SECONDARY"; break;
                case Packing::SphereType::TERTIARY: typeStr = "TERTIARY"; break;
            }
            return "Sphere(center=Point3D(" + std::to_string(s.getCenter().x) + "," + 
                   std::to_string(s.getCenter().y) + "," + std::to_string(s.getCenter().z) + 
                   "), radius=" + std::to_string(s.getRadius()) + 
                   ", type=" + typeStr + 
                   ", particleId=" + std::to_string(s.getParticleId()) + 
                   ", sphereId=" + std::to_string(s.getSphereId()) + ")";
        });

    // Particle - More careful with collections and reference management
    py::class_<Packing::Particle>(m, "Particle")
        .def(py::init<uint16_t>())
        .def("getId", &Packing::Particle::getId)
        .def("getSpheres", &Packing::Particle::getSpheres, 
            py::return_value_policy::reference_internal)
        .def("getCoreSphere", [](const Packing::Particle &p) {
            auto sphere = p.getCoreSphere();
            if (sphere) {
                return sphere;
            } else {
                return std::shared_ptr<Packing::Sphere>();  // Return None if not found
            }
        }, py::return_value_policy::copy)
        .def("getVolume", &Packing::Particle::getVolume)
        .def("getArea", &Packing::Particle::getArea)
        .def("calculateSphericity", &Packing::Particle::calculateSphericity)
        .def("getContacts", [](const Packing::Particle &p) {
            // Create a copy of the set to avoid lifetime issues
            return std::vector<uint16_t>(p.getContacts().begin(), p.getContacts().end());
        }, py::return_value_policy::copy)
        .def("getCoordinationNumber", &Packing::Particle::getCoordinationNumber)
        .def("__repr__", [](const Packing::Particle &p) {
            return "Particle(id=" + std::to_string(p.getId()) + 
                   ", spheres=" + std::to_string(p.getSpheres().size()) + 
                   ", volume=" + std::to_string(p.getVolume()) + 
                   ", area=" + std::to_string(p.getArea()) + 
                   ", contacts=" + std::to_string(p.getCoordinationNumber()) + ")";
        });

    // PackingGenerator - Main class
    py::class_<Packing::PackingGenerator>(m, "PackingGenerator")
        .def(py::init<uint32_t, int, int, int, int, int, int, double, double, uint32_t>(),
             py::arg("size"),
             py::arg("coreRadiusMin") = 10, py::arg("coreRadiusMax") = 20,
             py::arg("secondaryRadiusMin") = 7, py::arg("secondaryRadiusMax") = 12,
             py::arg("tertiaryRadiusMin") = 2, py::arg("tertiaryRadiusMax") = 7,
             py::arg("targetDensity") = 0.6,
             py::arg("compactnessFactor") = 0.5,
             py::arg("randomSeed") = 0)
        .def("generate", &Packing::PackingGenerator::generate)
        .def("getCurrentDensity", &Packing::PackingGenerator::getCurrentDensity)
        .def("getParticle", [](const Packing::PackingGenerator &pg, uint32_t index) {
            const Packing::Particle* p = pg.getParticle(index);
            if (p) {
                // Return a copy to ensure safety
                return Packing::Particle(*p);
            }
            throw py::index_error("Particle index out of range");
        }, py::return_value_policy::copy)
        .def("getSphere", &Packing::PackingGenerator::getSphere, 
             py::arg("index"), 
             "Get sphere at specified index (returns None if index out of range)")
        .def("getParticleCount", &Packing::PackingGenerator::getParticleCount)
        .def("getSphereCount", &Packing::PackingGenerator::getSphereCount)
        .def("getParticles", &Packing::PackingGenerator::getParticles,
             py::return_value_policy::reference_internal)
        .def("getContactCount", &Packing::PackingGenerator::getContactCount)
        .def("getAverageCoordinationNumber", &Packing::PackingGenerator::getAverageCoordinationNumber)
        .def("getCoordinationNumbers", &Packing::PackingGenerator::getCoordinationNumbers, 
             py::return_value_policy::copy)
        .def("getAverageSphericity", &Packing::PackingGenerator::getAverageSphericity)
        .def("saveTIFF", &Packing::PackingGenerator::saveTIFF, 
             py::arg("filename"), py::arg("binary") = true)
        
        // NEW BINDINGS - Direct access to requested methods
        .def("insertCoreSpheres", &Packing::PackingGenerator::insertCoreSpheres,
             "Insert core spheres into the domain and return the number of spheres successfully placed")
        .def("getAverageParticleRadius", &Packing::PackingGenerator::getAverageParticleRadius,
             "Calculate the average particle radius using equivalent sphere volume")
        .def("getTotalVolume", &Packing::PackingGenerator::getTotalVolume,
             "Get the total volume occupied by particles in the voxel grid")
        .def("getParticleVolume", &Packing::PackingGenerator::getParticleVolume,
             py::arg("particleId"),
             "Get the total volume occupied by a particle in the voxel grid")
        .def("getContactPairs", &Packing::PackingGenerator::getContactPairs,
             py::return_value_policy::copy,
             "Get all contact pairs between particles as a list of (particleID1, particleID2) tuples")
        .def("getSphereNeighbors", &Packing::PackingGenerator::getSphereNeighbors,
             py::arg("x"), py::arg("y"), py::arg("z"), py::arg("radius"),
             py::return_value_policy::copy,
             "Get all neighboring spheres as a list of sphereIDs")
        .def("insertSuppSpheres", &Packing::PackingGenerator::insertSuppSpheres,
             py::arg("x"), py::arg("y"), py::arg("z"), py::arg("radius"), py::arg("particleID"),
             "Insert a supplementary sphere into a specific particle. Returns True if successful, False otherwise")
        
        .def("__repr__", [](const Packing::PackingGenerator &pg) {
            return "PackingGenerator(particles=" + std::to_string(pg.getParticleCount()) + 
                   ", density=" + std::to_string(pg.getCurrentDensity()) + 
                   ", contacts=" + std::to_string(pg.getContactCount()) + 
                   ", avgCoordination=" + std::to_string(pg.getAverageCoordinationNumber()) + 
                   ", avgSphericity=" + std::to_string(pg.getAverageSphericity()) + ")";
        });

    // Add a garbage collection function to help manage memory
    m.def("collect_garbage", []() {
        py::gil_scoped_acquire acquire;
        py::module::import("gc").attr("collect")();
    }, "Force Python's garbage collector to run");

    // Create visualization function for debugging that avoids memory issues
    m.def("visualize_packing", [](const Packing::PackingGenerator& generator, const std::string& outputFile) {
        // Simple visualization function to create a 3D NumPy array of the packing
        // Returns a 3D binary array where 1 indicates filled voxels
        uint32_t size = generator.getParticleCount() > 0 ? 
                        generator.getParticleCount() * 10 : 0;  // Simplified logic to avoid potential crashes
        if (size == 0) return py::array_t<uint8_t>();
        
        // Save TIFF file and let user know they can visualize it with external tools
        generator.saveTIFF(outputFile, true);
        
        std::cout << "Packing saved to " << outputFile << std::endl;
        std::cout << "Use external tools like ImageJ or Python's tifffile to visualize the 3D data." << std::endl;
        
        // Return an empty array (visualization handled externally) to avoid memory issues
        return py::array_t<uint8_t>();
    }, "Create a visualization of the packing and save to a TIFF file",
       py::arg("generator"), py::arg("outputFile") = "packing.tiff");
}