/**
 * @file binding.cpp
 * @brief Python bindings for the Random Packing Generator library using pybind11
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "PackingGenerator.h"

namespace py = pybind11;

PYBIND11_MODULE(particle_packing, m) {
    m.doc() = "Python bindings for the Random Sequential Addition-based Particle Packing Generator";

    // Enums
    py::enum_<Packing::SphereType>(m, "SphereType")
        .value("CORE", Packing::SphereType::CORE)
        .value("SECONDARY", Packing::SphereType::SECONDARY)
        .value("TERTIARY", Packing::SphereType::TERTIARY)
        .export_values();

    py::enum_<Packing::VoxelType>(m, "VoxelType")
        .value("SURFACE", Packing::VoxelType::SURFACE)
        .value("BULK", Packing::VoxelType::BULK)
        .value("AIR", Packing::VoxelType::AIR)
        .value("GAP", Packing::VoxelType::GAP)
        .export_values();

    // Point3D
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

    // Sphere
    py::class_<Packing::Sphere>(m, "Sphere")
        .def(py::init<const Packing::Point3D&, int, Packing::SphereType, uint16_t, uint32_t>())
        .def("getCenter", &Packing::Sphere::getCenter)
        .def("getRadius", &Packing::Sphere::getRadius)
        .def("getType", &Packing::Sphere::getType)
        .def("getParticleId", &Packing::Sphere::getParticleId)
        .def("getSphereId", &Packing::Sphere::getSphereId)
        .def("intersectsWith", &Packing::Sphere::intersectsWith, py::arg("other"), py::arg("beta") = 0.5)
        .def("__repr__", [](const Packing::Sphere &s) {
            std::string typeStr;
            switch (s.getType()) {
                case Packing::SphereType::CORE: typeStr = "CORE"; break;
                case Packing::SphereType::SECONDARY: typeStr = "SECONDARY"; break;
                case Packing::SphereType::TERTIARY: typeStr = "TERTIARY"; break;
            }
            return "Sphere(center=" + std::to_string(s.getCenter().x) + "," + 
                   std::to_string(s.getCenter().y) + "," + std::to_string(s.getCenter().z) + 
                   ", radius=" + std::to_string(s.getRadius()) + 
                   ", type=" + typeStr + 
                   ", particleId=" + std::to_string(s.getParticleId()) + 
                   ", sphereId=" + std::to_string(s.getSphereId()) + ")";
        });

    // Particle
    py::class_<Packing::Particle>(m, "Particle")
        .def(py::init<uint16_t>())
        .def("getId", &Packing::Particle::getId)
        .def("getSpheres", &Packing::Particle::getSpheres)
        .def("getCoreSphere", &Packing::Particle::getCoreSphere)
        .def("getVolume", &Packing::Particle::getVolume)
        .def("getArea", &Packing::Particle::getArea)
        .def("incrementVolume", &Packing::Particle::incrementVolume)
        .def("decrementVolume", &Packing::Particle::decrementVolume)
        .def("incrementArea", &Packing::Particle::incrementArea)
        .def("decrementArea", &Packing::Particle::decrementArea)
        .def("calculateSphericity", &Packing::Particle::calculateSphericity)
        .def("getContacts", &Packing::Particle::getContacts)
        .def("addContact", &Packing::Particle::addContact)
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
        .def(py::init<uint32_t, int, int, int, int, int, int, double, double, double>(),
             py::arg("size"),
             py::arg("coreRadiusMin"), py::arg("coreRadiusMax"),
             py::arg("secondaryRadiusMin"), py::arg("secondaryRadiusMax"),
             py::arg("tertiaryRadiusMin"), py::arg("tertiaryRadiusMax"),
             py::arg("tertiaryVolumeFraction"),
             py::arg("targetDensity"),
             py::arg("compactnessFactor") = 0.5)
        .def("generate", &Packing::PackingGenerator::generate)
        .def("getCurrentDensity", &Packing::PackingGenerator::getCurrentDensity)
        .def("getParticle", &Packing::PackingGenerator::getParticle, py::return_value_policy::reference)
        .def("getParticleCount", &Packing::PackingGenerator::getParticleCount)
        .def("getParticles", &Packing::PackingGenerator::getParticles, py::return_value_policy::reference)
        .def("getContactCount", &Packing::PackingGenerator::getContactCount)
        .def("getAverageCoordinationNumber", &Packing::PackingGenerator::getAverageCoordinationNumber)
        .def("getCoordinationNumbers", &Packing::PackingGenerator::getCoordinationNumbers)
        .def("getAverageSphericity", &Packing::PackingGenerator::getAverageSphericity)
        .def("saveTIFF", &Packing::PackingGenerator::saveTIFF, 
             py::arg("filename"), py::arg("binary") = true)
        .def("__repr__", [](const Packing::PackingGenerator &pg) {
            return "PackingGenerator(particles=" + std::to_string(pg.getParticleCount()) + 
                   ", density=" + std::to_string(pg.getCurrentDensity()) + 
                   ", contacts=" + std::to_string(pg.getContactCount()) + 
                   ", avgCoordination=" + std::to_string(pg.getAverageCoordinationNumber()) + 
                   ", avgSphericity=" + std::to_string(pg.getAverageSphericity()) + ")";
        });

    // Create visualization function for debugging
    m.def("visualize_packing", [](const Packing::PackingGenerator& generator, const std::string& outputFile) {
        // Simple visualization function to create a 3D NumPy array of the packing
        // Returns a 3D binary array where 1 indicates filled voxels
        uint32_t size = generator.getParticles().empty() ? 0 : 
                        generator.getParticles()[0].getCoreSphere()->getRadius() * 10;
        if (size == 0) return py::array_t<uint8_t>();
        
        // Save TIFF file and let user know they can visualize it with external tools
        generator.saveTIFF(outputFile, true);
        
        std::cout << "Packing saved to " << outputFile << std::endl;
        std::cout << "Use external tools like ImageJ or Python's tifffile to visualize the 3D data." << std::endl;
        
        // Return an empty array (visualization handled externally)
        return py::array_t<uint8_t>();
    }, "Create a visualization of the packing and save to a TIFF file",
       py::arg("generator"), py::arg("outputFile") = "packing.tiff");
}