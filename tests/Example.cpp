/**
 * @file Example.cpp
 * @brief Example program demonstrating the Random Packing Generator
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This example shows how to use the PackingGenerator class to create
 * a random packing of non-spherical particles and analyze its properties.
 */

#include "PackingGenerator.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>

using namespace Packing;

/**
 * @brief Saves particle statistics to a CSV file
 * @param filename Output filename
 * @param particles Vector of particles to analyze
 * @return true if successful
 */
bool saveParticleStats(const std::string& filename, const std::vector<Particle>& particles) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to create statistics file: " << filename << std::endl;
        return false;
    }
    
    // Write header
    file << "ParticleID,NumSpheres,Volume,SurfaceArea,Sphericity,CoordinationNumber\n";
    
    // Write data for each particle
    for (const auto& particle : particles) {
        file << particle.getId() << ","
             << particle.getSpheres().size() << ","
             << particle.getVolume() << ","
             << particle.getArea() << ","
             << particle.calculateSphericity() << ","
             << particle.getCoordinationNumber() << "\n";
    }
    
    file.close();
    return true;
}

/**
 * @brief Main example program
 */
int main(int argc, char* argv[]) {
    std::cout << "Random Packing Generator Example" << std::endl;
    std::cout << "================================\n" << std::endl;
    
    // Set random seed based on current time
    srand(static_cast<unsigned int>(time(nullptr)));
    
    // Configuration parameters
    uint32_t domainSize = 300;
    int coreRadiusMin = 30;
    int coreRadiusMax = 40;
    int secondaryRadiusMin = 20;
    int secondaryRadiusMax = 30;
    int tertiaryRadiusMin = 5;
    int tertiaryRadiusMax = 10;
    double tertiaryVolumeFraction = 0.1;
    double targetDensity = 0.65;
    double compactnessFactor = 0.5;
    
    // Parse command line arguments if provided
    if (argc > 1) {
        domainSize = std::atoi(argv[1]);
        std::cout << "Using domain size from command line: " << domainSize << std::endl;
    }
    if (argc > 2) {
        targetDensity = std::atof(argv[2]);
        std::cout << "Using target density from command line: " << targetDensity << std::endl;
    }
    
    // Create the packing generator
    std::cout << "\nCreating packing generator with parameters:" << std::endl;
    std::cout << "  Domain size: " << domainSize << "×" << domainSize << "×" << domainSize << std::endl;
    std::cout << "  Target density: " << targetDensity << std::endl;
    std::cout << "  Core radius: [" << coreRadiusMin << ", " << coreRadiusMax << "]" << std::endl;
    std::cout << "  Secondary radius: [" << secondaryRadiusMin << ", " << secondaryRadiusMax << "]" << std::endl;
    std::cout << "  Tertiary radius: [" << tertiaryRadiusMin << ", " << tertiaryRadiusMax << "]" << std::endl;
    std::cout << "  Compactness factor: " << compactnessFactor << std::endl;
    std::cout << std::endl;
    
    PackingGenerator generator(
        domainSize,
        coreRadiusMin, coreRadiusMax,
        secondaryRadiusMin, secondaryRadiusMax,
        tertiaryRadiusMin, tertiaryRadiusMax,
        targetDensity,
        compactnessFactor
    );
    
    // Generate the packing
    std::cout << "Generating packing...\n" << std::endl;
    bool success = generator.generate();
    
    if (!success) {
        std::cerr << "Failed to generate packing!" << std::endl;
        return 1;
    }
    
    // Report statistics
    std::cout << "\nPacking Statistics:" << std::endl;
    std::cout << "==================" << std::endl;
    std::cout << "Number of particles: " << generator.getParticleCount() << std::endl;
    std::cout << "Packing density: " << generator.getCurrentDensity() << std::endl;
    std::cout << "Number of contacts: " << generator.getContactCount() << std::endl;
    std::cout << "Average coordination number: " << generator.getAverageCoordinationNumber() << std::endl;
    std::cout << "Average sphericity: " << generator.getAverageSphericity() << std::endl;
    
    // Save outputs
    std::cout << "\nSaving outputs..." << std::endl;
    
    // Save binary TIFF
    std::string binaryFilename = "packing_binary.tiff";
    if (generator.saveTIFF(binaryFilename, true)) {
        std::cout << "  Saved binary TIFF: " << binaryFilename << std::endl;
    } else {
        std::cerr << "  Failed to save binary TIFF" << std::endl;
    }
    
    // Save particle ID TIFF
    std::string idFilename = "packing_ids.tiff";
    if (generator.saveTIFF(idFilename, false)) {
        std::cout << "  Saved particle ID TIFF: " << idFilename << std::endl;
    } else {
        std::cerr << "  Failed to save particle ID TIFF" << std::endl;
    }
    
    // Save particle statistics
    std::string statsFilename = "particle_stats.csv";
    if (saveParticleStats(statsFilename, generator.getParticles())) {
        std::cout << "  Saved particle statistics: " << statsFilename << std::endl;
    } else {
        std::cerr << "  Failed to save particle statistics" << std::endl;
    }
    
    // Create coordination number histogram
    std::cout << "\nCoordination Number Distribution:" << std::endl;
    std::vector<uint32_t> coordNumbers = generator.getCoordinationNumbers();
    std::vector<int> histogram(20, 0);  // Assuming max coordination < 20
    
    for (uint32_t coord : coordNumbers) {
        if (coord < histogram.size()) {
            histogram[coord]++;
        }
    }
    
    // Display histogram
    for (size_t i = 0; i < histogram.size(); ++i) {
        if (histogram[i] > 0) {
            std::cout << "  " << i << " contacts: " << histogram[i] << " particles" << std::endl;
        }
    }
    
    std::cout << "\nExample completed successfully!" << std::endl;
    return 0;
}