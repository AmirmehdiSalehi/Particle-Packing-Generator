/**
 * @file CApi.cpp
 * @brief Implementation of the C API for the Random Packing Generator
 * @author Amirmehdi Salehi
 * @date 2024
 * @copyright MIT License
 *
 * This file implements the C-compatible interface functions that
 * wrap the C++ PackingGenerator class.
 */

#include "CApi.h"
#include "PackingGenerator.h"
#include <iostream>
#include <cstring>

extern "C" {

//=============================================================================
// Generator Management Functions
//=============================================================================

void* CreatePackingGenerator(
    uint32_t size,
    int coreRadiusMin, int coreRadiusMax,
    int secondaryRadiusMin, int secondaryRadiusMax,
    int tertiaryRadiusMin, int tertiaryRadiusMax,
    double tertiaryVolumeFraction,
    double targetDensity,
    double compactnessFactor
) {
    try {
        auto* generator = new Packing::PackingGenerator(
            size,
            coreRadiusMin, coreRadiusMax,
            secondaryRadiusMin, secondaryRadiusMax,
            tertiaryRadiusMin, tertiaryRadiusMax,
            targetDensity,
            compactnessFactor
        );
        return static_cast<void*>(generator);
    } catch (const std::exception& e) {
        std::cerr << "Error creating packing generator: " << e.what() << std::endl;
        return nullptr;
    }
}

void FreePackingGenerator(void* generator) {
    if (generator) {
        auto* gen = static_cast<Packing::PackingGenerator*>(generator);
        delete gen;
    }
}

//=============================================================================
// Generation Function
//=============================================================================

int Generate(void* generator) {
    if (!generator) {
        return 0;
    }
    
    try {
        auto* gen = static_cast<Packing::PackingGenerator*>(generator);
        return gen->generate() ? 1 : 0;
    } catch (const std::exception& e) {
        std::cerr << "Error generating packing: " << e.what() << std::endl;
        return 0;
    }
}

//=============================================================================
// Property Access Functions
//=============================================================================

double GetDensity(void* generator) {
    if (!generator) {
        return 0.0;
    }
    
    try {
        auto* gen = static_cast<Packing::PackingGenerator*>(generator);
        return gen->getCurrentDensity();
    } catch (const std::exception& e) {
        std::cerr << "Error getting density: " << e.what() << std::endl;
        return 0.0;
    }
}

uint32_t GetParticleCount(void* generator) {
    if (!generator) {
        return 0;
    }
    
    try {
        auto* gen = static_cast<Packing::PackingGenerator*>(generator);
        return gen->getParticleCount();
    } catch (const std::exception& e) {
        std::cerr << "Error getting particle count: " << e.what() << std::endl;
        return 0;
    }
}

uint32_t GetContactCount(void* generator) {
    if (!generator) {
        return 0;
    }
    
    try {
        auto* gen = static_cast<Packing::PackingGenerator*>(generator);
        return gen->getContactCount();
    } catch (const std::exception& e) {
        std::cerr << "Error getting contact count: " << e.what() << std::endl;
        return 0;
    }
}

double GetAverageCoordinationNumber(void* generator) {
    if (!generator) {
        return 0.0;
    }
    
    try {
        auto* gen = static_cast<Packing::PackingGenerator*>(generator);
        return gen->getAverageCoordinationNumber();
    } catch (const std::exception& e) {
        std::cerr << "Error getting coordination number: " << e.what() << std::endl;
        return 0.0;
    }
}

uint32_t* GetCoordinationNumbers(void* generator, uint32_t* result) {
    if (!generator || !result) {
        return nullptr;
    }
    
    try {
        auto* gen = static_cast<Packing::PackingGenerator*>(generator);
        std::vector<uint32_t> numbers = gen->getCoordinationNumbers();
        
        // Copy to pre-allocated array
        for (size_t i = 0; i < numbers.size(); ++i) {
            result[i] = numbers[i];
        }
        
        return result;
    } catch (const std::exception& e) {
        std::cerr << "Error getting coordination numbers: " << e.what() << std::endl;
        return nullptr;
    }
}

double GetAverageSphericity(void* generator) {
    if (!generator) {
        return 0.0;
    }
    
    try {
        auto* gen = static_cast<Packing::PackingGenerator*>(generator);
        return gen->getAverageSphericity();
    } catch (const std::exception& e) {
        std::cerr << "Error getting sphericity: " << e.what() << std::endl;
        return 0.0;
    }
}

//=============================================================================
// File I/O Function
//=============================================================================

int SaveTIFF(void* generator, const char* filename, int binary) {
    if (!generator || !filename) {
        return 0;
    }
    
    try {
        auto* gen = static_cast<Packing::PackingGenerator*>(generator);
        return gen->saveTIFF(filename, binary != 0) ? 1 : 0;
    } catch (const std::exception& e) {
        std::cerr << "Error saving TIFF: " << e.what() << std::endl;
        return 0;
    }
}

} // extern "C"