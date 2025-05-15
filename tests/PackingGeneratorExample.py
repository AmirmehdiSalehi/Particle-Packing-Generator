#!/usr/bin/env python3
"""
Example script for using the Random Packing Generator Python bindings.
This script demonstrates how to create, configure, and analyze particle packings.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Import the particle packing module
try:
    import particle_packing as pp
except ImportError:
    print("Error: particle_packing module not found.")
    print("Make sure you've built the Python bindings using CMake and added them to your Python path.")
    print("You can add them to your path with:")
    print("  import sys; sys.path.append('/path/to/build/python')")
    sys.exit(1)

def create_packing(
    size=300,
    core_radius_range=(30, 40),
    secondary_radius_range=(20, 30),
    tertiary_radius_range=(5, 10),
    tertiary_volume_fraction=0.1,
    target_density=0.65,
    compactness_factor=0.5,
    output_file="packing.tiff"
):
    """
    Create a random packing of non-spherical particles.
    
    Args:
        size: Size of the cubic domain
        core_radius_range: (min, max) range for core sphere radii
        secondary_radius_range: (min, max) range for secondary sphere radii
        tertiary_radius_range: (min, max) range for tertiary sphere radii
        tertiary_volume_fraction: Volume fraction of tertiary spheres
        target_density: Target packing density
        compactness_factor: Factor controlling sphere overlap (0-1)
        output_file: Filename for the output TIFF file
        
    Returns:
        generator: The PackingGenerator object with the generated packing
    """
    print(f"Creating packing with size {size} and target density {target_density}...")
    
    # Create generator
    generator = pp.PackingGenerator(
        size,
        core_radius_range[0], core_radius_range[1],
        secondary_radius_range[0], secondary_radius_range[1],
        tertiary_radius_range[0], tertiary_radius_range[1],
        tertiary_volume_fraction,
        target_density,
        compactness_factor
    )
    
    # Generate the packing
    print("Generating packing... (this may take a while)")
    success = generator.generate()
    
    if not success:
        print("Failed to generate packing.")
        return None
    
    # Print statistics
    print("\nPacking Statistics:")
    print(f"  Particle count: {generator.getParticleCount()}")
    print(f"  Actual density: {generator.getCurrentDensity():.4f}")
    print(f"  Contact count: {generator.getContactCount()}")
    print(f"  Average coordination number: {generator.getAverageCoordinationNumber():.2f}")
    print(f"  Average sphericity: {generator.getAverageSphericity():.4f}")
    
    # Save the packing
    print(f"\nSaving packing to {output_file}...")
    generator.saveTIFF(output_file, True)
    print("Done!")
    
    return generator

def analyze_packing(generator):
    """
    Analyze a generated packing and produce some plots.
    
    Args:
        generator: The PackingGenerator object to analyze
    """
    if generator is None or generator.getParticleCount() == 0:
        print("No valid packing to analyze.")
        return
    
    # Get coordination numbers
    coord_numbers = generator.getCoordinationNumbers()
    
    # Plot coordination number distribution
    plt.figure(figsize=(10, 6))
    bins = np.arange(min(coord_numbers), max(coord_numbers) + 2) - 0.5
    plt.hist(coord_numbers, bins=bins, alpha=0.7, color='blue', edgecolor='black')
    plt.xlabel('Coordination Number')
    plt.ylabel('Frequency')
    plt.title('Particle Coordination Number Distribution')
    plt.grid(alpha=0.3)
    plt.savefig('coordination_distribution.png')
    print("Saved coordination number distribution to coordination_distribution.png")
    
    # Get all particles for analysis
    particles = generator.getParticles()
    
    # Extract particle properties
    volumes = [p.getVolume() for p in particles]
    areas = [p.getArea() for p in particles]
    sphericities = [p.calculateSphericity() for p in particles]
    
    # Plot sphericity vs coordination number
    plt.figure(figsize=(10, 6))
    plt.scatter(sphericities, coord_numbers, alpha=0.5)
    plt.xlabel('Sphericity')
    plt.ylabel('Coordination Number')
    plt.title('Sphericity vs Coordination Number')
    plt.grid(alpha=0.3)
    plt.savefig('sphericity_vs_coordination.png')
    print("Saved sphericity vs coordination plot to sphericity_vs_coordination.png")
    
    # Plot volume vs area
    plt.figure(figsize=(10, 6))
    plt.scatter(volumes, areas, alpha=0.5)
    plt.xlabel('Particle Volume (voxels)')
    plt.ylabel('Particle Surface Area (voxels)')
    plt.title('Particle Volume vs Surface Area')
    plt.grid(alpha=0.3)
    plt.savefig('volume_vs_area.png')
    print("Saved volume vs area plot to volume_vs_area.png")
    
def visualize_particle(particle, ax=None):
    """
    Visualize a single particle by showing its constituent spheres.
    
    Args:
        particle: The Particle object to visualize
        ax: Optional matplotlib 3D axis to plot on
    """
    if ax is None:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
    
    # Draw each sphere in the particle
    for sphere in particle.getSpheres():
        # Get sphere properties
        center = sphere.getCenter()
        radius = sphere.getRadius()
        sphere_type = sphere.getType()
        
        # Color based on sphere type
        if sphere_type == pp.SphereType.CORE:
            color = 'red'
            alpha = 0.7
        elif sphere_type == pp.SphereType.SECONDARY:
            color = 'green'
            alpha = 0.5
        else:  # TERTIARY
            color = 'blue'
            alpha = 0.3
        
        # Create a wireframe sphere
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = center.x + radius * np.cos(u) * np.sin(v)
        y = center.y + radius * np.sin(u) * np.sin(v)
        z = center.z + radius * np.cos(v)
        
        # Plot the sphere
        ax.plot_wireframe(x, y, z, color=color, alpha=alpha)
    
    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Particle {particle.getId()} - Sphericity: {particle.calculateSphericity():.4f}')
    
    # Try to maintain aspect ratio
    max_range = max(
        [max(sphere.getCenter().x + sphere.getRadius() for sphere in particle.getSpheres()) - 
         min(sphere.getCenter().x - sphere.getRadius() for sphere in particle.getSpheres()),
         max(sphere.getCenter().y + sphere.getRadius() for sphere in particle.getSpheres()) - 
         min(sphere.getCenter().y - sphere.getRadius() for sphere in particle.getSpheres()),
         max(sphere.getCenter().z + sphere.getRadius() for sphere in particle.getSpheres()) - 
         min(sphere.getCenter().z - sphere.getRadius() for sphere in particle.getSpheres())]
    )
    
    mid_x = (max(sphere.getCenter().x for sphere in particle.getSpheres()) + 
             min(sphere.getCenter().x for sphere in particle.getSpheres())) / 2
    mid_y = (max(sphere.getCenter().y for sphere in particle.getSpheres()) + 
             min(sphere.getCenter().y for sphere in particle.getSpheres())) / 2
    mid_z = (max(sphere.getCenter().z for sphere in particle.getSpheres()) + 
             min(sphere.getCenter().z for sphere in particle.getSpheres())) / 2

    ax.set_xlim(mid_x - max_range/2, mid_x + max_range/2)
    ax.set_ylim(mid_y - max_range/2, mid_y + max_range/2)
    ax.set_zlim(mid_z - max_range/2, mid_z + max_range/2)
    
    return ax

def main():
    """
    Main function to demonstrate the usage of the packing generator.
    """
    print("Random Packing Generator - Python Demo")
    print("======================================")
    
    # Create a packing
    generator = create_packing(
        size=300,
        core_radius_range=(30, 40),
        secondary_radius_range=(20, 30),
        tertiary_radius_range=(5, 10),
        tertiary_volume_fraction=0.1,
        target_density=0.65,
        compactness_factor=0.5,
        output_file="packing.tiff"
    )
    
    if generator is None:
        return
    
    # Analyze the packing
    analyze_packing(generator)
    
    # Visualize a few particles
    if generator.getParticleCount() > 0:
        print("\nVisualizing example particles...")
        fig = plt.figure(figsize=(15, 5))
        
        # Get a few particles with different characteristics
        particles = generator.getParticles()
        
        # Sort particles by sphericity
        particles_by_sphericity = sorted(particles, key=lambda p: p.calculateSphericity())
        
        if len(particles) >= 3:
            # Plot lowest sphericity particle
            ax1 = fig.add_subplot(131, projection='3d')
            visualize_particle(particles_by_sphericity[0], ax1)
            
            # Plot median sphericity particle
            middle_idx = len(particles) // 2
            ax2 = fig.add_subplot(132, projection='3d')
            visualize_particle(particles_by_sphericity[middle_idx], ax2)
            
            # Plot highest sphericity particle
            ax3 = fig.add_subplot(133, projection='3d')
            visualize_particle(particles_by_sphericity[-1], ax3)
            
            plt.tight_layout()
            plt.savefig('example_particles.png')
            print("Saved example particle visualizations to example_particles.png")
        else:
            # Just plot one particle
            ax = fig.add_subplot(111, projection='3d')
            visualize_particle(particles[0], ax)
            plt.tight_layout()
            plt.savefig('example_particle.png')
            print("Saved example particle visualization to example_particle.png")
    
    print("\nDemo completed successfully!")

if __name__ == "__main__":
    main()