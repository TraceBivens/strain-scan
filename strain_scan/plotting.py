#!/usr/bin/env python3
"""
Plotting functionality for strain-scan energy analysis.

This module provides functions to plot energy vs strain from VASP calculations,
independent of the POSCAR generation functionality.
"""

import numpy as np
import matplotlib.pyplot as plt
import py4vasp as p4v
from pathlib import Path


def plot_energy_vs_strain(
    calc_dir,
    min_strain=-0.05,
    max_strain=0.05,
    steps=11,
    prefix="POSCAR_strain_",
    plot_file="energy_vs_strain.png"
):
    """
    Plot energy vs strain from VASP calculations in subdirectories.
    
    Args:
        calc_dir (Path): Directory containing VASP calculation subdirectories
        min_strain (float): Minimum strain value (default: -0.05)
        max_strain (float): Maximum strain value (default: 0.05)
        steps (int): Number of strain steps (default: 11)
        prefix (str): Prefix for subdirectory names (default: "POSCAR_strain_")
        plot_file (str or Path): Output file for the plot (default: "energy_vs_strain.png")
    
    Returns:
        tuple: (valid_strains, energies, optimal_strain, optimal_energy) or None if no data
    """
    calc_dir = Path(calc_dir)
    plot_file = Path(plot_file)
    
    # Generate strain values
    strains = np.linspace(min_strain, max_strain, steps)
    energies = []
    valid_strains = []
    
    print(f"Looking for VASP results in: {calc_dir}")
    
    # Read energies from vasprun.xml files
    for i, strain in enumerate(strains):
        dir_name = f"{prefix}{i+1:03d}"
        vasprun_path = calc_dir / dir_name / "vasprun.xml"
        
        if vasprun_path.exists():
            try:
                calc = p4v.Calculation.from_file(str(vasprun_path))
                energy = calc.energy.to_numpy()  # total energy (TOTEN) for final step
                energies.append(energy)
                valid_strains.append(strain)
                print(f"Read energy from {vasprun_path}: {energy:.6f} eV")
            except Exception as e:
                print(f"Error reading {vasprun_path}: {e}")
                print(f"  This may indicate incomplete VASP calculations or corrupted vasprun.xml")
                energies.append(np.nan)
                valid_strains.append(strain)
        else:
            print(f"Warning: {vasprun_path} not found")
            print(f"  Expected VASP calculation results in subdirectory: {calc_dir / dir_name}")
            energies.append(np.nan)
            valid_strains.append(strain)
    
    # Plot if we have any data
    if valid_strains:
        plt.figure(figsize=(8, 6))
        plt.plot(valid_strains, energies, 'o-', label='Potential Energy')
        
        optimal_strain = None
        optimal_energy = None
        
        # Fit quadratic and find minimum
        if len(valid_strains) >= 3:
            try:
                coeffs = np.polyfit(valid_strains, energies, 2)
                poly = np.poly1d(coeffs)
                strain_fit = np.linspace(min(valid_strains), max(valid_strains), 100)
                energy_fit = poly(strain_fit)
                plt.plot(strain_fit, energy_fit, '--', label='Quadratic Fit')
                
                # Find minimum: for ax^2 + bx + c, min at -b/(2a)
                a, b, c = coeffs
                if a > 0:
                    optimal_strain = -b / (2 * a)
                    optimal_energy = poly(optimal_strain)
                    print(f"Optimal strain: {optimal_strain:.6f}")
                    print(f"Optimal energy: {optimal_energy:.6f} eV")
                    plt.plot(optimal_strain, optimal_energy, 'rx', markersize=10, label='Minimum')
                else:
                    print("Warning: Quadratic fit is not convex; no minimum found")
            except Exception as e:
                print(f"Warning: Quadratic fit failed: {e}")
        
        plt.xlabel('Strain')
        plt.ylabel('Energy (eV)')
        plt.title('Energy vs Strain')
        plt.grid(True)
        plt.legend()
        
        # Ensure plot directory exists
        plot_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(plot_file)
        print(f"\nPlot saved to {plot_file}")
        
        return valid_strains, energies, optimal_strain, optimal_energy
    else:
        print("No valid vasprun.xml files found for plotting")
        return None