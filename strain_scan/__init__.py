#!/usr/bin/env python3
"""
strain-scan: CLI tool to generate strained POSCAR files from VASP POSCAR input.

Generates new POSCAR files at evenly spaced intervals of strain.
"""

import argparse
import sys
from pathlib import Path
import numpy as np
from .plotting import plot_energy_vs_strain

__all__ = ['parse_poscar', 'apply_strain', 'write_poscar', 'plot_energy_vs_strain', 'main']


def parse_poscar(filepath):
    """
    Parse a VASP POSCAR file and return a dictionary with structure data.

    Returns:
        dict: {
            'comment': str,
            'scaling': float or list,
            'lattice': np.array(3,3),
            'species': list[str],
            'counts': list[int],
            'selective_dynamics': bool,
            'coordinate_type': str,  # 'Direct' or 'Cartesian'
            'positions': np.array(n,3),
            'dynamics_flags': np.array(n,3) or None,  # if selective dynamics
        }
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Remove trailing whitespace
    lines = [line.strip() for line in lines]

    # Comment
    comment = lines[0]

    # Scaling factor
    scaling_line = lines[1].split()
    if len(scaling_line) == 1:
        scaling = float(scaling_line[0])
    else:
        scaling = [float(x) for x in scaling_line]

    # Lattice vectors
    lattice = np.array([
        [float(x) for x in lines[2].split()],
        [float(x) for x in lines[3].split()],
        [float(x) for x in lines[4].split()]
    ])

    # Species names (optional)
    line_idx = 5
    species_line = lines[line_idx].split()
    # Check if it's species or counts
    try:
        counts = [int(x) for x in species_line]
        species = None
        line_idx += 1
    except ValueError:
        species = species_line
        line_idx += 1
        counts = [int(x) for x in lines[line_idx].split()]
        line_idx += 1

    # Selective dynamics
    selective_dynamics = False
    if lines[line_idx].upper().startswith('S'):
        selective_dynamics = True
        line_idx += 1

    # Coordinate type
    coord_type = lines[line_idx]
    line_idx += 1

    # Positions
    total_atoms = sum(counts)
    positions = []
    dynamics_flags = []

    for i in range(total_atoms):
        parts = lines[line_idx + i].split()
        pos = [float(x) for x in parts[:3]]
        positions.append(pos)
        if selective_dynamics:
            flags = parts[3:6]
            dynamics_flags.append(flags)

    positions = np.array(positions)
    if not selective_dynamics:
        dynamics_flags = None
    else:
        dynamics_flags = np.array(dynamics_flags)

    # Basic validation
    if lattice.shape != (3, 3):
        raise ValueError("Lattice must be 3x3")

    if len(positions) != sum(counts):
        raise ValueError("Number of positions does not match total atoms")

    if coord_type.lower() not in ['direct', 'cartesian']:
        raise ValueError("Coordinate type must be 'Direct' or 'Cartesian'")

    return {
        'comment': comment,
        'scaling': scaling,
        'lattice': lattice,
        'species': species,
        'counts': counts,
        'selective_dynamics': selective_dynamics,
        'coordinate_type': coord_type,
        'positions': positions,
        'dynamics_flags': dynamics_flags,
    }


def apply_strain(structure, strain, strain_type):
    """
    Apply strain to the structure.

    Args:
        structure: dict from parse_poscar
        strain: float, strain value
        strain_type: str, 'volumetric', 'uniaxial-x', etc.

    Returns:
        dict: modified structure
    """
    new_structure = structure.copy()
    lattice = structure['lattice'].copy()
    positions = structure['positions'].copy()
    coord_type = structure['coordinate_type'].lower()

    if strain_type == 'volumetric':
        scale_factor = (1 + strain)**(1/3)
        lattice *= scale_factor
        if coord_type == 'cartesian':
            positions *= scale_factor

    elif strain_type.startswith('uniaxial-'):
        axis = {'x': 0, 'y': 1, 'z': 2}[strain_type[-1]]
        scale_factor = 1 + strain
        lattice[axis] *= scale_factor
        if coord_type == 'cartesian':
            positions[:, axis] *= scale_factor

    elif strain_type == 'biaxial':
        scale_factor = 1 + strain
        lattice[0] *= scale_factor  # a vector
        lattice[1] *= scale_factor  # b vector
        if coord_type == 'cartesian':
            positions[:, 0] *= scale_factor
            positions[:, 1] *= scale_factor

    new_structure['lattice'] = lattice
    new_structure['positions'] = positions
    # Update comment to include strain
    new_structure['comment'] = f"{structure['comment']} | Strain: {strain:.6f} ({strain_type})"

    return new_structure


def write_poscar(structure, filepath):
    """
    Write structure to POSCAR file.
    """
    with open(filepath, 'w') as f:
        # Comment
        f.write(f"{structure['comment']}\n")

        # Scaling
        scaling = structure['scaling']
        if isinstance(scaling, list):
            f.write(" ".join(f"{x:.6f}" for x in scaling) + "\n")
        else:
            f.write(f"{scaling:.6f}\n")

        # Lattice
        for vec in structure['lattice']:
            f.write(" ".join(f"{x:.9f}" for x in vec) + "\n")

        # Species
        if structure['species']:
            f.write(" ".join(structure['species']) + "\n")
        f.write(" ".join(str(c) for c in structure['counts']) + "\n")

        # Selective dynamics
        if structure['selective_dynamics']:
            f.write("Selective dynamics\n")

        # Coordinate type
        f.write(f"{structure['coordinate_type']}\n")

        # Positions
        for i, pos in enumerate(structure['positions']):
            line = " ".join(f"{x:.9f}" for x in pos)
            if structure['dynamics_flags'] is not None:
                flags = structure['dynamics_flags'][i]
                line += " " + " ".join(flags)
            f.write(line + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate strained POSCAR files at evenly spaced strain intervals",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s POSCAR --min-strain -0.05 --max-strain 0.05 --steps 11
  %(prog)s POSCAR --strain-type biaxial --output-dir strained_2d
  %(prog)s --plot --calc-dir /path/to/vasp/calculations
  %(prog)s --plot --calc-dir calculations --output-dir strained
        """
    )

    parser.add_argument(
        "poscar_file",
        type=Path,
        nargs='?',
        help="Input POSCAR file (required unless --plot)"
    )

    parser.add_argument(
        "--min-strain",
        type=float,
        default=-0.05,
        help="Minimum strain value (default: -0.05)"
    )

    parser.add_argument(
        "--max-strain",
        type=float,
        default=0.05,
        help="Maximum strain value (default: 0.05)"
    )

    parser.add_argument(
        "--steps",
        type=int,
        default=11,
        help="Number of strain steps (default: 11)"
    )

    parser.add_argument(
        "--strain-type",
        choices=["volumetric", "uniaxial-x", "uniaxial-y", "uniaxial-z", "biaxial"],
        default="volumetric",
        help="Type of strain to apply (default: volumetric)"
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("strained"),
        help="Output directory for strained POSCAR files (default: strained)"
    )

    parser.add_argument(
        "--calc-dir",
        type=Path,
        help="Directory containing VASP calculation subdirectories (default: --output-dir)"
    )

    parser.add_argument(
        "--prefix",
        type=str,
        default="POSCAR_strain_",
        help="Prefix for output filenames (default: POSCAR_strain_)"
    )

    parser.add_argument(
        "--plot",
        action="store_true",
        help="Plot final energies from VASP runs in subdirectories"
    )

    parser.add_argument(
        "--plot-file",
        type=Path,
        default="energy_vs_strain.png",
        help="Output file for the energy plot (default: energy_vs_strain.png)"
    )

    args = parser.parse_args()

    # Validate inputs
    if args.min_strain >= args.max_strain:
        parser.error("min-strain must be less than max-strain")

    if args.steps < 2:
        parser.error("steps must be at least 2")

    if not args.plot and not args.poscar_file:
        parser.error("POSCAR file is required unless --plot is specified")

    if args.poscar_file and not args.poscar_file.exists():
        parser.error(f"POSCAR file '{args.poscar_file}' does not exist")

    if args.plot:
        # Plotting mode - use separate plotting function
        calc_dir = args.calc_dir if args.calc_dir else args.output_dir
        result = plot_energy_vs_strain(
            calc_dir=calc_dir,
            min_strain=args.min_strain,
            max_strain=args.max_strain,
            steps=args.steps,
            prefix=args.prefix,
            plot_file=args.plot_file
        )
        
        if result is None:
            print("No valid data found for plotting")
            sys.exit(1)
    else:
        # Generation mode
        # Create output directory
        args.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Parse input POSCAR
        structure = parse_poscar(args.poscar_file)

        # Generate strain values
        strains = np.linspace(args.min_strain, args.max_strain, args.steps)

        # Process each strain
        for i, strain in enumerate(strains):
            strained_structure = apply_strain(structure, strain, args.strain_type)
            output_file = args.output_dir / f"{args.prefix}{i+1:03d}.poscar"
            write_poscar(strained_structure, output_file)
            print(f"Generated {output_file} (strain: {strain:.6f})")

        print(f"\nCompleted: {args.steps} strained POSCAR files generated in {args.output_dir}")


if __name__ == "__main__":
    main()