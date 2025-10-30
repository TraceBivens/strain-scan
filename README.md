# strain-scan

A CLI tool to generate strained POSCAR files from VASP POSCAR input at evenly spaced intervals of strain.

## Installation

The tool is packaged as a conda package using Pixi's build backends.

1. Ensure Pixi is installed.

2. Clone the repository and build the conda package:
```bash
git clone https://github.com/TraceBivens/strain-scan.git
cd strain-scan
pixi build
```

3. Install globally with Pixi:
```bash
pixi global install --path ./strain-scan-0.1.0-pyh4616a5c_0.conda
```

This installs the `strain-scan` command globally.

Requires Python 3.8+, numpy, py4vasp, and matplotlib.

## Usage

```bash
strain-scan [POSCAR] [options]
```

### Options

- `POSCAR`: Input POSCAR file (required unless --plot)
- `--min-strain FLOAT`: Minimum strain value (default: -0.05)
- `--max-strain FLOAT`: Maximum strain value (default: 0.05)
- `--steps INT`: Number of strain steps (default: 11)
- `--strain-type {volumetric,uniaxial-x,uniaxial-y,uniaxial-z,biaxial}`: Type of strain (default: volumetric)
- `--output-dir PATH`: Output directory (default: strained)
- `--prefix STR`: Filename prefix (default: POSCAR_strain_)
- `--plot`: Plot final energies from VASP runs in subdirectories
- `--plot-file PATH`: Output file for the energy plot (default: energy_vs_strain.png)

### Strain Types

- `volumetric`: Isotropic volume change
- `uniaxial-x/y/z`: Strain along specific axis
- `biaxial`: In-plane strain (x and y), suitable for 2D materials

### Examples

Generate 11 POSCAR files with volumetric strain from -5% to +5%:

```bash
strain-scan POSCAR
```

Biaxial strain for 2D material:

```bash
strain-scan POSCAR --strain-type biaxial --output-dir strained_2d
```

Uniaxial strain along z:

```bash
strain-scan POSCAR --strain-type uniaxial-z --min-strain 0.0 --max-strain 0.1 --steps 6
```

Plot energy vs strain from VASP calculations:

```bash
strain-scan --plot --output-dir vasp_runs
```

This reads `vasprun.xml` files in subdirectories (named with the prefix, e.g., POSCAR_strain_001, etc.), generates `energy_vs_strain.png`, and performs a quadratic fit to find the optimal strain with minimum energy.

## Output

Files are named as `POSCAR_strain_001.poscar`, `POSCAR_strain_002.poscar`, etc., in the specified output directory.

Each file includes the strain value in the comment line.

## Requirements

- Python 3.8+
- numpy
- py4vasp
- matplotlib
- scipy

## License

MIT