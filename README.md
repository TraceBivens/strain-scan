# strain-scan

A CLI tool to generate strained POSCAR files from VASP POSCAR input at evenly spaced intervals of strain.

## Installation

The tool is packaged as a conda package using Pixi's build backends.

1. Ensure Pixi is installed.

2. Clone the repository and build the conda package:
```bash
git clone <repo>
cd strain-scan
pixi build
```

3. Install globally with Pixi:
```bash
pixi global install --path ./strain-scan-0.1.0-pyh4616a5c_0.conda
```

This installs the `strain-scan` command globally.

Requires Python 3.8+ and numpy.

## Usage

```bash
strain-scan POSCAR [options]
```

### Options

- `POSCAR`: Input POSCAR file (required)
- `--min-strain FLOAT`: Minimum strain value (default: -0.05)
- `--max-strain FLOAT`: Maximum strain value (default: 0.05)
- `--steps INT`: Number of strain steps (default: 11)
- `--strain-type {volumetric,uniaxial-x,uniaxial-y,uniaxial-z,biaxial}`: Type of strain (default: volumetric)
- `--output-dir PATH`: Output directory (default: strained)
- `--prefix STR`: Filename prefix (default: POSCAR_strain_)

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

## Output

Files are named as `POSCAR_strain_001.poscar`, `POSCAR_strain_002.poscar`, etc., in the specified output directory.

Each file includes the strain value in the comment line.

## Requirements

- Python 3.6+
- numpy

## License

MIT