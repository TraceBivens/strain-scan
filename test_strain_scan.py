#!/usr/bin/env python3
"""
Unit tests for strain_scan.py
"""

import unittest
import tempfile
import os
from pathlib import Path
import numpy as np

from strain_scan import parse_poscar, apply_strain, write_poscar


class TestStrainScan(unittest.TestCase):

    def setUp(self):
        # Simple cubic BN POSCAR content
        self.poscar_content = """Cubic BN
3.57
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0
B N
1 1
Direct
0.00 0.00 0.00
0.25 0.25 0.25
"""

    def test_parse_poscar(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.poscar', delete=False) as f:
            f.write(self.poscar_content)
            f.flush()
            structure = parse_poscar(f.name)

        self.assertEqual(structure['comment'], 'Cubic BN')
        self.assertEqual(structure['scaling'], 3.57)
        self.assertEqual(structure['species'], ['B', 'N'])
        self.assertEqual(structure['counts'], [1, 1])
        self.assertFalse(structure['selective_dynamics'])
        self.assertEqual(structure['coordinate_type'], 'Direct')
        np.testing.assert_array_equal(structure['positions'], [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]])

    def test_apply_strain_volumetric(self):
        structure = {
            'comment': 'Test',
            'scaling': 1.0,
            'lattice': np.eye(3),
            'species': None,
            'counts': [1],
            'selective_dynamics': False,
            'coordinate_type': 'Direct',
            'positions': np.array([[0.0, 0.0, 0.0]]),
            'dynamics_flags': None,
        }

        strained = apply_strain(structure, 0.01, 'volumetric')
        expected_scale = (1 + 0.01)**(1/3)
        np.testing.assert_array_almost_equal(strained['lattice'], np.eye(3) * expected_scale)
        self.assertIn('Strain: 0.010000 (volumetric)', strained['comment'])

    def test_apply_strain_uniaxial_x(self):
        structure = {
            'comment': 'Test',
            'scaling': 1.0,
            'lattice': np.eye(3),
            'species': None,
            'counts': [1],
            'selective_dynamics': False,
            'coordinate_type': 'Cartesian',
            'positions': np.array([[1.0, 0.0, 0.0]]),
            'dynamics_flags': None,
        }

        strained = apply_strain(structure, 0.01, 'uniaxial-x')
        expected_lattice = np.eye(3)
        expected_lattice[0] *= 1.01
        expected_positions = np.array([[1.01, 0.0, 0.0]])
        np.testing.assert_array_almost_equal(strained['lattice'], expected_lattice)
        np.testing.assert_array_almost_equal(strained['positions'], expected_positions)

    def test_write_poscar(self):
        structure = {
            'comment': 'Test',
            'scaling': 1.0,
            'lattice': np.eye(3),
            'species': ['Si'],
            'counts': [1],
            'selective_dynamics': False,
            'coordinate_type': 'Direct',
            'positions': np.array([[0.0, 0.0, 0.0]]),
            'dynamics_flags': None,
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.poscar', delete=False) as f:
            write_poscar(structure, f.name)
            with open(f.name, 'r') as rf:
                content = rf.read()

        self.assertIn('Test', content)
        self.assertIn('1.000000', content)
        self.assertIn('Direct', content)


if __name__ == '__main__':
    unittest.main()