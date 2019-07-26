"""
Unit and regression test for the qm_project package.
"""

# Import package, test suite, and other packages as needed
from qm_project.semi_empirical_model import *
import pytest
import sys

def test_HF():
    """Sample test, will always pass so long as import statement worked"""
    # --------------------
    # Noble Gas Parameters
    # --------------------
    ionic_charge = 6
    orbital_types = ['s', 'px', 'py', 'pz']
    vec = {'px': [1, 0, 0], 'py': [0, 1, 0], 'pz': [0, 0, 1]}
    orbital_occupation = {'s': 0, 'px': 1, 'py': 1, 'pz': 1}
    model_parameters = {
        'r_hop': 3.1810226927827516,
        't_ss': 0.03365982238611262,
        't_sp': -0.029154833035109226,
        't_pp1': -0.0804163845390335,
        't_pp2': -0.01393611496959445,
        'r_pseudo': 2.60342991362958,
        'v_pseudo': 0.022972992186364977,
        'dipole': 2.781629275106456,
        'energy_s': 3.1659446174413004,
        'energy_p': -2.3926873325346554,
        'coulomb_s': 0.3603533286088998,
        'coulomb_p': -0.003267991835806299
    }
    # First Instantiation of Model class.
    model = Model(model_parameters, ionic_charge, orbital_types, orbital_occupation, vec)
    atomic_coordinates = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 5.0]])
    system = System(atomic_coordinates, model)
    energy = MP2(system, model)
    # print(energy.get_scf_energy())
    print(system.energy_ion)
    energy.run_scf()
    # print(system.density_matrix)
    # print(energy.fock_matrix)
    # print(system.hamiltonian_matrix)
    # print(system.chi_tensor) ## This prints correct.
    print(energy.get_hartree_fock_energy())
    # print(energy.get_partition_orbitals())
    #
    # # MP 2 - Fock matrix from SCF
    # occupied_energy, virtual_energy, occupied_matrix, virtual_matrix = partition_orbitals(fock_matrix)
    # interaction_tensor = transform_interaction_tensor(occupied_matrix, virtual_matrix, interaction_matrix, chi_tensor)
    # energy_mp2 = calculate_energy_mp2(fock_matrix, interaction_matrix, chi_tensor)
    # print(energy_mp2)
    #
    #
    assert True


def test_MP2():
    """Sample test, will always pass so long as import statement worked"""
    # --------------------
    # Noble Gas Parameters
    # --------------------
    ionic_charge = 6
    orbital_types = ['s', 'px', 'py', 'pz']
    vec = {'px': [1, 0, 0], 'py': [0, 1, 0], 'pz': [0, 0, 1]}
    orbital_occupation = {'s': 0, 'px': 1, 'py': 1, 'pz': 1}
    model_parameters = {
        'r_hop': 3.1810226927827516,
        't_ss': 0.03365982238611262,
        't_sp': -0.029154833035109226,
        't_pp1': -0.0804163845390335,
        't_pp2': -0.01393611496959445,
        'r_pseudo': 2.60342991362958,
        'v_pseudo': 0.022972992186364977,
        'dipole': 2.781629275106456,
        'energy_s': 3.1659446174413004,
        'energy_p': -2.3926873325346554,
        'coulomb_s': 0.3603533286088998,
        'coulomb_p': -0.003267991835806299
    }
    # First Instantiation of Model class.
    model = Model(model_parameters, ionic_charge, orbital_types, orbital_occupation, vec)
    atomic_coordinates = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 5.0]])
    system = System(atomic_coordinates, model)
    energy = MP2(system, model)
    # print(energy.get_scf_energy())
    print(system.energy_ion)
    energy.run_scf()
    # print(system.density_matrix)
    # print(energy.fock_matrix)
    # print(system.hamiltonian_matrix)
    # print(system.chi_tensor) ## This prints correct.
    print(energy.get_hartree_fock_energy())
    # print(energy.get_partition_orbitals())
    #
    # # MP 2 - Fock matrix from SCF
    # occupied_energy, virtual_energy, occupied_matrix, virtual_matrix = partition_orbitals(fock_matrix)
    # interaction_tensor = transform_interaction_tensor(occupied_matrix, virtual_matrix, interaction_matrix, chi_tensor)
    # energy_mp2 = calculate_energy_mp2(fock_matrix, interaction_matrix, chi_tensor)
    # print(energy_mp2)
    #
    #
    assert True

