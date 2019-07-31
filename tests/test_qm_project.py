"""
Unit and regression test for the qm_project package.
"""

# Import package, test suite, and other packages as needed
from qm_project.semi_empirical_model import *
import pytest
# import sys

@pytest.fixture
def system(request):
   """ Test the coulomb energy """
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
   return system

@pytest.fixture
def model(request):
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
    return model

def test_Model():
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
    assert True

def test_System(model):
    atomic_coordinates = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 5.0]])
    system = System(atomic_coordinates, model)
    assert True

def test_HF(system,model):
    """Sample test, will always pass so long as import statement worked"""
    energy = HartreeFock(system, model)
    print(energy.system.hamiltonian_matrix)
    print(energy.system.interaction_matrix)
    print(energy.system.density_matrix)
    print(energy.system.chi_tensor)
    print(energy.fock_matrix)
    print(system.energy_ion)
    energy.scf_cycle()
    print(energy.get_hartree_fock_energy())
    assert True

def test_MP2(system,model):
    """Sample test, will always pass so long as import statement worked"""
    energy = MP2(system, model)
    print(energy.system.hamiltonian_matrix)
    print(energy.system.interaction_matrix)
    print(energy.system.density_matrix)
    print(energy.system.chi_tensor)
    print(energy.fock_matrix)
    print(system.energy_ion)
    energy.scf_cycle()
    print(energy.get_hartree_fock_energy())
    assert True

def test_coulomb_energy(system):
       o1 = 's'
       o2 = 's'
       r12 = [0,0,1]
       expected_value = 1.0
       calculated_value = system.coulomb_energy(o1, o2, r12)
       assert expected_value == calculated_value

def test_chi_on_atom(system):
   """ Test for chi_on_atom function"""
   o1 = 'pz'
   o2 = 's'
   o3 = 'pz'
   expected_value = 2.781629275106456
   calculated_value = system.chi_on_atom(o1, o2, o3)
   assert expected_value == calculated_value

def test_hopping_energy(system):
   """ Test for hopping_energy function"""
   o1 = 's'
   o2 = 'px'
   r12 = np.array([3.1810226927827516,0.0,0.0])
   expected_value = -0.029154833035109226
   calculated_value = system.hopping_energy(o1, o2, r12)
   assert expected_value == calculated_value

def test_calculate_potential_vector(system):
   """ Test for calculate_potential_vector function"""
   atomic_coordinates = np.array([ [0.0,0.0,0.0] , [3.0,4.0,5.0] ])
   expected_value = np.array([-0.8, -0.1, -0.1, -0.1, -0.8,  0.1,  0.1,  0.1])
   calculated_value = system.calculate_potential_vector()
   assert expected_value.all() == calculated_value.all()

def test_hamiltonian(system):
    expected_hamil = np.array([[ 3.2e+00,  2.5e-04,  3.3e-04,  4.2e-04,  6.5e-04,  5.3e-04,  7.1e-04,  8.9e-04],
 [ 2.5e-04, -2.4e+00,  0.0e+00,  0.0e+0,  5.3e-04, 2.9e-04,  2.2e-03,  2.7e-03],
 [ 3.3e-04,  0.0e+00, -2.4e+00,  0.0e+00,  7.1e-04,  2.2e-03,  1.6e-03,  3.6e-03],
 [ 4.2e-04,  0.0e+00,  0.0e+00, -2.4e+00,  8.9e-04,  2.7e-03,  3.6e-03,  3.2e-03],
 [ 6.5e-04, -5.3e-04, -7.1e-04, -8.9e-04,  3.2e+00, -2.5e-04, -3.3e-04, -4.2e-04],
 [-5.3e-04,  2.9e-04,  2.2e-03,  2.7e-03, -2.5e-04, -2.4e+00,  0.0e+00,  0.0e+00],
 [-7.1e-04,  2.2e-03,  1.6e-03,  3.6e-03, -3.3e-04,  0.0e+00, -2.4e+00,  0.0e+00],
 [-8.9e-04,  2.7e-03,  3.6e-03,  3.2e-03, -4.2e-04,  0.0e+00,  0.0e+00, -2.4e+00]])
    calculated_hamil = system.hamiltonian_matrix
    assert calculated_hamil.all() == expected_hamil.all()

def test_interaction(system):
    expected_interaction = np.array([[ 3.6e-01,  0.0e+00,  0.0e+00,  0.0e+00,  7.4e-23, -4.4e-24, -5.9e-24, -7.4e-24],
 [ 0.0e+00, -3.3e-03,  0.0e+00,  0.0e+00, -4.4e-24,  6.8e-25, -1.1e-24, -1.3e-24],
 [ 0.0e+00,  0.0e+00, -3.3e-03,  0.0e+00, -5.9e-24, -1.1e-24,  5.9e-26, -1.8e-24],
 [ 0.0e+00,  0.0e+00,  0.0e+00, -3.3e-03, -7.4e-24, -1.3e-24, -1.8e-24, -7.4e-25],
 [ 7.4e-23,  4.4e-24,  5.9e-24,  7.4e-24,  3.6e-01,  0.0e+00,  0.0e+00,  0.0e+00],
 [ 4.4e-24,  6.8e-25, -1.1e-24, -1.3e-24,  0.0e+00, -3.3e-03,  0.0e+00,  0.0e+00],
 [ 5.9e-24, -1.1e-24,  5.9e-26, -1.8e-24,  0.0e+00,  0.0e+00, -3.3e-03,  0.0e+00],
 [ 7.4e-24, -1.3e-24, -1.8e-24, -7.4e-25,  0.0e+00,  0.0e+00,  0.0e+00, -3.3e-03]])
    calc_interaction = system.interaction_matrix
    assert calc_interaction.all() == expected_interaction.all()

def test_chi_tensor(system):
    expected_tensor = np.array([[
  [1.,  0.,  0.,  0.,  0.,  0.,  0.,  0., ],
  [0.,  2.8, 0.,  0.,  0.,  0.,  0.,  0., ],
  [0.,  0.,  2.8, 0.,  0.,  0.,  0.,  0., ],
  [0.,  0.,  0.,  2.8, 0.,  0.,  0.,  0., ],
  [0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., ],
  [0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., ],
  [0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., ],
  [0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., ]],

 [[0.  ,2.8 ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [1.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ]],

 [[0.  ,0.  ,2.8 ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [1.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ]],

 [[0.  ,0.  ,0.  ,2.8 ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [1.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ]],

 [[0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,1.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,2.8 ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,2.8 ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,2.8]],

 [[0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,2.8 ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,1.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ]],

 [[0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,2.8 ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,1.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ]],

 [[0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,2.8],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0. ],
  [0.  ,0.  ,0.  ,0.  ,1.  ,0.  ,0.  ,0. ]]])
    calc_chi_tensor = system.chi_tensor
    assert expected_tensor.all() == calc_chi_tensor.all()

def test_atomic_density_matrix(system):
    calc_density = np.array([
 [0. ,0. ,0. ,0. ,0. ,0. ,0. ,0.],
 [0. ,1. ,0. ,0. ,0. ,0. ,0. ,0.],
 [0. ,0. ,1. ,0. ,0. ,0. ,0. ,0.],
 [0. ,0. ,0. ,1. ,0. ,0. ,0. ,0.],
 [0. ,0. ,0. ,0. ,0. ,0. ,0. ,0.],
 [0. ,0. ,0. ,0. ,0. ,1. ,0. ,0.],
 [0. ,0. ,0. ,0. ,0. ,0. ,1. ,0.],
 [0. ,0. ,0. ,0. ,0. ,0. ,0. ,1.]])
    expected_density = system.density_matrix
    assert calc_density.all() == expected_density.all()

def test_fock_matrix_initial(system,model):
    expected_fock = np.array([
 [ 5.4e+00,  2.5e-04,  3.3e-04,  4.2e-04,  6.5e-04,  5.3e-04,  7.1e-04,  8.9e-04],
 [ 2.5e-04, -5.9e-01,  0.0e+00,  0.0e+00,  5.3e-04,  2.9e-04,  2.2e-03,  2.7e-03],
 [ 3.3e-04,  0.0e+00, -5.9e-01,  0.0e+00,  7.1e-04,  2.2e-03,  1.6e-03,  3.6e-03],
 [ 4.2e-04,  0.0e+00,  0.0e+00, -5.9e-01,  8.9e-04,  2.7e-03,  3.6e-03,  3.2e-03],
 [ 6.5e-04, -5.3e-04, -7.1e-04, -8.9e-04,  5.4e+00, -2.5e-04, -3.3e-04, -4.2e-04],
 [-5.3e-04,  2.9e-04,  2.2e-03,  2.7e-03, -2.5e-04, -5.9e-01,  0.0e+00,  0.0e+00],
 [-7.1e-04,  2.2e-03,  1.6e-03,  3.6e-03, -3.3e-04,  0.0e+00, -5.9e-01,  0.0e+00],
 [-8.9e-04,  2.7e-03,  3.6e-03,  3.2e-03, -4.2e-04,  0.0e+00,  0.0e+00, -5.9e-01]])
    energy = HartreeFock(system, model)
    calc_fock = energy.fock_matrix
    assert calc_fock.all() == expected_fock.all()

def test_fock_matrix_after_scf(system,model,capsys):
    energy = HartreeFock(system,model)
    energy.scf_cycle()
    calc_fock = energy.fock_matrix
    expected_fock = np.array([
 [ 5.4e+00,  2.7e-04,  3.6e-04,  4.5e-04,  6.5e-04,  5.3e-04,  7.1e-04,  8.9e-04],
 [ 2.7e-04, -5.9e-01,  4.8e-09,  6.0e-09,  5.3e-04,  2.9e-04,  2.2e-03,  2.7e-03],
 [ 3.6e-04,  4.8e-09, -5.9e-01,  8.0e-09,  7.1e-04,  2.2e-03,  1.6e-03,  3.6e-03],
 [ 4.5e-04,  6.0e-09,  8.0e-09, -5.9e-01,  8.9e-04,  2.7e-03,  3.6e-03,  3.2e-03],
 [ 6.5e-04, -5.3e-04, -7.1e-04, -8.9e-04,  5.4e+00, -2.7e-04, -3.6e-04, -4.5e-04],
 [-5.3e-04,  2.9e-04,  2.2e-03,  2.7e-03, -2.7e-04, -5.9e-01,  4.8e-09,  6.0e-09],
 [-7.1e-04,  2.2e-03,  1.6e-03,  3.6e-03, -3.6e-04,  4.8e-09, -5.9e-01,  8.0e-09],
 [-8.9e-04,  2.7e-03,  3.6e-03,  3.2e-03, -4.5e-04,  6.0e-09,  8.0e-09, -5.9e-01]])
    captured = capsys.readouterr()
    print(calc_fock, "\n\n", expected_fock )
    # captured = capsys.readouterr()
    # print(expected_fock)
    assert calc_fock.all() == expected_fock.all()

def test_calculate_energy_scf(system, model):
    """ Test for calcu function"""
    expected_value = -17.9011807466738
    energy = HartreeFock(system, model)
    energy.scf_cycle()
    calculated_value = energy.calculate_energy_scf()
    assert expected_value == calculated_value
