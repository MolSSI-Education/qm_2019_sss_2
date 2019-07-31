<<<<<<< HEAD:fock_matrix.hpp
#pragma once

#include <Eigen/Dense>
#include <vector>

double chi_on_atom(int o1, int o2, int o3, double dipole);

void calculate_fock_matrix_fast(Eigen::MatrixXd hamiltonian_matrix, Eigen::MatrixXd interaction_matrix, Eigen::MatrixXd density_matrix);
=======
#pragma once

#include <Eigen/Dense>
#include <vector>

double chi_on_atom(int o1, int o2, int o3, double dipole);

void calculate_fock_matrix_fast(int ndof, int orbitals_per_atom);
>>>>>>> 5c810edf8427ce24db386a4010bd9a3db45c29f7:qm_project/fock_matrix.hpp
