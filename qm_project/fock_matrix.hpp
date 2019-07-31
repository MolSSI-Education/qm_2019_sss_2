#pragma once

#include <Eigen/Dense>
#include <vector>

double chi_on_atom(int o1, int o2, int o3, double dipole);

void calculate_fock_matrix_fast(Eigen::MatrixXd hamiltonian_matrix, Eigen::MatrixXd interaction_matrix, Eigen::MatrixXd density_matrix);
