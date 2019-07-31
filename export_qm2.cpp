#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "fock_matrix.hpp"

PYBIND11_MODULE(sss_cpp, m)
{
  m.def("chi_on_atom", chi_on_atom, "calculate chi on atom");
  m.def("calculate_fock_matrix_fast", calculate_fock_matrix_fast, "calculate fock matrix fast");
}
