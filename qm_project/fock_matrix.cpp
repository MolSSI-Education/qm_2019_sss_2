#include <iostream>
#include "fock_matrix.hpp"
#include <vector>
#include <Eigen/Dense>

double dipole = 2.781629275106456;

double chi_on_atom(int o1, int o2, int o3, double dipole)
{
  if (o1 == o2 && o3 == 0)
      return 1.0;
  else if (o1 == o3 && o1!=0 && o2 == 0)
      return dipole;
  else if (o2 == o3 && o2!=0 && o1 == 0)
      return dipole;
  else
      return 0;
 }




void calculate_fock_matrix_fast(int ndof, int orbitals_per_atom )
{
/*def calculate_fock_matrix_fast(hamiltonian_matrix, interaction_matrix, density_matrix, model_parameters):
    '''Returns the Fock matrix defined by the input Hamiltonian, interaction, & density matrices.'''
*/

    int orbitals_per_atom = 4;
    int ndof = hamiltonian_matrix.rows()
    double fock_matrix = hamiltonian_matrix()



    // Hartree potential term


    //for p in range(ndof):
    for (int p=0; p<ndof; p++)
    {
        //for orb_q in orbital_types:
        for (int orb_q=0; orb_q<orbitals_per_atom; orb_q++)
        {
            q = ao_index(atom(p), orb_q) //p & q on same atom
            //for orb_t in orbital_types:
            for (int orb_t=0; orb_t<orbitals_per_atom; orb_t++)
            {
                t = ao_index(atom(p), orb_t) // p & t on same atom
                chi_pqt = chi_on_atom(orb(p), orb_q, orb_t, dipole)
                //for r in range(ndof):
                for (int r=0; r<ndof; r++)
                {
                    //for orb_s in orbital_types:
                    for (int orb_s=0; orb_s<orbitals_per_atom; orb_s++)
                    {
                        s = ao_index(atom(r), orb_s) // r & s on same atom
                        //for orb_u in orbital_types:
                        for (int orb_u=0; orb_u<orbitals_per_atom; orb_u++)
                        {
                            u = ao_index(atom(r), orb_u) // r & u on same atom
                            chi_rsu = chi_on_atom(orb(r), orb_s, orb_u, dipole)
                            fock_matrix[p,q] += 2.0 * chi_pqt * chi_rsu * interaction_matrix[t,u] * density_matrix[r,s]
                        }
                    }
                }
            }
         }
      }
    //Fock exchange term
    //for p in range(ndof):
    for (int p=0; p<ndof; p++)
    {
        //for orb_s in orbital_types:
        for (int orb_s=0; orb_s<orbitals_per_atom; orb_s++)
        {
            s = ao_index(atom(p), orb_s) // p & s on same atom
            //for orb_u in orbital_types:
            for (int orb_u=0; orb_u<orbitals_per_atom; orb_u++)
            {
                u = ao_index(atom(p), orb_u) // p & u on same atom
                chi_psu = chi_on_atom(orb(p), orb_s, orb_u, dipole)
                //for q in range(ndof):
                for (int q=0; q<ndof; q++)
                {
                    //for orb_r in orbital_types:
                    for (int orb_r=0; orb_r<orbitals_per_atom; orb_r++)
                    {
                        r = ao_index(atom(q), orb_r) // q & r on same atom
                        //for orb_t in orbital_types:
                        for (int orb_t=0; orb_t<orbitals_per_atom; orb_t++)
                        {
                            t = ao_index(atom(q), orb_t) // q & t on same atom
                            chi_rqt = chi_on_atom(orb_r, orb(q), orb_t, dipole)
                            fock_matrix[p,q] -= chi_rqt * chi_psu * interaction_matrix[t,u] * density_matrix[r,s]
                        }
                    }
                 }
              }
          }
      }
    return fock_matrix
}
