"""
atomicmodel.py
Refactoring of the QM project from MolSSI Software Summer School 2019

Handles the primary functions
"""
import numpy as np


class Model:

    def __init__(self, model_parameters, ionic_charge, orbital_types, orbital_occupation, vec):
        self.model_parameters = model_parameters
        self.ionic_charge = ionic_charge
        self.orbital_types = orbital_types
        self.orbital_occupation = orbital_occupation
        self.p_orbitals = orbital_types[1:]
        self.orbitals_per_atom = len(orbital_types)
        self.vec = vec

class System:

    def __init__(self, atomic_coordinates, model):
        self.atomic_coordinates = atomic_coordinates
        self.model = model
        self.ndof = len(self.atomic_coordinates) * model.orbitals_per_atom
        self.orbital = []
        for i in range(self.ndof):
            atom = int(np.floor(i / model.orbitals_per_atom))
            orbital_num = i % model.orbitals_per_atom
            self.orbital.append([atom, model.orbital_types[orbital_num]])

    @property
    def energy_ion(self):
        """Returns the ionic contribution to the total energy for an input list of atomic coordinates.

        Parameters
        ----------
        atomic_coordinates: numpy.ndarray
            array of atomic coordinates

        Returns
        -------
        energy_ion: float
            ionic contribution to the total energy
        """

        energy_ion = 0.0
        for i, r_i in enumerate(self.atomic_coordinates):
            for j, r_j in enumerate(self.atomic_coordinates):
                if i < j:
                    energy_ion += (self.model.ionic_charge**2) * self.coulomb_energy('s', 's', r_i - r_j)
        return energy_ion

    def orb(self, index):
        return self.orbital[index][1]

    def atom(self, index):
        return self.orbital[index][0]

    @property
    def interaction_matrix(self):
        """Returns the electron-electron interaction energy matrix for an input list of atomic coordinates.

        Parameters
        ----------
        atomic_coordinates: numpy.ndarray
            array of atomic coordinates

        model_parameters: dict,
            dictionary of parameters/constants for noble gases

        Returns
        -------
        interaction_matrix: numpy.ndarray
            electron-electron interaction energy matrix
        """

        interaction_matrix = np.zeros((self.ndof, self.ndof))
        for p in range(self.ndof):
            for q in range(self.ndof):
                if self.atom(p) != self.atom(q):
                    r_pq = self.atomic_coordinates[self.atom(p)] - self.atomic_coordinates[self.atom(q)]
                    interaction_matrix[p, q] = self.coulomb_energy(self.orb(p), self.orb(q), r_pq)
                if p == q and self.orb(p) == 's':
                    interaction_matrix[p, q] = self.model.model_parameters['coulomb_s']
                if p == q and self.orb(p) in self.model.p_orbitals:
                    interaction_matrix[p, q] = self.model.model_parameters['coulomb_p']
        return interaction_matrix

    def calculate_potential_vector(self):
        """Returns the electron-ion potential energy vector for an input list of atomic coordinates.

        Parameters
        ----------
        atomic_coordinates: numpy.ndarray
            array of atomic coordinates

        model_parameters: dict,
            dictionary of parameters/constants for noble gases

        Returns
        -------
        potential_vector: numpy.ndarray
            electron-ion potential energy vector
        """
        potential_vector = np.zeros(self.ndof)
        for p in range(self.ndof):
            potential_vector[p] = 0.0
            for atom_i, r_i in enumerate(self.atomic_coordinates):
                r_pi = self.atomic_coordinates[self.atom(p)] - r_i
                if atom_i != self.atom(p):
                    potential_vector[p] += (self.pseudopotential_energy(self.orb(p), r_pi) -
                                            self.model.ionic_charge * self.coulomb_energy(self.orb(p), 's', r_pi))
        return potential_vector

    def coulomb_energy(self, o1, o2, r12):
        '''Returns the Coulomb matrix element for a pair of multipoles of type o1 & o2 separated by a vector r12.
        Parameters
        ----------
        o1: str
            type of the orbital 1
        o2: str
            type of the orbital 1
        r12: float
            distance between two orbitals before rescaling

        Returns
        -------
        ans: float
            Gives coulomb energy between orbitals

        '''
        r12_length = np.linalg.norm(r12)
        if o1 == 's' and o2 == 's':
            ans = 1.0 / r12_length
        if o1 == 's' and o2 in self.model.p_orbitals:
            ans = np.dot(self.model.vec[o2], r12) / r12_length**3
        if o2 == 's' and o1 in self.model.p_orbitals:
            ans = -1 * np.dot(self.model.vec[o1], r12) / r12_length**3
        if o1 in self.model.p_orbitals and o2 in self.model.p_orbitals:
            ans = (np.dot(self.model.vec[o1], self.model.vec[o2]) / r12_length**3 -
                   3.0 * np.dot(self.model.vec[o1], r12) * np.dot(self.model.vec[o2], r12) / r12_length**5)
        return ans

    def pseudopotential_energy(self, o, r):
        """Returns the energy of a pseudopotential between a multipole of type o and an atom separated by a vector r.

        Parameters
        ----------
        o: str
            orbital type of 1st orbital

        r: float
            separation between atom and multipole

        model_parameters: dict,
            dictionary of parameters/constants for noble gases

        Returns
        -------
        ans: float
            energy of pseudopotential between a multipole of type o and an atom separated by a vector r
        """
        ans = self.model.model_parameters['v_pseudo']
        r_rescaled = r / self.model.model_parameters['r_pseudo']
        ans *= np.exp(1.0 - np.dot(r_rescaled, r_rescaled))
        if o in self.model.p_orbitals:
            ans *= -2.0 * np.dot(self.model.vec[o], r_rescaled)
        return ans

    def hopping_energy(self, o1, o2, r12):
        '''Returns the hopping matrix element for a pair of orbitals of type o1 & o2 separated by a vector r12.
        Parameters
        ----------
        o1: str
            type of the orbital 1
        o2: str
            type of the orbital 1
        r12: float
            distance between two orbitals before rescaling

        model_parameters: float
            constants and required unit conversions

        Returns
        -------
        ans: float
            Gives hopping energy between orbitals
        '''
        r12_rescaled = r12 / self.model.model_parameters['r_hop']
        r12_length = np.linalg.norm(r12_rescaled)
        ans = np.exp(1.0 - r12_length**2)
        if o1 == 's' and o2 == 's':
            ans *= self.model.model_parameters['t_ss']
        if o1 == 's' and o2 in self.model.p_orbitals:
            ans *= np.dot(self.model.vec[o2], r12_rescaled) * self.model.model_parameters['t_sp']
        if o2 == 's' and o1 in self.model.p_orbitals:
            ans *= -1 * np.dot(self.model.vec[o1], r12_rescaled) * self.model.model_parameters['t_sp']
        if o1 in self.model.p_orbitals and o2 in self.model.p_orbitals:
            ans *= ((r12_length**2) * np.dot(self.model.vec[o1], self.model.vec[o2]) * self.model.model_parameters['t_pp2'] -
                    np.dot(self.model.vec[o1], r12_rescaled) * np.dot(self.model.vec[o2], r12_rescaled) *
                    (self.model.model_parameters['t_pp1'] + self.model.model_parameters['t_pp2']))
        return ans

    @property
    def chi_tensor(self):
        '''
        Returns the chi tensor for an input list of atomic coordinates

        Parameters
        ----------
        atomic_coordinates : np.array
            Set of atomic coordinates size (n,3) where n is the number of particles.
        model_parameters : dictionary
            dictionary of known semi-empirical constants to be used in calculations

        Returns
        -------
        chi_tensor : np.array
            Size (n,n,n) array where n is the number of degrees of freedom in the set of coordinates.
        '''

        chi_tensor = np.zeros((self.ndof, self.ndof, self.ndof))
        for p in range(self.ndof):
            for orb_q in self.model.orbital_types:
                q = p % self.model.orbitals_per_atom + self.model.orbital_types.index(orb_q)
                for orb_r in self.model.orbital_types:
                    r = p % self.model.orbitals_per_atom + self.model.orbital_types.index(orb_r)
                    chi_tensor[p, q, r] = self.chi_on_atom(self.orb(p), self.orb(q), self.orb(r))
        return chi_tensor

    @property
    def hamiltonian_matrix(self):
        '''
        Returns the 1-body Hamiltonian matrix for an input list of atomic coordinates.

        Parameters
        ----------
        atomic_coordinates : np.array
            Set of atomic coordinates size (n,3) where n is the number of particles.
        model_parameters : dictionary
            dictionary of known semi-empirical constants to be used in calculations

        Returns
        -------
        hamiltonian_matrix : np.array
            The Hamiltonian Matrix of size (n,n) where n is the number of degrees of freedom in the system

        '''
        hamiltonian_matrix = np.zeros((self.ndof, self.ndof))
        potential_vector = self.calculate_potential_vector()
        for p in range(self.ndof):
            for q in range(self.ndof):
                if self.atom(p) != self.atom(q):
                    r_pq = self.atomic_coordinates[self.atom(p)] - self.atomic_coordinates[self.atom(q)]
                    hamiltonian_matrix[p, q] = self.hopping_energy(self.orb(p), self.orb(q), r_pq)
                if self.atom(p) == self.atom(q):
                    if p == q and self.orb(p) == 's':
                        hamiltonian_matrix[p, q] += self.model.model_parameters['energy_s']
                    if p == q and self.orb(p) in self.model.p_orbitals:
                        hamiltonian_matrix[p, q] += self.model.model_parameters['energy_p']
                    for orb_r in self.model.orbital_types:
                        r = p % self.model.orbitals_per_atom + self.model.orbital_types.index(orb_r)
                        hamiltonian_matrix[p, q] += (self.chi_on_atom(self.orb(p), self.orb(q), orb_r) *
                                                     potential_vector[r])
        return hamiltonian_matrix

    @property
    def density_matrix(self):
        '''
        Returns a trial 1-electron density matrix for an input list of atomic coordinates.

        Parameters
        ----------
        atomic_coordinates : np.array
            Set of atomic coordinates size (n,3) where n is the number of particles.

        Returns
        -------
        density_matrix : np.array
            The density matrix of size (n,n) where n is the number of degrees of freedom in the system.

        '''
        density_matrix = np.zeros((self.ndof, self.ndof))
        print(range(self.ndof))
        for p in range(self.ndof):
            density_matrix[p, p] = self.model.orbital_occupation[self.orb(p)]
            print(F"({p},{p}): {self.orb(p)}, {self.model.orbital_occupation[self.orb(p)]}")
        return density_matrix

    def chi_on_atom(self, o1, o2, o3):
        """Returns the value of the chi tensor for 3 orbital indices on the same atom.

        Parameters
        ----------
        o1: str
            orbital type of 1st orbital

        o2: str
            orbital type of 2nd orbital

        o3: str
            orbital type of 3rd orbital

        model_parameters: dict
            dictionary of parameters/constants for noble gases

        Returns
        -------
        chi: float
            value of chi tensor
        """

        if o1 == o2 and o3 == 's':
            return 1.0
        if o1 == o3 and o3 in self.model.p_orbitals and o2 == 's':
            return self.model.model_parameters['dipole']
        if o2 == o3 and o3 in self.model.p_orbitals and o1 == 's':
            return self.model.model_parameters['dipole']
        return 0.0

class HartreeFock:
    def __init__(self, system, model):
        self.system = system
        self.model = model
        self.max_scf_iterations = 100
        self.convergence_tolerance = 1e-4
        self.mixing_fraction = 0.25
        self._fock_matrix = self.calculate_fock_matrix(system.hamiltonian_matrix, system.interaction_matrix, system.density_matrix, system.chi_tensor)
        self.density_matrix = system.density_matrix
    def calculate_fock_matrix(self, hamiltonian_matrix, interaction_matrix, density_matrix, chi_tensor):

        '''
        Returns the Fock matrix defined by the input Hamiltonian, interaction, & density matrices.


        Parameters
        ----------
        hamiltonian_matrix : np.array
            Hamiltonian Matrix of size (n,n) where n is the number of degrees of freedom in the system,
            produced by calculate_hamiltonian_matrix()
        interaction_matrix : np.array
            Interaction Matrix of size (n,n) where n is the number of degrees of freedom in the system,
            produced by calculate_interaction_matrix()
        density_matrix : np.array
            Density Matrix of size (n,n) where n is the number of degrees of freedom in the system,
            produced by calculate_atomic_density_matrix()
        chi_tensor : np.array
            Chi Tensor produced by calculate_chi_tensor(). Size (n,n,n) array where n is the number
            of degrees of freedom in the set of coordinates.

        Returns
        -------
        fock_matrix : np.array
            Fock Matrix of the same size as the input Hamiltonian Matrix.
        '''
        fock_matrix = hamiltonian_matrix.copy()
        fock_matrix += 2.0 * np.einsum('pqt,rsu,tu,rs',
                                       chi_tensor,
                                       chi_tensor,
                                       interaction_matrix,
                                       density_matrix,
                                       optimize=True)
        fock_matrix -= np.einsum('rqt,psu,tu,rs',
                                 chi_tensor,
                                 chi_tensor,
                                 interaction_matrix,
                                 density_matrix,
                                 optimize=True)
        return fock_matrix



    @property
    def fock_matrix(self):
        self._fock_matrix = self.system.hamiltonian_matrix.copy()
        self._fock_matrix += 2.0 * np.einsum(
            'pqt,rsu,tu,rs', self.system.chi_tensor, self.system.chi_tensor, self.system.interaction_matrix, self.system.density_matrix, optimize=True)
        self._fock_matrix -= np.einsum('rqt,psu,tu,rs',
                                self.system.chi_tensor,
                                self.system.chi_tensor,
                                self.system.interaction_matrix,
                                self.system.density_matrix,
                                optimize=True)

        num_occ = (self.model.ionic_charge // 2) * np.size(self._fock_matrix, 0) // self.model.orbitals_per_atom
        _, orbital_matrix = np.linalg.eigh(self._fock_matrix)
        occupied_matrix = orbital_matrix[:, :num_occ]
        self.density_matrix = occupied_matrix @ occupied_matrix.T

    @property
    def fock_matrix(self):
        return self._fock_matrix

    @fock_matrix.setter
    def fock_matrix(self, density_matrix):
        self._fock_matrix = self.system.hamiltonian_matrix.copy()
        self._fock_matrix += 2.0 * np.einsum(
            'pqt,rsu,tu,rs', self.system.chi_tensor, self.system.chi_tensor, self.system.interaction_matrix, density_matrix, optimize=True)
        self._fock_matrix -= np.einsum('rqt,psu,tu,rs',
                                self.system.chi_tensor,
                                self.system.chi_tensor,
                                self.system.interaction_matrix,
                                density_matrix,
                                optimize=True)

        num_occ = (self.model.ionic_charge // 2) * np.size(self.fock_matrix, 0) // self.model.orbitals_per_atom
        _, orbital_matrix = np.linalg.eigh(self.fock_matrix)
        occupied_matrix = orbital_matrix[:, :num_occ]
        self.density_matrix = occupied_matrix @ occupied_matrix.T

    def scf_cycle(self):
        """Calculate the density & Fock matrices

        Parameters
        ----------
        hamiltonian_matrix : np.array
        	The size is [n,n]
        interaction_matrix : np.array
        	 The size is [n,n]
        density_matrix : np.array
        	 The size is [n,n]
        chi_tensor : np.array
        max_scf_iterations : integer, default 100
        	The maximum number of scf cycles.
        mixing_fraction : float, default 0.25
        	The fraction of the previous density matrix incorporated into the next guess for the scf cycle
        convergence_tolerance : float, default 0.0001 (Hartree)
        	The energy difference between the current and previous scf cycle for which the convergence criterion
        	is satisfied

        Returns
        -------
        scf_cycle : np.array
        	Returns converged density & Fock matrices defined by the input Hamiltonian, interaction, & density matrices.

        """
        new_density_matrix = self.density_matrix.copy()
        for _ in range(self.max_scf_iterations):
            self.fock_matrix = new_density_matrix
            error_norm = np.linalg.norm(new_density_matrix - self.density_matrix)
            if error_norm < self.convergence_tolerance:
                return self.density_matrix, self.fock_matrix

            new_density_matrix = (self.mixing_fraction * self.density_matrix + (1.0 - self.mixing_fraction) * new_density_matrix)
        print("WARNING: SCF cycle didn't converge")
        return new_density_matrix, self.fock_matrix

    def calculate_energy_scf(self):
        '''the Hartree-Fock total energy defined by the input Hamiltonian, Fock, & density matrices.

        Parameters
        ----------
        hamiltonian_matrix : np.array
        	The size is [n,n]
        fock_matrix : np.array
        density_matrix : np.array

        Returns
        -------
        energy_scf : np.array
        	the Hartree-Fock total energy defined by the input Hamiltonian, Fock, & density matrices
        '''
        energy_scf = np.einsum('pq,pq', self.system.hamiltonian_matrix + self.fock_matrix, self.density_matrix)
        return energy_scf

    def get_hartree_fock_energy(self):
        return self.calculate_energy_scf() + self.system.energy_ion

class MP2(HartreeFock):
    def __init__(self, system, model):
        super().__init__(system, model)

    def get_partition_orbitals(self):
        self.occupied_energy, self.virtual_energy, self.occupied_matrix, self.virtual_matrix = self.partition_orbitals(
            self.fock_matrix, self.model)
        return self.occupied_energy, self.virtual_energy, self.occupied_matrix, self.virtual_matrix

    def get_mp2_energy(self):
        self.mp2_energy = self.calculate_energy_mp2(self.fock_matrix, self.system.interaction_matrix, self.system.chi_tensor,
                                                    self.model)
        return self.mp2_energy

    def partition_orbitals(self, fock_matrix, model):
        '''Returns a list with the occupied/virtual energies & orbitals defined by the input Fock matrix.

        Parameters
        ----------
        fock_matrix : np.array

        Returns
        -------
        occupied_energy : np.array
        virtual_energy : np.array
        occupies_matrix : np.array
        virtual_matrix : np.array
        '''
        num_occ = (model.ionic_charge // 2) * np.size(fock_matrix, 0) // model.orbitals_per_atom
        orbital_energy, orbital_matrix = np.linalg.eigh(fock_matrix)
        occupied_energy = orbital_energy[:num_occ]
        virtual_energy = orbital_energy[num_occ:]
        occupied_matrix = orbital_matrix[:, :num_occ]
        virtual_matrix = orbital_matrix[:, num_occ:]

        return occupied_energy, virtual_energy, occupied_matrix, virtual_matrix

    def transform_interaction_tensor(self, occupied_matrix, virtual_matrix, interaction_matrix, chi_tensor):
        ''''Calculates an interaction tensor.

        Parameters
        ----------
        occupied_energy : np.array
        virtual_energy : np.array
        interaction_matrix : np.array
        chi_tensor : np.array

        Returns
        --------
        A a transformed V tensor defined by the input occupied, virtual, & interaction matrices
        '''
        chi2_tensor = np.einsum('qa,ri,qrp', virtual_matrix, occupied_matrix, chi_tensor, optimize=True)
        interaction_tensor = np.einsum('aip,pq,bjq->aibj', chi2_tensor, interaction_matrix, chi2_tensor, optimize=True)
        return interaction_tensor

    def calculate_energy_mp2(self, fock_matrix, interaction_matrix, chi_tensor, model):
        '''Calculate the MP2 contribution to the total energy defined by the input Fock & interaction matrices.

        Parameters
        -----------
        fock_matrix : np.array
        interaction_matrix : np.array
        chi_tensor : np.array

        Returns
        -------
        the MP2 contribution to the total energy defined by the input Fock & interaction matrices
        '''
        E_occ, E_virt, occupied_matrix, virtual_matrix = self.partition_orbitals(fock_matrix, model)
        V_tilde = self.transform_interaction_tensor(occupied_matrix, virtual_matrix, interaction_matrix, chi_tensor)

        energy_mp2 = 0.0
        num_occ = len(E_occ)
        num_virt = len(E_virt)
        for a in range(num_virt):
            for b in range(num_virt):
                for i in range(num_occ):
                    for j in range(num_occ):
                        energy_mp2 -= ((2.0 * V_tilde[a, i, b, j]**2 - V_tilde[a, i, b, j] * V_tilde[a, j, b, i]) /
                                       (E_virt[a] + E_virt[b] - E_occ[i] - E_occ[j]))
        return energy_mp2
