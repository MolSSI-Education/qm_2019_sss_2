"""
atomicmodel.py
Refactoring of the QM project from MolSSI Software Summer School 2019

Handles the primary functions
"""
import numpy as np
class Model:
    def __init__(self,model_parameters,ionic_charge,orbital_types, orbital_occupation,vec):
        self.model_parameters = model_parameters
        self.ionic_charge = ionic_charge
        self.orbital_types = orbital_types
        self.orbital_occupation = orbital_occupation
        self.p_orbitals = orbital_types[1:]
        self.orbitals_per_atom = len(orbital_types)
        self.vec = vec

class System:
    def __init__(self,atomic_coordinates,model):
        self.atomic_coordinates = atomic_coordinates
        self.ndof = len(self.atomic_coordinates)*model.orbitals_per_atom
        self.orbital = []
        for i in range(self.ndof):
            atom= int(np.floor(i / model.orbitals_per_atom))
            orbital_num = i % model.orbitals_per_atom
            self.orbital.append([atom,model.orbital_types[orbital_num]])
        self.interaction_matrix = self.calculate_interaction_matrix()
        self.chi_tensor = self.calculate_chi_tensor(model)
        self.density_matrix = self.calculate_atomic_density_matrix(model)
        self.hamiltonian_matrix = self.calculate_hamiltonian_matrix(model)
        self.energy_ion = self.calculate_energy_ion(model)
    def calculate_energy_ion(self,model):
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
                    energy_ion += (model.ionic_charge**2) * self.coulomb_energy(
                        's', 's', r_i - r_j,model)
        return energy_ion

    def orb(self,index):
        return self.orbital[index][1]
    def atom(self,index):
        return self.orbital[index][0]
    def calculate_interaction_matrix(self):
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

        interaction_matrix = np.zeros( (self.ndof,self.ndof) )
        for p in range(self.ndof):
            for q in range(self.ndof):
                if self.atom(p) != self.atom(q):
                    r_pq = self.atomic_coordinates[self.atom(p)] - self.atomic_coordinates[self.atom(q)]
                    interaction_matrix[p,q] = self.coulomb_energy(self.orb(p), self.orb(q), r_pq, model)
                if p == q and self.orb(p) == 's':
                    interaction_matrix[p,q] = model.model_parameters['coulomb_s']
                if p == q and self.orb(p) in model.p_orbitals:
                    interaction_matrix[p,q] = model.model_parameters['coulomb_p']
        return interaction_matrix
    def calculate_potential_vector(self, model):
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
                    potential_vector[p] += (
                        self.pseudopotential_energy(self.orb(p), r_pi, model) -
                        ionic_charge * self.coulomb_energy(self.orb(p), 's', r_pi, model))
        return potential_vector
    def coulomb_energy(self, o1, o2, r12, model):
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
        if o1 == 's' and o2 in model.p_orbitals:
            ans = np.dot(model.vec[o2], r12) / r12_length**3
        if o2 == 's' and o1 in model.p_orbitals:
            ans = -1 * np.dot(model.vec[o1], r12) / r12_length**3
        if o1 in model.p_orbitals and o2 in model.p_orbitals:
            ans = (
                np.dot(model.vec[o1], model.vec[o2]) / r12_length**3 -
                3.0 * np.dot(model.vec[o1], r12) * np.dot(model.vec[o2], r12) / r12_length**5)
        return ans
    def pseudopotential_energy(self, o, r, model):
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
        ans = model.model_parameters['v_pseudo']
        r_rescaled = r / model.model_parameters['r_pseudo']
        ans *= np.exp(1.0 - np.dot(r_rescaled, r_rescaled))
        if o in model.p_orbitals:
            ans *= -2.0 * np.dot(model.vec[o], r_rescaled)
        return ans
    def hopping_energy(self,o1, o2, r12, model):
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
        r12_rescaled = r12 / model.model_parameters['r_hop']
        r12_length = np.linalg.norm(r12_rescaled)
        ans = np.exp( 1.0 - r12_length**2 )
        if o1 == 's' and o2 == 's':
            ans *= model.model_parameters['t_ss']
        if o1 == 's' and o2 in model.p_orbitals:
            ans *= np.dot(model.vec[o2], r12_rescaled) * model.model_parameters['t_sp']
        if o2 == 's' and o1 in model.p_orbitals:
            ans *= -np.dot(model.vec[o1], r12_rescaled)* model.model_parameters['t_sp']
        if o1 in model.p_orbitals and o2 in model.p_orbitals:
            ans *= ( (r12_length**2) * np.dot(vec[o1], vec[o2]) * model.model_parameters['t_pp2']
                     - np.dot(model.vec[o1], r12_rescaled) * np.dot(model.vec[o2], r12_rescaled)
                     * ( model.model_parameters['t_pp1'] + model.model_parameters['t_pp2'] ) )
        return ans

    def calculate_chi_tensor(self, model):
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
            for orb_q in model.orbital_types:
                q = p % model.orbitals_per_atom + model.orbital_types.index(orb_q)
                for orb_r in model.orbital_types:
                    r = p % model.orbitals_per_atom + model.orbital_types.index(orb_r)
                    chi_tensor[p, q, r] = self.chi_on_atom(self.orb(p), self.orb(q), self.orb(r), model)
        return chi_tensor
    def calculate_hamiltonian_matrix(self, model):
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
        potential_vector = self.calculate_potential_vector(model)
        for p in range(self.ndof):
            for q in range(self.ndof):
                if self.atom(p) != self.atom(q):
                    r_pq = self.atomic_coordinates[self.atom(p)] - self.atomic_coordinates[self.atom(
                        q)]
                    hamiltonian_matrix[p, q] = self.hopping_energy(
                        self.orb(p), self.orb(q), r_pq, model)
                if self.atom(p) == self.atom(q):
                    if p == q and self.orb(p) == 's':
                        hamiltonian_matrix[p, q] += model.model_parameters['energy_s']
                    if p == q and self.orb(p) in model.p_orbitals:
                        hamiltonian_matrix[p, q] += model.model_parameters['energy_p']
                    for orb_r in model.orbital_types:
                        r = p % model.orbitals_per_atom + model.orbital_types.index(orb_r)
                        hamiltonian_matrix[p, q] += (
                            self.chi_on_atom(self.orb(p), self.orb(q), orb_r, model) *
                            potential_vector[r])
        return hamiltonian_matrix
    def calculate_atomic_density_matrix(self,model):
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
            density_matrix[p, p] = model.orbital_occupation[self.orb(p)]
            print(F"({p},{p}): {self.orb(p)}, {model.orbital_occupation[self.orb(p)]}")
        return density_matrix
    def chi_on_atom(self, o1, o2, o3, model):
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
        if o1 == o3 and o3 in model.p_orbitals and o2 == 's':
            return model.model_parameters['dipole']
        if o2 == o3 and o3 in model.p_orbitals and o1 == 's':
            return model.model_parameters['dipole']
        return 0.0

class HartreeFock:
    def __init__(self,system,model):
        self.fock_matrix= self.calculate_fock_matrix(system.hamiltonian_matrix, system.interaction_matrix, system.density_matrix, system.chi_tensor)
        self.density_matrix= self.calculate_density_matrix(self.fock_matrix)
        self.max_scf_iterations = 100
        self.convergence_tolerance = 1e-4
        self.mixing_fraction = 0.25
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
    def calculate_density_matrix(self, fock_matrix):
        '''Returns the 1-electron density matrix defined by the input Fock matrix.'''
        num_occ = (model.ionic_charge // 2) * np.size(fock_matrix,
                                                0) // model.orbitals_per_atom
        orbital_energy, orbital_matrix = np.linalg.eigh(fock_matrix)
        occupied_matrix = orbital_matrix[:, :num_occ]
        density_matrix = occupied_matrix @ occupied_matrix.T
        return density_matrix
    def set_max_iterations(self,max_iter):
        self.max_scf_iterations = max_iter
    def set_convergence_tolerance(self, tolerance):
        self.convergence_tolerance = tolerance
    def set_mixing_fraction(self,fraction):
        self.mixing_fraction = fraction
    def run_scf(self):
        self.density_matrix, self.fock_matrix = self.scf_cycle( system.hamiltonian_matrix, system.interaction_matrix, self.density_matrix, system.chi_tensor, self.max_scf_iterations, self.mixing_fraction, self.convergence_tolerance)
    def scf_cycle(self, hamiltonian_matrix, interaction_matrix, density_matrix, chi_tensor, max_scf_iterations, mixing_fraction, convergence_tolerance):
        '''Calculate the density & Fock matrices

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

        '''
        old_density_matrix = density_matrix.copy()
        for iteration in range(max_scf_iterations):
            new_fock_matrix = self.calculate_fock_matrix(hamiltonian_matrix, interaction_matrix, old_density_matrix, chi_tensor)
            new_density_matrix = self.calculate_density_matrix(new_fock_matrix)

            error_norm = np.linalg.norm( old_density_matrix - new_density_matrix )
            if error_norm < convergence_tolerance:
                return new_density_matrix, new_fock_matrix

            old_density_matrix = (mixing_fraction * new_density_matrix
                                  + (1.0 - mixing_fraction) * old_density_matrix)
        print("WARNING: SCF cycle didn't converge")
        return new_density_matrix, new_fock_matrix
    def get_scf_energy(self):
        return self.calculate_energy_scf(system.hamiltonian_matrix,self.fock_matrix,self.density_matrix)
    def calculate_energy_scf(self,hamiltonian_matrix, fock_matrix, density_matrix):
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
        energy_scf = np.einsum('pq,pq', hamiltonian_matrix + fock_matrix,
                               density_matrix)
        return energy_scf
    def get_hartree_fock_energy(self):
        return self.get_scf_energy() + system.energy_ion
class MP2(HartreeFock):
    def __init__(self,system,model):
        super().__init__(system,model)
    def get_partition_orbitals(self):
        self.occupied_energy, self.virtual_energy, self.occupied_matrix, self.virtual_matrix=self.partition_orbitals(self.fock_matrix,model)
        return self.occupied_energy, self.virtual_energy, self.occupied_matrix, self.virtual_matrix
    def get_mp2_energy(self):
        self.mp2_energy = self.calculate_energy_mp2(self.fock_matrix,system.interaction_matrix, system.chi_tensor,model)
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
        num_occ = (model.ionic_charge // 2) * np.size(fock_matrix,
                                                0) // model.orbitals_per_atom
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
        chi2_tensor = np.einsum('qa,ri,qrp',
                                virtual_matrix,
                                occupied_matrix,
                                chi_tensor,
                                optimize=True)
        interaction_tensor = np.einsum('aip,pq,bjq->aibj',
                                       chi2_tensor,
                                       interaction_matrix,
                                       chi2_tensor,
                                       optimize=True)
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
        E_occ, E_virt, occupied_matrix, virtual_matrix = self.partition_orbitals(
            fock_matrix, model)
        V_tilde = self.transform_interaction_tensor(occupied_matrix, virtual_matrix,
                                               interaction_matrix, chi_tensor)

        energy_mp2 = 0.0
        num_occ = len(E_occ)
        num_virt = len(E_virt)
        for a in range(num_virt):
            for b in range(num_virt):
                for i in range(num_occ):
                    for j in range(num_occ):
                        energy_mp2 -= (
                            (2.0 * V_tilde[a, i, b, j]**2 -
                             V_tilde[a, i, b, j] * V_tilde[a, j, b, i]) /
                            (E_virt[a] + E_virt[b] - E_occ[i] - E_occ[j]))
        return energy_mp2

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    # pass
    # --------------------
    # Noble Gas Parameters
    # --------------------
    ionic_charge = 6
    orbital_types = ['s', 'px', 'py', 'pz']
    vec = {'px': [1, 0, 0], 'py': [0, 1, 0], 'pz': [0, 0, 1]}
    orbital_occupation = { 's':0, 'px':1, 'py':1, 'pz':1 }
    model_parameters = {
    'r_hop' : 3.1810226927827516,
    't_ss' : 0.03365982238611262,
    't_sp' : -0.029154833035109226,
    't_pp1' : -0.0804163845390335,
    't_pp2' : -0.01393611496959445,
    'r_pseudo' : 2.60342991362958,
    'v_pseudo' : 0.022972992186364977,
    'dipole' : 2.781629275106456,
    'energy_s' : 3.1659446174413004,
    'energy_p' : -2.3926873325346554,
    'coulomb_s' : 0.3603533286088998,
    'coulomb_p' : -0.003267991835806299
    }
    # First Instantiation of Model class.
    model = Model(model_parameters,ionic_charge,orbital_types, orbital_occupation,vec)
    atomic_coordinates = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 5.0]])
    system = System(atomic_coordinates,model)
    energy = MP2(system,model)
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
