import numpy as np

class AtomicSystem:
    def __init__(self,ionic_charge,orbital_types, orbital_occupation,vec):
        self.ionic_charge = ionic_charge
        self.orbital_types = orbital_types
        self.orbital_occupation = orbital_occupation
        self.p_orbitals = orbital_types[1:]
        self.orbitals_per_atom = len(orbital_types)

class StaticValues(AtomicSystem):
    def __init__(self,atomic_coordinates,model_parameters,ionic_charge,orbital_types, orbital_occupation,vec):
        .super()__init__(ionic_charge, orbital_types, orbital_occupation,vec)
        self.atomic_coordinates = atomic_coordinates
        self.model_parameters = model_parameters
        self.hamiltonian_matrix = self.calculate_hamiltonian_matrix()
        self.interaction_matrix = self.calculate_interaction_matrix()
        self.chi_tensor = self.calculate_chi_tensor()
        self.atomic_density = self.calculate_atomic_density_matrix()
        self.ndof = len(self.atomic_coordinates)*self.orbitals_per_atom
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
        for p in range(ndof):
            for q in range(ndof):
                if atom(p) != atom(q):
                    r_pq = self.atomic_coordinates[atom(p)] - self.atomic_coordinates[atom(q)]
                    interaction_matrix[p,q] = coulomb_energy(orb(p), orb(q), r_pq)
                if p == q and orb(p) == 's':
                    interaction_matrix[p,q] = self.model_parameters['coulomb_s']
                if p == q and orb(p) in self.p_orbitals:
                    interaction_matrix[p,q] = self.model_parameters['coulomb_p']
        return interaction_matrix

    def calculate_chi_tensor(atomic_coordinates, model_parameters):
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
            for orb_q in self.orbital_types:
                q = ao_index(atom(p), orb_q)
                for orb_r in self.orbital_types:
                    r = ao_index(atom(p), orb_r)
                    chi_tensor[p, q, r] = chi_on_atom(orb(p), orb(q), orb(r),
                                                      self.model_parameters)
        return chi_tensor

    def calculate_hamiltonian_matrix(atomic_coordinates, model_parameters):
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
        potential_vector = calculate_potential_vector(self.atomic_coordinates,
                                                      self.model_parameters)
        for p in range(self.ndof):
            for q in range(self.ndof):
                if atom(p) != atom(q):
                    r_pq = self.atomic_coordinates[atom(p)] - self.atomic_coordinates[atom(
                        q)]
                    hamiltonian_matrix[p, q] = hopping_energy(
                        orb(p), orb(q), r_pq, self.model_parameters)
                if atom(p) == atom(q):
                    if p == q and orb(p) == 's':
                        hamiltonian_matrix[p, q] += self.model_parameters['energy_s']
                    if p == q and orb(p) in p_orbitals:
                        hamiltonian_matrix[p, q] += self.model_parameters['energy_p']
                    for orb_r in self.orbital_types:
                        r = ao_index(atom(p), orb_r)
                        hamiltonian_matrix[p, q] += (
                            chi_on_atom(orb(p), orb(q), orb_r, self.model_parameters) *
                            potential_vector[r])
        return hamiltonian_matrix

    def calculate_atomic_density_matrix(atomic_coordinates):
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
        for p in range(self.ndof):
            density_matrix[p, p] = self.orbital_occupation[orb(p)]
        return density_matrix
