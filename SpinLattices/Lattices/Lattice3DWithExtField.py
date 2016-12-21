# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
import numpy as np
import scipy.linalg as lg
#import matplotlib.pyplot as plt
from random import randint as rand
from random import uniform as randu
import Lattice3D
import Lattice



class  SpinLatticeInExternalField3D(Lattice3D.SpinLattice3D):
    def __init__(self, interaction_type, 
                 symmetry,
                 boundary_type, 
                 n_column,  
                 n_row,  
                 n_height,  
                 j_column_column,  
                 j_row_row,  
                 j_height_height, 
                 total_spin,
                 regularization_param,
                 permitivity,
                 external_field_magnitude):
        Lattice3D.SpinLattice3D.__init__(self, interaction_type, 
                 symmetry,
                 boundary_type, 
                 n_column,  
                 n_row,  
                 n_height,  
                 j_column_column,  
                 j_row_row,  
                 j_height_height,
                 total_spin,
                 regularization_param)
        self.external_field_magnitude = external_field_magnitude
        self.permitivity = permitivity
#   
    def get_lattice_energy(self):
        #
        energy = 0
        n_column = self.n_column
        n_row = self.n_row
        n_height = self.n_height        
#        # external field
        for idx_column in range(n_column):
            for idx_row in range(n_row):
                for idx_height in range(n_height):
                    energy = energy + \
                        self.external_field_magnitude * self.permitivity * \
                        self.lattice[idx_column][idx_row][idx_height].item()
        #                        
        # internal struct
        #
        #
        energy = energy + Lattice3D.SpinLattice3D.get_lattice_energy(self)
        #
        return energy
#
#
    def get_local_energy_change(self, idx_column,idx_row, idx_height):
    # for given i and j calc flipping energy -> \delta E
        # return 2 \delta E
        dE = 0
        #
        n_column = self.n_column
        n_row = self.n_row   
        n_height = self.n_height
        j_row_row = self.j_row_row
        j_height_height = self.j_height_height
        j_column_column = self.j_column_column
        dE = -1 * Lattice3D.SpinLattice3D.get_local_energy_change(self, idx_column,idx_row, idx_height)
        dE = dE + 2 * self.lattice[idx_column][idx_row][idx_height] * self.external_field_magnitude * self.permitivity
        return -dE
