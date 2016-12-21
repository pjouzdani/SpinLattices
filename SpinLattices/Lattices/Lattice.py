
# -*- coding: utf-8 -*-
import numpy as np
import scipy.linalg as lg
#import matplotlib.pyplot as plt
from random import randint as rand
from random import uniform as randu

class Lattice():
    def  __init__(self, 
                  dimension):
      #
      self.dimension =dimension
      self.lattice=[]
      
class SpinLattice(Lattice):
    def  __init__(self,interaction_type ,symmetry,boundary_type,  n_column ,n_row, j_column_column,j_row_row ):
        self.interaction_type = interaction_type
        self.symmetry = symmetry
        self.n_column = n_column
        self.n_row = n_row 
        self.j_column_column  =j_column_column
        self.j_row_row  =j_row_row
        self.boundary_type = boundary_type
        self.lattice = []
        for idx_row in range(self.n_row):
            self.lattice.append(np.ones(self.n_column))
        def num_bonds():
            if self.symmetry == 'Z2':
                if self.interaction_type == 'nn':
                    if self.boundary_type == 'np':
                        return (self.n_row)*(self.n_column -1)  + (self.n_column) * (self.n_row-1)
        self.num_bonds = num_bonds() 
        
    def get_local_energy_change(self, idx_column,idx_row):
    # for given i and j calc flipping energy -> \delta E
        # return 2 \delta E
        n_column = self.n_column
        n_row = self.n_row   
        j_row_row = self.j_row_row
        j_column_column = self.j_column_column
        if idx_column>=0 and idx_column<n_column and idx_row>=0 and idx_row<n_row:
            dE = 0
            if idx_column != n_column-1:
                dE = dE + \
                self.lattice[idx_row][idx_column].item() * self.lattice[idx_row][idx_column+1].item() * j_column_column
            if idx_column != 0:
                dE = dE + \
            self.lattice[idx_row][idx_column].item() * self.lattice[idx_row][idx_column-1].item() * j_column_column
            if idx_row != n_row-1:
                dE = dE + \
                self.lattice[idx_row][idx_column].item() * self.lattice[idx_row+1][idx_column].item() * j_row_row
            if idx_row != 0:
                dE = dE + \
                self.lattice[idx_row][idx_column].item() * self.lattice[idx_row-1][idx_column].item() * j_row_row
            return -2*dE
        else:
            return 'Error out of bound indixes'
        #
    def flip_dipole(self, idx_column,idx_row):
        if self.lattice[idx_row][idx_column]== 1.:
            self.lattice[idx_row][idx_column] = -1.
        else:
            self.lattice[idx_row][idx_column] = +1.
        
                     
                        
    #       methods ##
    def get_lattice_energy(self):
        n_column = self.n_column
        n_row = self.n_row
        j_row_row = self.j_row_row
        j_column_column = self.j_column_column
        
        energy = 0        
        for idx_row in range(n_row-1):
            energy = energy + j_row_row * self.lattice[idx_row][n_column-1].item()*\
                   self.lattice[idx_row+1][n_column-1]
                   
        for idx_column in range(n_column-1):
            energy = energy + j_column_column * self.lattice[n_row-1][idx_column].item()*\
            self.lattice[n_row-1][idx_column+1]
            
        for idx_column in range(n_column-1):
            for idx_row in range(n_row-1):
                energy = j_column_column * self.lattice[idx_row][idx_column].item() * \
                               self.lattice[idx_row][idx_column+1].item()+\
                        j_row_row  * self.lattice[idx_row][idx_column].item()*\
                               self.lattice[idx_row+1][idx_column+1].item() + \
                               energy
        return energy
#####################
#
#
# Spin lattice in an external Field        
#
#
#####################        
class SpinLatticeInExternalField(SpinLattice):
    def __init__(self, interaction_type,
                             symmetry,boundary_type,  
                             n_column ,
                             n_row, 
                             j_column_column,
                             j_row_row, 
                             external_field_magnitude):
        SpinLattice.__init__(self, 
                             interaction_type,
                             symmetry,boundary_type,  
                             n_column ,
                             n_row, 
                             j_column_column,
                             j_row_row)
        self.external_field_magnitude = external_field_magnitude
        self.external_field_profile = []
        for idx_row in range(self.n_row):
            self.external_field_profile.append(np.ones(self.n_column))
        for idx_row in range(self.n_row):
            for idx_column in range(self.n_column):
                self.external_field_profile[idx_row][idx_column] =\
                        self.external_field_profile[idx_row][idx_column]* \
                        self.external_field_magnitude
    #
    def get_lattice_energy(self):
        #
        n_column = self.n_column
        n_row = self.n_row
        #
        energy = SpinLattice.get_lattice_energy(self)
        #
        for idx_row in range(n_row):
            for idx_column in range(n_column):
                energy = energy + self.external_field_magnitude * self.lattice[idx_row][idx_column].item()
        return energy
    #            
    def get_local_energy_change(self, idx_column,idx_row):
    # for given i and j calc flipping energy -> \delta E
        # return 2 \delta E
        #
        dE = SpinLattice.get_local_energy_change(self, idx_column, idx_row) + \
                self.lattice[idx_row][idx_column] * self.external_field_magnitude
        
        return dE        
        
        
        