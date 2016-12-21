# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
import numpy as np
import scipy.linalg as lg
#import matplotlib.pyplot as plt
from random import randint as rand
from random import uniform as randu
import Lattice



class  SpinLattice3D(Lattice.SpinLatticeInExternalField):
    def __init__(self, interaction_type, 
                 symmetry,
                 boundary_type, 
                 n_column,  
                 n_row,  
                 n_height,  
                 j_column_column,  
                 j_row_row,  
                 j_height_height,
                 total_spin_constant,
                 regularization_param):
        self.interaction_type = interaction_type
        self.symmetry = symmetry
        self.n_column = n_column
        self.n_row = n_row
        self.n_height = n_height
        self.j_column_column  =j_column_column
        self.j_row_row  =j_row_row
        self.j_height_height  =j_height_height
        self.boundary_type = boundary_type
        self.regularization_param = regularization_param
        self.lattice = []
        #############
        self.lattice = np.ones((self.n_column, self.n_row, self.n_height))
#        lattice_2D = []
#        for idx_col in range(self.n_column):
#            lattice_2D.append(np.ones(self.n_height))
#        #
#        for idx_row in range(self.n_row):
#            self.lattice.append(lattice_2D)
        def num_bonds():
            if self.symmetry == 'Z2':
                if self.interaction_type == 'nn':
                    if self.boundary_type == 'np':
                        return self.n_height* ( (self.n_row)*(self.n_column -1)  +\
                         (self.n_column) * (self.n_row-1) ) + \
                         self.n_row * self.n_column * (self.n_height-1)
        self.num_bonds = num_bonds()
        #####  constraint on the S_total
        self.total_spin_constant = total_spin_constant
        self.total_spin = self.get_total_spin()
        ############
#        ############
#        self.external_field_profile = []
#        external_field_profile_2D = []
#        for idx_col in range(self.n_column):
#            external_field_profile_2D.append(np.ones(self.n_height))
#        #
#        for idx_row in range(self.n_row):
#            self.external_field_profile.append(external_field_profile_2D)
#        #
#        #
#        for idx_row in range(self.n_row):
#            for idx_column in range(self.n_column):
#                for idx_height in range(self.n_height):
#                    self.external_field_profile[idx_row][idx_column][idx_height] =\
#                        self.external_field_profile[idx_row][idx_column][idx_height]* \
#                        self.external_field_magnitude
#        ##########
    #
    #
    #
    #
    ###############
    def get_total_spin(self):
        total_spin  = 0
        for idx_column in range(self.n_column):
            for idx_row in range(self.n_row):
                for idx_height in range(self.n_height):
                    total_spin = total_spin + self.lattice[idx_column][idx_row][idx_height]
        return total_spin
#        ##########
    #
    #
    #
    #
    ###############    
    
    def get_lattice_energy(self):
        n_column = self.n_column
        n_row = self.n_row
        n_height = self.n_height        
        #
        energy = 0

        # internal struct
        #
        #
        if self.boundary_type =='np' and self.interaction_type=='nn':
            for idx_row in range(n_row):
                for idx_column in range(n_column):
                    for idx_height in range(n_height):
                    
                        #      Bulk
                        #
                        #    h    row                      
                        #     \  /
                        #      \/
                        #      0___col
                        #
                        if idx_row !=n_row-1 and \
                            idx_column !=n_column-1 and \
                            idx_height != n_height-1:
                                energy = energy + \
                                self.j_row_row *self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row+1][idx_height] + \
                                self.j_column_column *self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column+1][idx_row][idx_height] + \
                                self.j_height_height *self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row][idx_height+1] 
            #
            # Boundary conditions:
                #     x = L 
                        #    h                          
                        #     \  
                        #      \
                        #      0___col
                        #
                        #
                        if idx_row ==n_row-1 and \
                            idx_column !=n_column-1 and \
                            idx_height != n_height-1:
                                energy = energy + \
                                self.j_column_column * self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column+1][idx_row][idx_height] + \
                                self.j_height_height * self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row][idx_height+1] 
                        #
                        #    h                      
                        #     \  
                        #      \
                        #      0
                        #
                        #
                        if idx_row ==n_row-1 and \
                            idx_column ==n_column-1 and \
                            idx_height != n_height-1:
                                energy = energy + \
                                self.j_height_height * self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row][idx_height+1] 
                        #
                        #                          
                        #     
                        #      
                        #      0___col
                        #
                        #
                        if idx_row ==n_row-1 and \
                            idx_column !=n_column-1 and \
                            idx_height == n_height-1:
                                energy = energy + \
                                self.j_column_column * self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column+1][idx_row][idx_height] 
                #     y = L 
                        #
                        #
                        #    h    row                      
                        #     \  /
                        #      \/
                        #      0
                        #
                        #
                        if idx_row !=n_row-1 and \
                            idx_column ==n_column-1 and \
                            idx_height != n_height-1:
                                energy = energy + \
                                self.j_row_row *self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row+1][idx_height] + \
                                self.j_height_height * self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row][idx_height+1] 
                        #
                        #         row                      
                        #        /
                        #       /
                        #      0
                        #
                        #
                        if idx_row !=n_row-1 and \
                            idx_column ==n_column-1 and \
                            idx_height == n_height-1:                        
                                energy = energy + \
                                self.j_row_row * self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row+1][idx_height]
                                #
                #     z = L 
                        #
                        #
                        #         row                      
                        #        /
                        #       /
                        #      0___col
                        #
                        #
                        if idx_row != n_row-1 and \
                            idx_column != n_column-1 and \
                            idx_height == n_height-1:
                                energy = energy + \
                                self.j_row_row * self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row+1][idx_height] + \
                                self.j_column_column *self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column+1][idx_row][idx_height] 
                                #
                        #
                        #
        return energy     
    #
    #
    def get_local_energy_change(self, idx_column,idx_row, idx_height):
    # for given i and j , and k calc flipping energy -> \delta E
        # return         - 2 \delta E
        n_column = self.n_column
        n_row = self.n_row   
        n_height = self.n_height
        j_row_row = self.j_row_row
        j_height_height = self.j_height_height
        j_column_column = self.j_column_column
        if idx_column>=0 and idx_column<n_column and \
            idx_row>=0 and idx_row<n_row and\
            idx_height>=0 and idx_height<n_height:
            dE = 0
            if idx_column != n_column-1:
                dE = dE + \
                self.lattice[idx_column][idx_row][idx_height].item() * \
                self.lattice[idx_column+1][idx_row][idx_height].item() * j_column_column
            if idx_column != 0:
                dE = dE + \
                self.lattice[idx_column][idx_row][idx_height].item() *\
                self.lattice[idx_column-1][idx_row][idx_height].item() * j_column_column
            if idx_row != n_row-1:
                dE = dE + \
                self.lattice[idx_column][idx_row][idx_height].item() *\
                self.lattice[idx_column][idx_row+1][idx_height].item() * j_row_row
            if idx_row != 0:
                dE = dE + \
                self.lattice[idx_column][idx_row][idx_height].item() *\
                self.lattice[idx_column][idx_row-1][idx_height].item() * j_row_row
            if idx_height != n_height-1:
                dE = dE + \
                self.lattice[idx_column][idx_row][idx_height].item() *\
                self.lattice[idx_column][idx_row][idx_height+1].item() * j_height_height
            if idx_height != 0:
                dE = dE + \
                self.lattice[idx_column][idx_row][idx_height].item() *\
                self.lattice[idx_column][idx_row][idx_height-1].item() * j_height_height
            #
            #
            
            return -2*dE
        else:
            return 'Error out of bound indixes'
##
#    def get_local_spin_change_effect(self, idx_column,idx_row, idx_height):
#        total_spin=self.total_spin
#        total_spin_constant = self.total_spin_constant
#        # takes into account the conservation of the total spin (\sum_i S_i )^2 ~ const
##        return self.regularization_param * \
##        Sigma_i = (total_spin*total_spin - total_spin_constant*total_spin_constant )
#        delta_k = 2*self.lattice[idx_column][idx_row][idx_height]
#        dS = 0
#        dS = self.regularization_param *(np.sqrt( ((total_spin - delta_k)* (total_spin - delta_k) - \
#                                                        total_spin_constant*total_spin_constant) * \
#                                                 ((total_spin - delta_k)* (total_spin - delta_k) - \
#                                                 total_spin_constant*total_spin_constant) )-\
#                                         np.sqrt((total_spin*total_spin - \
#                                                 total_spin_constant*total_spin_constant) * \
#                                         (total_spin*total_spin - total_spin_constant*total_spin_constant)))
#        return dS
                
    #
    #
    def flip_dipole(self, idx_column,idx_row, idx_height):
        if self.lattice[idx_column][idx_row][idx_height]== 1.:
            self.lattice[idx_column][idx_row][idx_height] = -1.
        else:
            self.lattice[idx_column][idx_row][idx_height] = +1.
        
    
         