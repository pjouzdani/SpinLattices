# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

import numpy as np
import scipy.linalg as lg
#import matplotlib.pyplot as plt
from random import randint as rand
from random import uniform as randu
import Lattice3DWithExtField, Lattice3D


class  SpinLatticeInExternalField3DWithShape(Lattice3DWithExtField.SpinLatticeInExternalField3D, Lattice3D.SpinLattice3D):
    def __init__(self, interaction_type, 
                 symmetry,
                 boundary_type, 
                 n_column,  
                 n_row,  
                 n_height,  
                 j_field_tensor,  
                 total_spin_constant,
                 regularization_param,
                 permitivity,
                 external_field_tensor,
                 radius_on_top,
                 radius_on_bottom):
        self.interaction_type = interaction_type
        self.symmetry = symmetry
        self.n_column = n_column
        self.n_row = n_row
        self.n_height = n_height
        self.j_field_tensor =j_field_tensor
        self.boundary_type = boundary_type
        self.regularization_param = regularization_param
        self.total_spin_constant = total_spin_constant
        self.lattice = []
        #############
        self.lattice = np.ones((self.n_column, self.n_row, self.n_height))
        self.external_field_tensor = external_field_tensor
        self.permitivity = permitivity
        #        if radius_on_top<= np.sqrt(n_row*n_row + n_column * n_column)/2 and radius_on_top<= np.sqrt(n_row*n_row + n_column * n_column)/2:
        self.radius_on_top = radius_on_top
        self.radius_on_bottom = radius_on_bottom
        self.shape()
        self.num_bonds = self.get_lattice_num_bonds()
        # this will effectively calculate the valume of the sample - since at constructor call the lattice is all +1
        def _norm_fact():
            norm_fact = 0
            for idx_height in range(n_height):
                for idx_row in range(n_row):
                    for idx_column in range(n_column):
                        norm_fact = norm_fact + self.lattice[idx_column][idx_row][idx_height]
            return norm_fact
        #
        def _norm_fact_cylinder():
            n_column = self.n_column
            n_row = self.n_row
            n_height = self.n_height
            R_top = self.radius_on_top
            norm_fact_cylinder = 0
            for idx_height in range(n_height):
                for idx_row in range(n_row):
                    for idx_column in range(n_column):
                        Idx_row = np.sqrt(((idx_row+.5) - int(n_row/2.)) * ((idx_row+.5) -int( n_row/2.)) )
                        Idx_column = np.sqrt(((idx_column+.5) - int(n_column/2.)) * ((idx_column+.5) -int( n_column/2.)) )
                        r = np.sqrt(Idx_row*Idx_row + Idx_column*Idx_column) 
                        if(r<= R_top):
                            norm_fact_cylinder +=1
            return norm_fact_cylinder
        #
        self.norm_fact = _norm_fact()
        self.norm_fact_cylinder=_norm_fact_cylinder()
#
    def get_cylinder_polarization(self):
        n_column = self.n_column
        n_row = self.n_row
        n_height = self.n_height
        R_top = self.radius_on_top
        cylinder_polarization = 0
        for idx_height in range(n_height):
            for idx_row in range(n_row):
                for idx_column in range(n_column):
                    Idx_row = np.sqrt(((idx_row+.5) - int(n_row/2.)) * ((idx_row+.5) -int( n_row/2.)) )
                    Idx_column = np.sqrt(((idx_column+.5) - int(n_column/2.)) * ((idx_column+.5) -int( n_column/2.)) )
                    r = np.sqrt(Idx_row*Idx_row + Idx_column*Idx_column) 
                    if(r<= R_top):
                        cylinder_polarization = cylinder_polarization + self.lattice[idx_column][idx_row][idx_height]
        return cylinder_polarization
#   
    def shape(self):
        n_column = self.n_column
        n_row = self.n_row
        n_height = self.n_height
        R_top = self.radius_on_top
        #        print R_top
        R_bottom = self.radius_on_bottom
        for idx_height in range(n_height):
            for idx_row in range(n_row):
                for idx_column in range(n_column):
                    Idx_row = np.sqrt(((idx_row+.5) - int(n_row/2.)) * ((idx_row+.5) -int( n_row/2.)) )
                    Idx_column = np.sqrt(((idx_column+.5) - int(n_column/2.)) * ((idx_column+.5) -int( n_column/2.)) )
                    r = np.sqrt(Idx_row*Idx_row + Idx_column*Idx_column) 
                    r_max = R_bottom - (idx_height+1) * 1./n_height * (R_bottom - R_top)
                    if(r> r_max):
                        self.lattice[idx_column][idx_row][idx_height] = 0
#
    def get_lattice_energy(self):
        #
        energy = 0
        n_column = self.n_column
        n_row = self.n_row
        n_height = self.n_height        
        for idx_column in range(n_column):
            for idx_row in range(n_row):
                for idx_height in range(n_height):
                    Idx_row = np.sqrt(((idx_row+.5) - int(n_row/2.)) * ((idx_row+.5) -int( n_row/2.)) )
                    Idx_column = np.sqrt(((idx_column+.5) - int(n_column/2.)) * ((idx_column+.5) -int( n_column/2.)) )
                    r = np.sqrt(Idx_row*Idx_row + Idx_column*Idx_column) 
                    if r<=self.radius_on_bottom and r<=self.radius_on_top:
                        energy = energy + \
                            self.external_field_tensor[idx_column][idx_row][idx_height] * self.permitivity * \
                            self.lattice[idx_column][idx_row][idx_height].item()
        #                        
        # internal struct
        #
        #
        #############################################################################################################
        energy = energy + self.get_lattice_bare_energy()
        #############################################################################################################
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
        
        #############################################################################################################
        dE = -1 * self.get_local_bare_energy_change(idx_column,idx_row, idx_height)
        #############################################################################################################
        # contribution from external field:
        Idx_row = np.sqrt(((idx_row+.5) - int(n_row/2.)) * ((idx_row+.5) -int( n_row/2.)) )
        Idx_column = np.sqrt(((idx_column+.5) - int(n_column/2.)) * ((idx_column+.5) -int( n_column/2.)) )
        r = np.sqrt(Idx_row*Idx_row + Idx_column*Idx_column) 
        #
        if (r<=self.radius_on_bottom and r<=self.radius_on_top):
            dE = dE + 2 * self.lattice[idx_column][idx_row][idx_height] * \
            self.external_field_tensor[idx_column][idx_row][idx_height] * \
            self.permitivity
        return -dE

    def get_lattice_num_bonds(self):
        n_column = self.n_column
        n_row = self.n_row
        n_height = self.n_height        
        #
        num_bonds = 0
  
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
                                num_bonds = num_bonds + \
                                np.abs(self.lattice[idx_column][idx_row][idx_height])*np.abs(self.lattice[idx_column][idx_row+1][idx_height]) + \
                                np.abs(self.lattice[idx_column][idx_row][idx_height])*np.abs(self.lattice[idx_column+1][idx_row][idx_height]) + \
                                np.abs(self.lattice[idx_column][idx_row][idx_height])*np.abs(self.lattice[idx_column][idx_row][idx_height+1] )
            #
            # Boundary conditions:
                #     x = L 
                        #    h                          
                        #     \  
                        #      \
                        #      0___col
                        #
#                            #
                        if idx_row ==n_row-1 and \
                            idx_column !=n_column-1 and \
                            idx_height != n_height-1:
                                num_bonds = num_bonds + \
                                np.abs(self.lattice[idx_column][idx_row][idx_height])*np.abs(self.lattice[idx_column+1][idx_row][idx_height]) + \
                                np.abs(self.lattice[idx_column][idx_row][idx_height])*np.abs(self.lattice[idx_column][idx_row][idx_height+1]) 
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
                                num_bonds = num_bonds + \
                                np.abs( self.lattice[idx_column][idx_row][idx_height])*\
                                        np.abs(self.lattice[idx_column][idx_row][idx_height+1] )
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
                                num_bonds = num_bonds + \
                                np.abs( self.lattice[idx_column][idx_row][idx_height])*\
                                        np.abs(self.lattice[idx_column+1][idx_row][idx_height] )
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
                                num_bonds = num_bonds + \
                                np.abs(self.lattice[idx_column][idx_row][idx_height])*\
                                    np.abs(self.lattice[idx_column][idx_row+1][idx_height]) + \
                                np.abs(self.lattice[idx_column][idx_row][idx_height])*\
                                    np.abs(self.lattice[idx_column][idx_row][idx_height+1] )
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
                                num_bonds = num_bonds + \
                                np.abs(self.lattice[idx_column][idx_row][idx_height])*\
                                    np.abs(self.lattice[idx_column][idx_row+1][idx_height])
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
                                num_bonds = num_bonds + \
                                np.abs(self.lattice[idx_column][idx_row][idx_height])*\
                                    np.abs(self.lattice[idx_column][idx_row+1][idx_height]) + \
                                np.abs(self.lattice[idx_column][idx_row][idx_height])*\
                                    np.abs(self.lattice[idx_column+1][idx_row][idx_height]) 
    ##                            
        return num_bonds 
   
    def get_lattice_bare_energy(self):
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
                                self.j_field_tensor[idx_column][idx_row][idx_height][0] *self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row+1][idx_height] + \
                                self.j_field_tensor[idx_column][idx_row][idx_height][1] *self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column+1][idx_row][idx_height] + \
                                self.j_field_tensor[idx_column][idx_row][idx_height][2] *self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row][idx_height+1] 
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
                                self.j_field_tensor[idx_column][idx_row][idx_height][0] * self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column+1][idx_row][idx_height] + \
                                self.j_field_tensor[idx_column][idx_row][idx_height][2] * self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row][idx_height+1] 
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
                                self.j_field_tensor[idx_column][idx_row][idx_height][2] * self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row][idx_height+1] 
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
                                self.j_field_tensor[idx_column][idx_row][idx_height][0] * self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column+1][idx_row][idx_height] 
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
                                self.j_field_tensor[idx_column][idx_row][idx_height][1] *self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row+1][idx_height] + \
                                self.j_field_tensor[idx_column][idx_row][idx_height][2]* self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row][idx_height+1] 
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
                                self.j_field_tensor[idx_column][idx_row][idx_height][1] * self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row+1][idx_height]
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
                                self.j_field_tensor[idx_column][idx_row][idx_height][1] * self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column][idx_row+1][idx_height] + \
                                self.j_field_tensor[idx_column][idx_row][idx_height][0] *self.lattice[idx_column][idx_row][idx_height]*self.lattice[idx_column+1][idx_row][idx_height] 
                                #
                        #
                        #
        return energy
    #
    #        
    def get_local_bare_energy_change(self, idx_column,idx_row, idx_height):
    # for given i and j , and k calc flipping energy -> \delta E
        # return         - 2 \delta E
        n_column = self.n_column
        n_row = self.n_row   
        n_height = self.n_height
        if idx_column>=0 and idx_column<n_column and \
            idx_row>=0 and idx_row<n_row and\
            idx_height>=0 and idx_height<n_height:
            dE = 0
            if idx_column != n_column-1:
                dE = dE + \
                self.lattice[idx_column][idx_row][idx_height].item() * \
                self.lattice[idx_column+1][idx_row][idx_height].item() * self.j_field_tensor[idx_column][idx_row][idx_height][0]
            if idx_column != 0:
                dE = dE + \
                self.lattice[idx_column][idx_row][idx_height].item() *\
                self.lattice[idx_column-1][idx_row][idx_height].item() * self.j_field_tensor[idx_column-1][idx_row][idx_height][0]
            if idx_row != n_row-1:
                dE = dE + \
                self.lattice[idx_column][idx_row][idx_height].item() *\
                self.lattice[idx_column][idx_row+1][idx_height].item() * self.j_field_tensor[idx_column][idx_row][idx_height][1]
            if idx_row != 0:
                dE = dE + \
                self.lattice[idx_column][idx_row][idx_height].item() *\
                self.lattice[idx_column][idx_row-1][idx_height].item() * self.j_field_tensor[idx_column][idx_row-1][idx_height][1]
            if idx_height != n_height-1:
                dE = dE + \
                self.lattice[idx_column][idx_row][idx_height].item() *\
                self.lattice[idx_column][idx_row][idx_height+1].item() * self.j_field_tensor[idx_column][idx_row][idx_height][2]
            if idx_height != 0:
                dE = dE + \
                self.lattice[idx_column][idx_row][idx_height].item() *\
                self.lattice[idx_column][idx_row][idx_height-1].item() * self.j_field_tensor[idx_column][idx_row][idx_height-1][2]
            #
            #
            
            return -2*dE
        else:
            return 'Error out of bound indixes'
    #
    #
    def flip_dipole(self, idx_column,idx_row, idx_height):
        return Lattice3D.SpinLattice3D.flip_dipole(self, idx_column,idx_row, idx_height)
    #
    #
    def get_local_spin_change_effect(self, idx_column,idx_row, idx_height):
        total_spin = self.get_total_spin()
        total_spin_constant = self.total_spin_constant
        # takes into account the conservation of the total spin (\sum_i S_i )^2 ~ const
#        return self.regularization_param * \
#        Sigma_i = (total_spin*total_spin - total_spin_constant*total_spin_constant )
        delta_k = 2*self.lattice[idx_column][idx_row][idx_height]
        dS = 0
        dS = self.regularization_param *(np.sqrt( ((total_spin - delta_k)* (total_spin - delta_k) - \
                                                        total_spin_constant*total_spin_constant) * \
                                                 ((total_spin - delta_k)* (total_spin - delta_k) - \
                                                 total_spin_constant*total_spin_constant) )-\
                                         np.sqrt((total_spin*total_spin - \
                                                 total_spin_constant*total_spin_constant) * \
                                         (total_spin*total_spin - total_spin_constant*total_spin_constant)))
        return dS