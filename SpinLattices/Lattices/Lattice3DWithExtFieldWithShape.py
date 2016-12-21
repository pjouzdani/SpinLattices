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
                 j_column_column,  
                 j_row_row,  
                 j_height_height, 
                 total_spin,
                 regularization_param,
                 permitivity,
                 external_field_magnitude,
                 radius_on_top,
                 radius_on_bottom):
        Lattice3DWithExtField.SpinLatticeInExternalField3D.__init__(self, interaction_type, 
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
                 external_field_magnitude)
        self.external_field_magnitude = external_field_magnitude
        self.permitivity = permitivity
        if radius_on_top<= np.sqrt(n_row*n_row + n_column * n_column)/2 and radius_on_top<= np.sqrt(n_row*n_row + n_column * n_column)/2:
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
        self.norm_fact = _norm_fact()
        
#   
    def shape(self):
        n_column = self.n_column
        n_row = self.n_row
        n_height = self.n_height
        R_top = self.radius_on_top
#        print R_top
        R_bottom = self.radius_on_bottom
#        print R_bottom
#        r_max = 0
#        r = 0
        for idx_height in range(n_height):
#            print 'r_max = %0.2f' %r_max
            for idx_row in range(n_row):
                for idx_column in range(n_column):
                    Idx_row = np.sqrt(((idx_row+.5) - int(n_row/2.)) * ((idx_row+.5) -int( n_row/2.)) )
                    Idx_column = np.sqrt(((idx_column+.5) - int(n_column/2.)) * ((idx_column+.5) -int( n_column/2.)) )
                    r = np.sqrt(Idx_row*Idx_row + Idx_column*Idx_column) 
#                    print r
                    r_max = R_bottom - (idx_height+1) * 1./n_height * (R_bottom - R_top)
                    if(r> r_max):
#                        print 'r = %0.2f' %r
#                        print 'r_max = %0.2f' %r_max
                        self.lattice[idx_column][idx_row][idx_height] = 0
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
                    Idx_row = np.sqrt(((idx_row+.5) - int(n_row/2.)) * ((idx_row+.5) -int( n_row/2.)) )
                    Idx_column = np.sqrt(((idx_column+.5) - int(n_column/2.)) * ((idx_column+.5) -int( n_column/2.)) )
                    r = np.sqrt(Idx_row*Idx_row + Idx_column*Idx_column) 
                    if r<=self.radius_on_bottom and r<=self.radius_on_top:
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
        #
#        n_height = self.n_height
#        j_row_row = self.j_row_row
#        j_height_height = self.j_height_height
#        j_column_column = self.j_column_column
        #
        dE = -1 * Lattice3D.SpinLattice3D.get_local_energy_change(self, idx_column,idx_row, idx_height)
        #
        Idx_row = np.sqrt(((idx_row+.5) - int(n_row/2.)) * ((idx_row+.5) -int( n_row/2.)) )
        Idx_column = np.sqrt(((idx_column+.5) - int(n_column/2.)) * ((idx_column+.5) -int( n_column/2.)) )
        r = np.sqrt(Idx_row*Idx_row + Idx_column*Idx_column) 
        #
        if (r<=self.radius_on_bottom and r<=self.radius_on_top):
            dE = dE + 2 * self.lattice[idx_column][idx_row][idx_height] * self.external_field_magnitude * self.permitivity
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
   