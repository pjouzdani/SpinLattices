# -*- coding: utf-8 -*-
import numpy as np
import scipy.linalg as lg
#import matplotlib.pyplot as plt
from random import randint as rand
from random import uniform as randu

import Lattice
#interaction_type ,symmetry,boundary_type,  n_column ,n_row, j_column_column,j_row_row
class SpinLattice2D(Lattice.SpinLattice):
    def  __init__(self, settings_list ):
        self.boundary_type =settings_list[0]['general_settings']['boundary_type']
        self.interaction_type =settings_list[0]['general_settings']['interaction_type']
        self.lattice_geometry = settings_list[0]['general_settings']['lattice_geometry']
        self.quench_disorder = settings_list[0]['general_settings']['quench_disorder']
        self.symmetry = settings_list[0]['general_settings']['symmetry']
        if (self.lattice_geometry == 'square'):
          self.n_column = settings_list[1]['geometry_settings']['n_column']
          self.n_row = settings_list[1]['geometry_settings']['n_row']
          self.j_column_column  = settings_list[1]['geometry_settings']['j_column_column']
          self.j_row_row  = settings_list[1]['geometry_settings']['j_row_row']
          if self.quench_disorder == 'yes':
              self.num_quenches = settings_list[2]['quenches_settings']['num_quenches']
              self.quenches_locations = []
              for item in range(self.num_quenches):
                  # 2D lattice
                  q_row = settings_list[2]['quenches_settings']['quenches_locations'][item][0]
                  q_column = settings_list[2]['quenches_settings']['quenches_locations'][item][1]
                  q_bond = settings_list[2]['quenches_settings']['quenches_locations'][item][2]
                  if q_row <= self.n_row-1 and q_column<=self.n_column and q_bond <=1:                   
                      self.quenches_locations.append(settings_list[2]['quenches_settings']['quenches_locations'][item])
                  else:
                      self.quenches_locations.append([-2,-2,-2])
        #                    
        def num_bonds():
            if self.symmetry == 'Z2':
                if self.interaction_type == 'nn':
                    if self.boundary_type == 'np':
                        if self.lattice_geometry=='square':
                            return (self.n_row)*(self.n_column -1)  + (self.n_column) * (self.n_row-1)
        #
        def build_lattice():
            lattice = []
            if self.lattice_geometry == 'square':
                for idx_row in range(self.n_row):
                    lattice.append(np.ones(self.n_column))
                return lattice
        #
        self.lattice = build_lattice()
        #
        self.num_bonds = num_bonds() 
    #
    #flip_dipole the same
    #
    def get_lattice_energy(self):
        if self.lattice_geometry=='square':
            return Lattice.SpinLattice.get_lattice_energy(self)
    
    def get_local_energy_change(self, idx_column,idx_row):
        if self.quench_disorder == 'yes':
            dE = 0
    #        count = 0
            j_row_row = self.j_row_row
            j_column_column = self.j_column_column      
            n_column = self.n_column
            n_row = self.n_row
            switch1 = 0
            switch1row = 0
            switch1col = 0
            switch2 = 0 
            switch3 = 0         
            for quench in self.quenches_locations:
                #coincides with a disorder
                if idx_row == quench[0] and idx_column == quench[1]:
                    switch1 = 1 
                    if quench[2]==0:
                        if switch1row != 1:
                            switch1row = 1
                    if quench[2]==1:
                        if switch1col !=1:
                            switch1col = 1
                # below to a disorder
                if idx_row-1 == quench[0] and idx_column == quench[1] and quench[2]==0:
                    switch2 = 1 
                # left
                if idx_row == quench[0] and idx_column-1 == quench[1] and quench[2] == 1:
                    switch3 = 1 
            #
            #
            if switch1 ==1:
                if idx_row != n_row-1:
                    if switch1row == 1 :
                        dE = dE + \
                        self.lattice[idx_row][idx_column].item() * \
                        self.lattice[idx_row+1][idx_column].item() *\
                        -1 * j_row_row
                    else:
                        dE = dE + \
                        self.lattice[idx_row][idx_column].item() * \
                        self.lattice[idx_row+1][idx_column].item() *\
                        +1 * j_row_row
                #
                if idx_column != n_column-1:
                    if switch1col == 1 :
                        dE = dE + \
                        self.lattice[idx_row][idx_column].item() * \
                        self.lattice[idx_row][idx_column+1].item() * \
                        -1 * j_column_column
                    else:
                        dE = dE + \
                        self.lattice[idx_row][idx_column].item() * \
                        self.lattice[idx_row][idx_column+1].item() * \
                        +1 * j_column_column
            else:
                if idx_row != n_row-1:
                    dE = dE + \
                    self.lattice[idx_row][idx_column].item() * \
                    self.lattice[idx_row+1][idx_column].item() *\
                    +1 * j_row_row
                if idx_column != n_column-1:
                    dE = dE + \
                    self.lattice[idx_row][idx_column].item() * \
                    self.lattice[idx_row][idx_column+1].item() * \
                    +1 * j_column_column
            #
            if switch2 == 1:
                if idx_row != 0:
                    dE = dE + \
                    self.lattice[idx_row][idx_column].item() * \
                    self.lattice[idx_row-1][idx_column].item() *\
                    -1 * j_row_row
            else:
                dE = dE + \
                self.lattice[idx_row][idx_column].item() * \
                self.lattice[idx_row-1][idx_column].item() *\
                +1 * j_row_row
            #
            if switch3 == 1:
                if idx_column != 0:
                    dE = dE + \
                    self.lattice[idx_row][idx_column].item() * \
                    self.lattice[idx_row][idx_column-1].item() *\
                    -1 * j_column_column
            else:
                if idx_column != 0:
                    dE = dE + \
                    self.lattice[idx_row][idx_column].item() * \
                    self.lattice[idx_row][idx_column-1].item() *\
                    +1 * j_column_column
            return -2*dE
        else:
            return Lattice.SpinLattice.get_local_energy_change(self, idx_column, idx_row)
        
        
        
        
        
        
        
        
        
        
        