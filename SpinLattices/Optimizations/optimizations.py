# -*- coding: utf-8 -*-
import numpy as np
import scipy.linalg as lg
#import matplotlib.pyplot as plt
from random import randint as rand
from random import uniform as randu

class Optimizers:
    #    
    def make_moves_v1(self,  lattice, n_warmup, temperature, moves, lattice_enery):                       
        n_column = lattice.n_column
        n_row = lattice.n_row
        #prob_record = []
        for s in range(n_warmup):
            idx_column=rand(0,n_column-1)
            idx_row=rand(0,n_row-1)
            dE = lattice.get_local_energy_change(idx_column,idx_row)                    
            probability_ratio = np.exp(-dE/(1.0*temperature))
            #prob_record.append(probability_ratio)
            # if energetically preferred             
            if probability_ratio > 1:
                lattice.flip_dipole(idx_column,idx_row)
                moves.append(1)
                lattice_enery = lattice_enery + dE
            # if not energetically preferred, by random...
            if probability_ratio <= 1 and probability_ratio >= 0:
                r = randu(0,1)
                if r< probability_ratio: 
                    lattice.flip_dipole(idx_column,idx_row)
                    moves.append(1)
                    lattice_enery = lattice_enery + dE
                else:
                    moves.append(0)
            if probability_ratio < 0:
                moves.append(0)
        return moves, lattice_enery
    #
    # minimizes with respect to the energy 
    def optimize_spinlattice_at_energy_warmup(self, lattice, temperature):
        _max_num_warmup_loops = 200
        _tolerance = 0.01
        _num_warmp = 30
        _averaging_intervales = 10
        # Constructor
        # energy init
        energy_array = []
        energy_mean_array=[]
        energy_var_array = []
        energy_2power2_mean_array = []
        #observables init
        magnetization_array = []
        magnetization_mean_array = []
        # init energy
        # normalization factor
        energy_normalization_factor=lattice.num_bonds
        #total energy of the lattice
        lattice_energy = lattice.get_lattice_energy()
        # save energy
        energy_array.append(lattice_energy/energy_normalization_factor)
        # mean of the energy up to the num_warmup index
        energy_mean_array.append(np.mean(energy_array))
        # E^2 
        energy_2power2_array=self.get_squre_list(energy_array)
        # mean of E^2 : <E^2>
        energy_2power2_mean_array.append(np.mean(energy_2power2_array))
        # < E^2 > - <E>^2
        energy_var_array.append( energy_2power2_mean_array[-1] - \
                        energy_mean_array[-1]*energy_mean_array[-1] )
        
        #
        #    
        #magnetic
        magnetization_norm_fact = lattice.n_column * lattice.n_row * 1.
        magnetization_array.append(sum(sum(lattice.lattice))/magnetization_norm_fact)
        magnetization_mean_array.append(np.mean(magnetization_array))
        #
        #
        error = []
        # init
        num_warmup = 1 
        error.append(1)
        # end of constructor
        #
        while(error[num_warmup-1] >= _tolerance and num_warmup < _max_num_warmup_loops):
            if num_warmup <= _averaging_intervales:
                moves, lattice_energy = self.make_moves_v1(lattice, _num_warmp,
                                                      temperature, 
                                                      [], lattice_energy)
                energy_array.append(lattice_energy/energy_normalization_factor)
                energy_mean_array.append(np.mean(energy_array))
                energy_2power2_array=self.get_squre_list(energy_array)
                energy_2power2_mean_array.append(np.mean(energy_2power2_array))
                energy_var_array.append( energy_2power2_mean_array[-1] - \
                        energy_mean_array[-1]*energy_mean_array[-1] )
                #
                magnetization_array.append(sum(sum(lattice.lattice))/magnetization_norm_fact)
                magnetization_mean_array.append(np.mean(magnetization_array))
                #
                #print num_warmup
                num_warmup = num_warmup + 1
                error.append(1)
            #                
            if num_warmup > _averaging_intervales:
                moves, lattice_energy = self.make_moves_v1(lattice, _num_warmp,
                                                      temperature, 
                                                      [], lattice_energy)
                
                energy_array.append(lattice_energy/energy_normalization_factor)
                energy_mean_array.append(np.mean(energy_array))
                energy_2power2_array=self.get_squre_list(energy_array)
                energy_2power2_mean_array.append(np.mean(energy_2power2_array))
                energy_var_array.append( energy_2power2_mean_array[-1] - \
                        energy_mean_array[-1]*energy_mean_array[-1] )
                #
                magnetization_array.append(sum(sum(lattice.lattice))/magnetization_norm_fact)
                magnetization_mean_array.append(np.mean(magnetization_array))
                #
                error.append( np.abs(np.mean(energy_mean_array[-_averaging_intervales/2:-1])\
                        - np.mean(energy_mean_array[-_averaging_intervales:-_averaging_intervales/2])) )
                #print num_warmup
                num_warmup = num_warmup +1
                #print error
                
        #
        return energy_mean_array, energy_var_array, magnetization_mean_array
    #                
    #
    def optimize_spinlattice_at_energy_measurement(self,  lattice, temperature, num_sampling, measurement_settings):
        energy_mean_array, energy_var_array, magnetization_mean_array = self.optimize_spinlattice_at_energy_warmup(lattice, temperature)
        # constructor
        lattice_avrg = []
        for idx_row in range(lattice.n_row):
             lattice_avrg.append(np.zeros(lattice.n_column))
        energy_avrg = 0
        energy_var_avg = 0
        polarization_avrg = 0
        #
        for n in range(num_sampling):
            energy_mean_array, energy_var_array,magnetization_mean_array = self.optimize_spinlattice_at_energy_warmup(lattice, temperature)
            #
            energy_avrg = energy_avrg + energy_mean_array[-1]/num_sampling
            polarization_avrg = polarization_avrg+ magnetization_mean_array[-1]/num_sampling
            energy_var_avg   = energy_var_avg + energy_var_array[-1]/num_sampling
            # !!! ENCAPSULATION VIOLATION !!!! ##############################################
            for idx_column in range(lattice.n_column):
                for idx_row in range(lattice.n_row):
                    lattice_avrg[idx_row][idx_column]  = lattice_avrg[idx_row][idx_column] + lattice.lattice[idx_row][idx_column]/num_sampling
            #!!! ENCAPSULATION VIOLATION !!!! ##############################################
        #
        return lattice_avrg, energy_avrg, polarization_avrg, energy_var_avg
    #                
    #
    ##########  To Be transfered
    def get_squre_list(self, list):
        rslt=[]
        for i in range(len(list)):
            rslt.append(list[i]*list[i])
        return rslt
    #
    #
class Optimizers3D:
    def make_moves_v13D(self, lattice, n_warmup, temperature, moves, lattice_enery):                       
        n_column = lattice.n_column
        n_row = lattice.n_row
        n_height = lattice.n_height
        #prob_record = []
        for s in range(n_warmup):
            idx_column=rand(0,n_column-1)
            idx_row=rand(0,n_row-1)
            idx_height=rand(0,n_height-1)
            dE = lattice.get_local_energy_change(idx_column,idx_row, idx_height)                    
            probability_ratio = np.exp(-dE/(1.0*temperature))
            #prob_record.append(probability_ratio)
            # if energetically preferred             
            if probability_ratio > 1:
                lattice.flip_dipole(idx_column,idx_row, idx_height)
                moves.append(1)
                lattice_enery = lattice_enery + dE
            # if not energetically preferred, by random...
            if probability_ratio <= 1 and probability_ratio >= 0:
                r = randu(0,1)
                if r< probability_ratio: 
                    lattice.flip_dipole(idx_column,idx_row, idx_height)
                    moves.append(1)
                    lattice_enery = lattice_enery + dE
                else:
                    moves.append(0)
            if probability_ratio < 0:
                moves.append(0)
        return moves, lattice_enery
    
    # minimizes with respect to the energy 
    def optimize_spinlattice_at_energy_warmup3D(self, lattice, temperature):
        _max_num_warmup_loops = 200
        _tolerance = 0.01
        _num_warmp = 30
        _averaging_intervales = 10
        # Constructor
        # energy init
        energy_array = []
        energy_mean_array=[]
        energy_var_array = []
        energy_2power2_mean_array = []
        #observables init
        magnetization_array = []
        magnetization_mean_array = []
        # init energy
        # normalization factor
        energy_normalization_factor=lattice.num_bonds
        #total energy of the lattice
        lattice_energy = lattice.get_lattice_energy()
        # save energy
        energy_array.append(lattice_energy/energy_normalization_factor)
        # mean of the energy up to the num_warmup index
        energy_mean_array.append(np.mean(energy_array))
        # E^2 
        energy_2power2_array= self.get_squre_list(energy_array)
        # mean of E^2 : <E^2>
        energy_2power2_mean_array.append(np.mean(energy_2power2_array))
        # < E^2 > - <E>^2
        energy_var_array.append( energy_2power2_mean_array[-1] - \
                        energy_mean_array[-1]*energy_mean_array[-1] )
        
        #
        #    
        #magnetic
        magnetization_norm_fact = lattice.n_column * lattice.n_row *  lattice.n_height * 1.
        # contraction onto a z-dimension
        aux = []
        for idx_height in range(lattice.n_height):
            aux.append(sum(sum(lattice.lattice[idx_height])))
        magnetization_array.append(sum(aux)/magnetization_norm_fact)
        magnetization_mean_array.append(np.mean(magnetization_array))
        #
        #
        error = []
        # init
        num_warmup = 1 
        error.append(1)
        # end of constructor
        #
        while(error[num_warmup-1] >= _tolerance and num_warmup < _max_num_warmup_loops):
            if num_warmup <= _averaging_intervales:
                moves, lattice_energy = self.make_moves_v13D(lattice, _num_warmp,
                                                      temperature, 
                                                      [], lattice_energy)
                energy_array.append(lattice_energy/energy_normalization_factor)
                energy_mean_array.append(np.mean(energy_array))
                energy_2power2_array= self.get_squre_list(energy_array)
                energy_2power2_mean_array.append(np.mean(energy_2power2_array))
                energy_var_array.append( energy_2power2_mean_array[-1] - \
                        energy_mean_array[-1]*energy_mean_array[-1] )
                #
                aux = []
                for idx_height in range(lattice.n_height):
                    aux.append(sum(sum(lattice.lattice[idx_height])))
                magnetization_array.append(sum(aux)/magnetization_norm_fact)
                magnetization_mean_array.append(np.mean(magnetization_array))
                #
                #print num_warmup
                num_warmup = num_warmup + 1
                error.append(1)
            #                
            if num_warmup > _averaging_intervales:
                moves, lattice_energy = self.make_moves_v13D(lattice, _num_warmp,
                                                      temperature, 
                                                      [], lattice_energy)
                
                energy_array.append(lattice_energy/energy_normalization_factor)
                energy_mean_array.append(np.mean(energy_array))
                energy_2power2_array=self.get_squre_list(energy_array)
                energy_2power2_mean_array.append(np.mean(energy_2power2_array))
                energy_var_array.append( energy_2power2_mean_array[-1] - \
                        energy_mean_array[-1]*energy_mean_array[-1] )
                #
                aux = []
                for idx_height in range(lattice.n_height):
                    aux.append(sum(sum(lattice.lattice[idx_height])))
                magnetization_array.append(sum(aux)/magnetization_norm_fact)
                magnetization_mean_array.append(np.mean(magnetization_array))
                #
                error.append( np.abs(np.mean(energy_mean_array[-_averaging_intervales/2:-1])\
                        - np.mean(energy_mean_array[-_averaging_intervales:-_averaging_intervales/2])) )
                #print num_warmup
                num_warmup = num_warmup +1
                #print error
                
        #
        return energy_mean_array, energy_var_array, magnetization_mean_array
    #            
    #
    def optimize_spinlattice_at_energy_measurement3D(self,  lattice, temperature, num_sampling, measurement_settings):
        energy_mean_array, energy_var_array, magnetization_mean_array = self.optimize_spinlattice_at_energy_warmup3D(lattice, temperature)
        # constructor
        lattice_avrg = []
        lattice_avrg_2D = []
        for idx_column in range(lattice.n_column):
            lattice_avrg_2D.append(np.zeros(lattice.n_height))
        for idx_row in range(lattice.n_row):
             lattice_avrg.append(lattice_avrg_2D)
        energy_avrg = 0
        energy_var_avg = 0
        polarization_avrg = 0
        #
        for n in range(num_sampling):
            energy_mean_array, energy_var_array,magnetization_mean_array = self.optimize_spinlattice_at_energy_warmup3D(lattice, temperature)
            #
            energy_avrg = energy_avrg + energy_mean_array[-1]/num_sampling
            polarization_avrg = polarization_avrg+ magnetization_mean_array[-1]/num_sampling
            energy_var_avg   = energy_var_avg + energy_var_array[-1]/num_sampling
            # !!! ENCAPSULATION VIOLATION !!!! ##############################################
            for idx_column in range(lattice.n_column):
                for idx_row in range(lattice.n_row):
                    for idx_height in range(lattice.n_height):                
                        lattice_avrg[idx_row][idx_column][idx_height]  = lattice_avrg[idx_row][idx_column][idx_height] +\
                        lattice.lattice[idx_row][idx_column][idx_height]/num_sampling
            #!!! ENCAPSULATION VIOLATION !!!! ##############################################
        #
        return lattice_avrg, energy_avrg, polarization_avrg, energy_var_avg
    #                
    #
    def get_squre_list(self, list):
        rslt=[]
        for i in range(len(list)):
            rslt.append(list[i]*list[i])
        return rslt

    
    
    