# -*- coding: utf-8 -*-
import numpy as np
from random import randint as rand
from random import uniform as randu

def set_impurity( j_field_settings):
#    1)
    tpe = j_field_settings['type']
#    2)
    if tpe == 'random':
        percentage=j_field_settings['percentage']
    n_row = j_field_settings['n_row']
    n_column = j_field_settings['n_column']
    n_height = j_field_settings['n_height']
    magnitude_row = j_field_settings['magnitude_row']
    magnitude_column = j_field_settings['magnitude_column']
    magnitude_height = j_field_settings['magnitude_height']
    random_magnitude_row = j_field_settings['random_magnitude_row']
    random_magnitude_column = j_field_settings['random_magnitude_column']
    random_magnitude_height = j_field_settings['random_magnitude_height']
    j_field_tensor = np.ones((n_column, n_row, n_height, 3))
    for idx_row in range(n_row):
        for idx_column in range(n_column):
            for idx_height in range(n_height):
                # row_row interaction
                j_field_tensor[idx_column][idx_row][idx_height][0] =j_field_tensor[idx_column][idx_row][idx_height][0] * magnitude_row
                # column_column interaction
                j_field_tensor[idx_column][idx_row][idx_height][1] = j_field_tensor[idx_column][idx_row][idx_height][1] * magnitude_column
                # height_height interaction
                j_field_tensor[idx_column][idx_row][idx_height][2] = j_field_tensor[idx_column][idx_row][idx_height][2] * magnitude_height        
    
    for n in range(int(percentage/100. * n_row*n_column*n_height)):
        r_row = rand(0, n_row-1)
#        print 'r_row %d' %r_row
        r_column = rand(0, n_column-1)
#        print 'r_column %d' %r_column
        r_height = rand(0, n_height-1)
#        print 'r_height %d' %r_height
        j_field_tensor[r_column][r_row][r_height][0] =  random_magnitude_row
        j_field_tensor[r_column][r_row][r_height][1] =  random_magnitude_column
        j_field_tensor[r_column][r_row][r_height][2] =  random_magnitude_height
    return j_field_tensor
def set_external_field( E_field_settings):
#    1)
    tpe = E_field_settings['type']
#   dimenssions
    n_row = E_field_settings['n_row']
    n_column = E_field_settings['n_column']
    n_height = E_field_settings['n_height']
#    2)
    magnitude = E_field_settings['magnitude']
    external_field_tensor = np.ones((n_column, n_row, n_height)) * magnitude    
    return external_field_tensor