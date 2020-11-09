
import sys
import argparse
from os import makedirs
from os import remove

import numpy as np
 
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--db_path', help="The path to the database.", type=str)
    parser.add_argument('--model_path', help="The path to the model.", type=str)
    parser.add_argument('--purpose', help="The purpsoe, train or load.", type=str)
    parser.add_argument('--selected_hashes', help="The file name of the saved selected hashes.", type=str)
    
    args = parser.parse_args()
    db_path = args.db_path
    model_path = args.model_path
    purpose = args.purpose
    selected_hashes = args.selected_hashes
    selected_hashes = np.load(selected_hashes)
    
    # load the fingerprint part 
    afp = tuple(np.load('%s/afp/%s' % (db_path, item)) for item in selected_hashes)
    species_afp = tuple(np.load('%s/species_afp/%s' % (db_path, item)) for item in selected_hashes)
    nAtoms_per_image = np.array([item.shape[0] for item in afp])
    inds_mask = np.concatenate([np.ones(item) * ind for ind, item in enumerate(nAtoms_per_image)])
    
    fp_dim = np.unique([item.shape[1] for item in afp])
    try: assert len(fp_dim) == 1
    except AssertionError: print('Multiple fingerprint dimensions detected! ')
    fp_dim = fp_dim[0]
    
    afp = np.vstack(afp)
    species_afp = np.concatenate(species_afp)
    
    # load the derafp part 
    ind_center = tuple(np.load('%s/ind_center/%s' % (db_path, item)) for item in selected_hashes)
    ind_neighbor = tuple(np.load('%s/ind_neighbor/%s' % (db_path, item)) for item in selected_hashes)
    species_derafp = tuple(np.load('%s/species_derafp/%s' % (db_path, item)) for item in selected_hashes)
    # load the big guy 
    num_neighbors = tuple(item.shape[0] for item in ind_center)
    der_mat = np.empty((sum(num_neighbors), fp_dim, 3)).astype(np.float32)
    offsets = np.cumsum(num_neighbors)
    for ind in range(len(offsets)):
        if ind == 0: 
            container = der_mat[0: offsets[ind]]
        else:
            container = der_mat[offsets[ind - 1]: offsets[ind]]
        container[:] = np.load('%s/der_mat/%s' % (db_path, selected_hashes[ind]))
        
    der_mask = np.concatenate([np.ones(len(item)) * index for index, item in enumerate(ind_center)])
    ind_center = np.concatenate(ind_center)
    ind_neighbor = np.concatenate(ind_neighbor)
    species_derafp = np.concatenate(species_derafp)

    if purpose == 'train': 
        afp_max = np.array([np.max(afp[:, ind]) for ind in range(fp_dim)])
        afp_min = np.array([np.min(afp[:, ind]) for ind in range(fp_dim)])
        diff = afp_max - afp_min
        small_inds = diff < 1e-8
        afp_max[small_inds], afp_min[small_inds] = 1, 0
        np.save('%s/afp_scalar' % model_path, np.array([afp_min, afp_max]))
    elif purpose == 'validation':
        afp_min, afp_max = np.load('%s/afp_scalar.npy' % model_path)
    afp = (afp - afp_min) / (afp_max - afp_min)
    diff = afp_max.reshape(fp_dim, 1) - afp_min.reshape(fp_dim, 1)
    der_mat /= diff
    
    makedirs('tmp_data')
    np.save('tmp_data/afp', afp)
    np.save('tmp_data/species_afp', species_afp)
    np.save('tmp_data/ind_center', ind_center)
    np.save('tmp_data/ind_neighbor', ind_neighbor)
    np.save('tmp_data/species_derafp', species_derafp)
    np.save('tmp_data/der_mat', der_mat)
    np.save('tmp_data/nAtoms_per_image', nAtoms_per_image)
    np.save('tmp_data/inds_mask', inds_mask)
    np.save('tmp_data/der_mask', der_mask)    
    
    remove('./selected.npy')
    
    sys.exit()