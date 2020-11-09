
import numpy as np
import pickle

from ..train.model import neural_network

from ase.data import atomic_numbers as o_dict
atomic_numbers = o_dict.copy()
atomic_numbers[-1] = -1


def model_extractor(model_paths, des_path):
        
    with open('%s/des_summary' % des_path, 'rb') as source:
        des_dict = pickle.load(source)
    cutoff = des_dict['cutoff']
    cutofffn_code = des_dict['cutofffn_code']
    p_gamma = des_dict['p_gamma']
    elements = des_dict['elements']
    atomic_descriptors = des_dict['atomic_descriptors']
    fp_dim = des_dict['fp_dim']
        
    with open('model.param', 'w') as source:
        source.write('%i\n' % len(elements))
        source.write('%s\n' % ' '.join(elements))
        source.write('%s\n' % ' '.join([str(atomic_numbers[ele]) for ele in elements]))
        source.write('%i\n' % len(atomic_descriptors))
        
        for ind, des in enumerate(atomic_descriptors):
            des_type = des[0]
            source.write('%s\n' % des_type)
            if des_type == 'g2':
                atomic_num = atomic_numbers[des[1]]
                source.write('%i\n' % atomic_num)
                # r_scalars
                r_scalars = des[2].copy()
                source.write('%i\n' % len(r_scalars))
                r_scalars = [str(item) for item in r_scalars]
                source.write('%s\n' % ' '.join(r_scalars))
                # r_types
                r_types = des[3].copy()
                source.write('%i\n' % len(r_types))
                r_types = [str(item) for item in r_types]
                source.write('%s\n' % ' '.join(r_types))
                
            if des_type == 'g3':
                atomicn_one = atomic_numbers[des[1][0]]
                atomicn_two = atomic_numbers[des[1][1]]
                source.write('%i %i\n' % (atomicn_one, atomicn_two))
                # r_scalars
                r_scalars = des[2].copy()
                source.write('%i\n' % len(r_scalars))
                r_scalars = [str(item) for item in r_scalars]
                source.write('%s\n' % ' '.join(r_scalars))
                # r_types
                r_types = des[3].copy()
                source.write('%i\n' % len(r_types))
                r_types = [str(item) for item in r_types]
                source.write('%s\n' % ' '.join(r_types))
                
            if des_type == 'g3_angle':
                atomicn_one = atomic_numbers[des[3][0]]
                atomicn_two = atomic_numbers[des[3][1]]
                source.write('%i %i\n' % (atomicn_one, atomicn_two))
                source.write('%i %i\n' % (des[1], des[2]))
                # r_scalars
                r_scalars = des[4].copy()
                source.write('%i\n' % len(r_scalars))
                r_scalars = [str(item) for item in r_scalars]
                source.write('%s\n' % ' '.join(r_scalars))
                # r_types
                r_types = des[5].copy()
                source.write('%i\n' % len(r_types))
                r_types = [str(item) for item in r_types]
                source.write('%s\n' % ' '.join(r_types))
            
            if des_type == 'g3_angle_partial':
                atomicn_one = atomic_numbers[des[3][0]]
                atomicn_two = atomic_numbers[des[3][1]]
                source.write('%i %i\n' % (atomicn_one, atomicn_two))
                source.write('%i %i\n' % (des[1], des[2]))
                # r_scalars
                r_scalars = des[4].copy()
                source.write('%i\n' % len(r_scalars))
                r_scalars = [str(item) for item in r_scalars]
                source.write('%s\n' % ' '.join(r_scalars))
                # r_types
                r_types = des[5].copy()
                source.write('%i\n' % len(r_types))
                r_types = [str(item) for item in r_types]
                source.write('%s\n' % ' '.join(r_types))
                
            if des_type == 'g4':
                atomicn_one = atomic_numbers[des[1][0]]
                atomicn_two = atomic_numbers[des[1][1]]
                atomicn_three = atomic_numbers[des[1][2]]
                source.write('%i %i %i\n' % (atomicn_one, atomicn_two, atomicn_three))
                # r_scalars
                r_scalars = des[2].copy()
                source.write('%i\n' % len(r_scalars))
                r_scalars = [str(item) for item in r_scalars]
                source.write('%s\n' % ' '.join(r_scalars))
                # r_types
                r_types = des[3].copy()
                source.write('%i\n' % len(r_types))
                r_types = [str(item) for item in r_types]
                source.write('%s\n' % ' '.join(r_types))
                            
        source.write('%f\n' % cutoff)
        # 1 for cosine, 2 for polynomial
        source.write('%i\n' % cutofffn_code)
        source.write('%i\n' % p_gamma)
        
        source.write('%i\n' % len(model_paths))
        
        for model_path in model_paths:
            
            afp_min, afp_max = np.load('%s/afp_scalar.npy' % model_path)
            with open('%s/model_config' % model_path, 'rb') as model_source:
                model_dict = pickle.load(model_source)
            activation = model_dict['activation']
            species_partitions = model_dict['species_partitions']
            hidden_layers = model_dict['hidden_layers']
                
            nn_model = neural_network(
                    species_partitions=species_partitions,
                    fp_dim=fp_dim,
                    hidden_layers=hidden_layers,
                    activation=activation,
                    )
            flattened_parameters = nn_model.model_extractor(model_path)
            
            atomicn_partitions = [[atomic_numbers[inner] for inner in item] for item in species_partitions]
            navigator = []
            for ind, item in enumerate(atomicn_partitions):
                temp_nav = [(item_inner, ind + 1) for item_inner in item]
                navigator.extend(temp_nav)
            
            # write the models 
            source.write('%i %i\n' % (len(species_partitions), max([atomic_numbers[ele] for ele in elements])))
            for item in navigator: source.write('%i %i\n' % item)
            
            for nn_ind, nn_config in enumerate(hidden_layers):
                all_layers = [fp_dim, *nn_config, 1]
                source.write('%i\n' % len(all_layers))
                no_of_nodes_layer = ' '.join([str(item) for item in all_layers])
                source.write('%s\n' % no_of_nodes_layer)
                source.write('%s\n' % ' '.join(str(_) for _ in flattened_parameters[nn_ind]))
            # write the min and max 
            source.write('%s\n' % ' '.join(str(_) for _ in afp_min))
            source.write('%s\n' % ' '.join(str(_) for _ in afp_max))