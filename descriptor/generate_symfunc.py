import numpy as np

def radii_convertor(radii):
    radii_names = dict(atomic=1, ionic=2, covalent=3, vdw=4, crystal=5)
    return list(map(lambda name: radii_names[name], radii))

def find_pairs(elements):
    found_pairs = []
    for ind, ele1 in enumerate(elements):
        for ele2 in elements[ind: ]:
            found_pairs.append([ele1, ele2])
    return found_pairs

def find_triples(elements):
    found_triples = []
    for ind1, ele1 in enumerate(elements):
        for ind2, ele2 in enumerate(elements[ind1: ]):
            for el3 in elements[ind1 + ind2: ]:
                found_triples.append([ele1, ele2, el3])
    return found_triples

def des_generator(des_dict):
    keys = list(des_dict)
    
    if 'g2' in keys:
        g2_des = [['g2', env, 
                   des_dict['g2']['r_scalars'], 
                   radii_convertor(des_dict['g2']['r_types'])] for env in des_dict['g2']['env']]
    else:
        g2_des = []
        
    if 'g3' in keys:
        g3_des = [['g3', env, 
                   des_dict['g3']['r_scalars'], 
                   radii_convertor(des_dict['g3']['r_types'])] for env in des_dict['g3']['env']]
    else:
        g3_des = []
        
    if 'g3_angle' in keys:
        g3_angle_des = [
                ['g3_angle', gamma, zeta, env, 
                 des_dict['g3_angle']['r_scalars'], 
                 radii_convertor(des_dict['g3_angle']['r_types'])] 
                for zeta in des_dict['g3_angle']['zetas']
                for gamma in des_dict['g3_angle']['gammas'] 
                for env in des_dict['g3_angle']['env']]
    else:
        g3_angle_des = []
        
    if 'g3_angle_partial' in keys:
        g3_angle_partial_des = [
                ['g3_angle_partial', gamma, zeta, env, 
                 des_dict['g3_angle_partial']['r_scalars'], 
                 radii_convertor(des_dict['g3_angle_partial']['r_types'])] 
                for zeta in des_dict['g3_angle_partial']['zetas'] 
                for gamma in des_dict['g3_angle_partial']['gammas']  
                for env in des_dict['g3_angle_partial']['env'] ]
    else:
        g3_angle_partial_des = []
        
    if 'g4' in keys:
        g4_des = [['g4', env, 
                   des_dict['g4']['r_scalars'], 
                   radii_convertor(des_dict['g4']['r_types'])] for env in des_dict['g4']['env']]
    else:
        g4_des = []
    
    des_fp  = [*g2_des, *g3_des, *g3_angle_des, *g3_angle_partial_des, *g4_des]
    lengths = np.sum([len(item[-1]) * len(item[-2]) for item in des_fp])
    
    return des_fp, lengths

#elements = ['Li', 'N']
#des_dict = {
#       'g2': 
#           {'env': [-1], 
#            'r_scalars': [0.5, 0.75, 1.], 
#            'r_types': ['atomic', 'ionic', 'crystal']}, 
#       'g3': 
#           {'env': [[-1, -1]], 
#            'r_scalars': [0.5, 0.75, 1.], 
#            'r_types': ['atomic', 'ionic', 'crystal']}, 
#       'g3_angle':
#           {'env': find_pairs(elements), 
#            'gammas': [1, -1], 
#            'zetas': [1, 4], 
#            'r_scalars': [0.5, 0.75, 1.], 
#            'r_types': ['atomic', 'ionic', 'crystal']}, 
#       'g3_angle_partial':
#           {'env': find_pairs(elements), 
#            'gammas': [1, -1], 
#            'zetas': [1, 4], 
#            'r_scalars': [0.5, 0.75, 1.], 
#            'r_types': ['atomic', 'ionic', 'crystal']}, 
#        'g4': 
#           {'env': [['Li', 'Li', 'Li']], 
#            'r_scalars': [0.5, 0.75, 1.], 
#            'r_types': ['atomic', 'ionic', 'crystal']}
#       }
#
#des_fp, fp_dim = des_generator(des_dict)