import sys

from kyzhang.descriptor.gaussian import compute
from kyzhang.descriptor.generate_symfunc import find_pairs
from kyzhang.descriptor.generate_symfunc import find_triples
from kyzhang.descriptor.generate_symfunc import des_generator

from kyzhang.utilities import get_nodes
#nodes = sys.argv[-1]
#node_info = get_nodes(nodes)

node_info = {'th-hpc2-ln0': 10}

images = '/THL7/home/lcyin/KYZhang/bp_nn_pes/solid_electrolytes/alpha-Li3N/gamma_only/r_scalars-1-5-10-no_g4/run_lmp/nvt_runs/bad_trajs/bad_strus.traj'
db_path = '/THL7/home/lcyin/KYZhang/bp_nn_pes/playground/create_db/database'
elements = ['Li', 'N']
des_dict = {
       'g2': 
           {'env': elements, 
            'r_scalars': [1, 5, 10],
            'r_types': ['atomic', 'ionic']}, 
       'g3': 
           {'env': [[-1, -1]], 
            'r_scalars': [1, 5, 10],
            'r_types': ['atomic', 'ionic']}, 
       'g3_angle':
           {'env': find_pairs(elements), 
            'gammas': [1, -1], 
            'zetas': [1, 4], 
            'r_scalars': [1, 5, 10],
            'r_types': ['atomic', 'ionic']}, 
       'g3_angle_partial':
           {'env': find_pairs(elements), 
            'gammas': [1, -1], 
            'zetas': [1, 4], 
            'r_scalars': [1, 5, 10], 
            'r_types': ['atomic', 'ionic']}, 
       'g4':
           {'env': [['Li', 'Li', 'Li']],
            'r_scalars': [1, 5, 10],
            'r_types': ['atomic', 'ionic']}
       }
cutoff = 4.5
cutofffn_code = 2
p_gamma = 4
atomic_descriptors, fp_dim = des_generator(des_dict)
worker = compute(
        elements=elements, 
        des_fp=atomic_descriptors, fp_dim=fp_dim,
        cutoff=cutoff, cutofffn_code=cutofffn_code, p_gamma=p_gamma)
worker.create_db(db_path, images, node_info)

