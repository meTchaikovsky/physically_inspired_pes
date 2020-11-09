
import sys
import pickle 

from compare.gaussian import compute

from kyzhang.utilities import get_nodes
nodes = sys.argv[-1]
node_info = get_nodes(nodes)

images = 'ligeps-every_50-high_accuracy.traj'
db_path = '/THL7/home/lcyin/KYZhang/bp_nn_pes/solid_electrolytes/Li10GeP2S12/r_scalars-5_20/cutoff-4.5/original_acsfs/database'

elements = ['Li', 'Ge', 'P', 'S']
cutoff = 4.5
cutofffn_code = 2
p_gamma = 4

with open('descriptors','rb') as source:
    atomic_descriptors, fp_dim = pickle.load(source)
    
worker = compute(
        elements=elements,
        des_fp=atomic_descriptors, fp_dim=fp_dim,
        cutoff=cutoff, cutofffn_code=cutofffn_code, p_gamma=p_gamma)
worker.create_db(db_path, images, node_info)
                                                 