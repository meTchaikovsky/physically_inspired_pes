
import hashlib
import numpy as np
from numpy.linalg import det


def get_hash(atoms):
    string = str(atoms.pbc)
    for number in atoms.cell.flatten():
        string += '%.15f' % number
    for number in atoms.get_atomic_numbers():
        string += '%3d' % number
    for number in atoms.get_positions().flatten():
        string += '%.15f' % number

    md5 = hashlib.md5(string.encode('utf-8'))
    hash = md5.hexdigest()
    return hash


def get_nodes(nodes, n_processes=28):
    # designed specifically for TH1A
    node_info = {}
    if "[" in nodes:
        nodes = nodes[3:-1]
        if ',' in nodes:
            # 'cn[634,639,793-794]'
            nodes = nodes.split(',')
            for partition in nodes:
                if '-' in partition:
                    nodes = [
                            'cn' + str(item) for item in \
                                range(int(partition.split('-')[0]), int(partition.split('-')[1]) + 1)
                                ]
                    node_info.update({node: n_processes for node in nodes})
                else:
                    node_info.update({'cn' + str(int(partition)): n_processes})
        else:
            nodes = [
                    'cn' + str(item) for item in \
                        range(int(nodes.split('-')[0]), int(nodes.split('-')[1]) + 1)
                        ]
            node_info.update({node: n_processes for node in nodes})
    else:
        node_info = {nodes: n_processes}

    no_of_workers = 0
    for _, ntasks in node_info.items():
        no_of_workers += ntasks
    return node_info


def frac2cart(cell, fracCoords):
    cartCoords = fracCoords @ cell
    return cartCoords


def cart2frac(cell, cartCoords):
    trans_cell = cell.T
    det_trans_cell = det(trans_cell)
    fracCoords = []
    for i in cartCoords:
        aPos, bPos, cPos = trans_cell.copy(), trans_cell.copy(), trans_cell.copy()
        aPos[:, 0] = i; bPos[:, 1] = i; cPos[:, 2] = i
        aPos = 1 / det_trans_cell * det(aPos)
        bPos = 1 / det_trans_cell * det(bPos)
        cPos = 1 / det_trans_cell * det(cPos)
        fracCoords.append([aPos, bPos, cPos])
    return np.array(fracCoords)
        

def traj_creator(cell, name, convert_dict, timestep):
    
    from pymatgen.core.trajectory import Trajectory
    
    with open(name, 'r') as source:
        trajs = source.readlines()
    posts = [ind for ind, item in enumerate(trajs) if 'ITEM: TIMESTEP' in item]
    num_trajs = len(posts)
    
    # take into account the on the fly nature
    try:
        assert len(trajs) // num_trajs == len(trajs) / num_trajs
    except AssertionError:
        num_trajs -= 1
        posts.pop(-1)  
    read_length = len(trajs) // len(posts)
    
    def extractor(post):
        one_stru = trajs[post: post + read_length]
        box = np.array([[float(inner) for inner in item.split()] for item in one_stru[5:8]])
        positions = one_stru[9:]
        inds = [int(item.split()[0]) for item in positions]
        corrected_inds = sorted(range(len(inds)), key=lambda x: inds[x])
        positions = [positions[ind] for ind in corrected_inds]
        
        atomic_species = [int(item.split()[1]) for item in positions]
        atomic_species = np.array([convert_dict[item] for item in atomic_species])
        cart_positions = np.array([[float(inner) for inner in item.split()[2:]] for item in positions])
        frac_positions = cart2frac(cell, cart_positions)

        return [atomic_species, frac_positions, box]
    
    traj_source = [extractor(post) for post in posts]
    species, frac_pos, boxes = zip(*traj_source)
    frac_pos = np.array(frac_pos)[::timestep]
    
    box_initial = boxes[0]
    unique_species = species[0]
    try: assert len(boxes) == sum([np.all(np.abs(box_initial - item) < 1e-5) for item in boxes])
    except AssertionError: print('Simulation box changed during the process!'); return 
    try: assert np.all([unique_species == item for item in species]) == True
    except AssertionError: print('Atomic species changed during the process'); return 
    
    np.save('cell', cell)
    np.save('frac_pos', frac_pos)
    np.save('species', unique_species)
    
    native_trajs = Trajectory(
            lattice=cell, 
            species=unique_species, 
            frac_coords=frac_pos, 
            time_step=timestep
            )        
    
    return native_trajs