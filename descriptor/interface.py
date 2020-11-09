
import numpy as np

from ase.io import Trajectory
from ase.neighborlist import NeighborList
from ase.data import atomic_numbers as o_dict
'''Need -1 if no surrounding element is specified.'''
atomic_numbers = o_dict.copy()
atomic_numbers[-1] = -1

from kyzhang.descriptor.f_calc import char_zeta
from kyzhang.utilities import get_hash
        
def calc(
        des_fp, cutoff, 
        ind_center, 
        ind_neighbor, symbol_neighbor, pos_neighbor, 
        ind_neighbor_neighbor, symbols_neighbor_neighbor, pos_neighbor_neighbor, 
        cutofffn_code, p_gamma
        ):
    
    piles = []
    for des in des_fp:
        des_type = des[0]
        if des_type == 'g2':
            des_company, r_scalars, r_types = des[1:]
            piles.append(
                    calculate_g2_prime(
                            r_scalars, 
                            symbol_neighbor, ind_neighbor, pos_neighbor, 
                            symbols_neighbor_neighbor, ind_neighbor_neighbor, pos_neighbor_neighbor, 
                            des_company,  r_types, 
                            cutoff, 
                            ind_center, 
                            cutofffn_code, 
                            p_gamma
                            )
                    )
        if des_type == 'g3':
            des_company, r_scalars, r_types = des[1:]
            piles.append(
                    calculate_g3_prime(
                            r_scalars, 
                            symbol_neighbor, ind_neighbor, pos_neighbor,
                            symbols_neighbor_neighbor, ind_neighbor_neighbor, pos_neighbor_neighbor, 
                            des_company, r_types, 
                            cutoff, 
                            ind_center, 
                            cutofffn_code, 
                            p_gamma
                            )
                    )
        if des_type == 'g3_angle':
            gamma, zeta, des_company, r_scalars, r_types = des[1:]
            piles.append(
                    calculate_g3_angle_prime(
                            r_scalars, 
                            symbol_neighbor, ind_neighbor, pos_neighbor,
                            symbols_neighbor_neighbor, ind_neighbor_neighbor, pos_neighbor_neighbor, 
                            des_company, gamma, zeta,
                            r_types, 
                            cutoff, 
                            ind_center, 
                            cutofffn_code, 
                            p_gamma
                            )
                    )
        if des_type == 'g3_angle_partial':
            gamma, zeta, des_company, r_scalars, r_types = des[1:]
            piles.append(
                    calculate_g3_angle_partial_prime(
                            r_scalars, 
                            symbol_neighbor, ind_neighbor, pos_neighbor,
                            symbols_neighbor_neighbor, ind_neighbor_neighbor, pos_neighbor_neighbor, 
                            des_company, gamma, zeta,
                            r_types, 
                            cutoff, 
                            ind_center, 
                            cutofffn_code, 
                            p_gamma
                            )
                    )
        if des_type == 'g4':
            des_company, r_scalars, r_types = des[1:]
            piles.append(
                    calculate_g4_prime(
                            r_scalars, 
                            symbol_neighbor, ind_neighbor, pos_neighbor, 
                            symbols_neighbor_neighbor, ind_neighbor_neighbor, pos_neighbor_neighbor, 
                            des_company, r_types, 
                            cutoff, 
                            ind_center, 
                            cutofffn_code, 
                            p_gamma
                            )
                    ) 
    der_mat = np.vstack(piles)
    
    return der_mat

def calculate_g2(
        r_scalars, 
        atom_symbol, 
        symbol_neighbors, pos_neighbors,
        des_company, r_types, 
        cutoff, 
        Ri, 
        cutofffn_code, 
        p_gamma=None):
        
    atom_number = atomic_numbers[atom_symbol]
    atomicn_company = atomic_numbers[des_company]
    atomicn_neighbors = [atomic_numbers[symbol] for symbol in symbol_neighbors]

    if len(atomicn_neighbors) == 0:
        return 0.
    else:
        args_calculate_g2 = dict(
                num_rscalars=len(r_scalars), r_scalars=r_scalars, 
                num_n=len(atomicn_neighbors), 
                atomic_num=atom_number, ri=Ri, 
                atomic_num_n=atomicn_neighbors, pos_n=pos_neighbors,
                atomicn_first=atomicn_company, 
                num_rtypes=len(r_types), 
                r_types=r_types, 
                rc=cutoff, cutofffn_code=cutofffn_code
                )
        if p_gamma != None:
            args_calculate_g2['p_gamma'] = p_gamma
            
        ridge = char_zeta.calculate_g2(**args_calculate_g2)
        return ridge

def calculate_g3(
        r_scalars, 
        atom_symbol, 
        symbol_neighbors, pos_neighbors,
        des_company, r_types, 
        cutoff,
        Ri, 
        cutofffn_code,
        p_gamma=None):
    atom_number = atomic_numbers[atom_symbol]
    atomicn_company = sorted([atomic_numbers[item] for item in des_company])
    atomicn_neighbors = [atomic_numbers[symbol] for symbol in symbol_neighbors]
    if len(pos_neighbors) == 0:
        return 0.
    else:
        args_calculate_g3 = dict(
                num_rscalars=len(r_scalars), r_scalars=r_scalars, 
                num_n=len(atomicn_neighbors), 
                atomic_num=atom_number, ri=Ri, 
                atomic_num_n=atomicn_neighbors, pos_n=pos_neighbors, 
                atomicn_first=atomicn_company[0],
                atomicn_second=atomicn_company[1],
                num_rtypes=len(r_types), 
                r_types=r_types, 
                rc=cutoff, cutofffn_code=cutofffn_code
                )
        if p_gamma != None:
            args_calculate_g3['p_gamma'] = p_gamma
            
        ridge = char_zeta.calculate_g3(**args_calculate_g3)
        return ridge
    
def calculate_g3_angle(
        r_scalars, 
        atom_symbol, 
        symbol_neighbors, pos_neighbors,
        des_company, gamma, zeta,
        r_types, 
        cutoff,
        Ri, 
        cutofffn_code,
        p_gamma=None):
    atom_number = atomic_numbers[atom_symbol]
    atomicn_company = sorted([atomic_numbers[item] for item in des_company])
    atomicn_neighbors = [atomic_numbers[symbol] for symbol in symbol_neighbors]
    if len(pos_neighbors) == 0:
        return 0.
    else:
        args_calculate_g3_angle = dict(
                num_rscalars=len(r_scalars), r_scalars=r_scalars, 
                num_n=len(atomicn_neighbors), 
                atomic_num=atom_number, ri=Ri, 
                atomic_num_n=atomicn_neighbors, pos_n=pos_neighbors, 
                atomicn_first=atomicn_company[0],
                atomicn_second=atomicn_company[1],
                g_gamma=gamma, g_zeta=zeta, 
                num_rtypes=len(r_types), 
                r_types=r_types, 
                rc=cutoff, cutofffn_code=cutofffn_code
                )
        if p_gamma != None:
            args_calculate_g3_angle['p_gamma'] = p_gamma
            
        ridge = char_zeta.calculate_g3_angle(**args_calculate_g3_angle)
        return ridge
    
def calculate_g3_angle_partial(
        r_scalars, 
        atom_symbol, 
        symbol_neighbors, pos_neighbors,
        des_company, gamma, zeta,
        r_types, 
        cutoff,
        Ri, 
        cutofffn_code,
        p_gamma=None):
    atom_number = atomic_numbers[atom_symbol]
    atomicn_company = sorted([atomic_numbers[item] for item in des_company])
    atomicn_neighbors = [atomic_numbers[symbol] for symbol in symbol_neighbors]
    if len(pos_neighbors) == 0:
        return 0.
    else:
        args_calculate_g3_angle_partial = dict(
                num_rscalars=len(r_scalars), r_scalars=r_scalars, 
                num_n=len(atomicn_neighbors), 
                atomic_num=atom_number, ri=Ri, 
                atomic_num_n=atomicn_neighbors, pos_n=pos_neighbors, 
                atomicn_first=atomicn_company[0],
                atomicn_second=atomicn_company[1],
                g_gamma=gamma, g_zeta=zeta, 
                num_rtypes=len(r_types), 
                r_types=r_types, 
                rc=cutoff, cutofffn_code=cutofffn_code
                )
        if p_gamma != None:
            args_calculate_g3_angle_partial['p_gamma'] = p_gamma
            
        ridge = char_zeta.calculate_g3_angle_partial(**args_calculate_g3_angle_partial)
        return ridge
    
def calculate_g4(
        r_scalars, 
        atom_symbol, 
        symbol_neighbors, pos_neighbors,
        des_company, r_types, 
        cutoff,
        Ri, 
        cutofffn_code,
        p_gamma=None):
    atom_number = atomic_numbers[atom_symbol]
    atomicn_company = sorted([atomic_numbers[item] for item in des_company])    
    atomicn_neighbors = [atomic_numbers[symbol] for symbol in symbol_neighbors]
    if len(pos_neighbors) == 0:
        return 0.
    else:
        # this dict function is actually special
        args_calculate_g4 = dict(
                num_rscalars=len(r_scalars), r_scalars=r_scalars, 
                num_n=len(atomicn_neighbors), 
                atomic_num=atom_number, ri=Ri, 
                atomic_num_n=atomicn_neighbors, pos_n=pos_neighbors, 
                atomicn_first=atomicn_company[0], 
                atomicn_second=atomicn_company[1], 
                atomicn_third=atomicn_company[2], 
                num_rtypes=len(r_types), 
                r_types=r_types, 
                rc=cutoff, cutofffn_code=cutofffn_code
                )
        if p_gamma != None:
            args_calculate_g4['p_gamma'] = p_gamma
        
        ridge = char_zeta.calculate_g4(**args_calculate_g4)
        return ridge

def calculate_g2_prime(
        r_scalars, 
        symbol_neighbor, ind_neighbor, pos_neighbor, 
        symbol_neighbor_neighbor, ind_neighbor_neighbor, pos_neighbor_neighbor,
        des_company, r_types, 
        cutoff,
        ind_center, 
        cutofffn_code, 
        p_gamma=None):
    atomic_number_n = atomic_numbers[symbol_neighbor]
    atomicn_company = atomic_numbers[des_company]
    atomic_number_n_n = [atomic_numbers[symbol] for symbol in symbol_neighbor_neighbor]
    
    ridge = np.zeros((3, len(r_scalars) * len(r_types))).astype(np.float64)
    if len(atomic_number_n_n) != 0:
        for d in range(3):
            args_calculate_g2_prime=dict(
                    num_rscalars=len(r_scalars), r_scalars=r_scalars, 
                    num_n_n=len(atomic_number_n_n), 
                    ind_n_n=ind_neighbor_neighbor,
                    atomic_num_n_n=atomic_number_n_n,
                    pos_n_n=pos_neighbor_neighbor, 
                    ind_n=ind_neighbor, 
                    atomic_num_n=atomic_number_n, 
                    pos_n=pos_neighbor, 
                    ind_center=ind_center, d=d, 
                    atomicn_first=atomicn_company,
                    num_rtypes=len(r_types), 
                    r_types=r_types, 
                    rc=cutoff, cutofffn_code=cutofffn_code
                    )
            if p_gamma != None:
                args_calculate_g2_prime['p_gamma'] = p_gamma
            
            ridge[d] += char_zeta.calculate_g2_prime(**args_calculate_g2_prime)
    ridge = ridge.T
    return ridge
    
def calculate_g3_prime(
        r_scalars, 
        symbol_neighbor, ind_neighbor, pos_neighbor, 
        symbol_neighbor_neighbor, ind_neighbor_neighbor, pos_neighbor_neighbor,
        des_company, r_types, 
        cutoff,
        ind_center, 
        cutofffn_code, 
        p_gamma=None):
    atomic_number_n = atomic_numbers[symbol_neighbor]
    atomicn_company = sorted([atomic_numbers[item] for item in des_company])
    atomic_number_n_n = [atomic_numbers[symbol] for symbol in symbol_neighbor_neighbor]
    
    ridge = np.zeros((3, len(r_scalars) * len(r_types))).astype(np.float64)
    if len(atomic_number_n_n) != 0:
        for d in range(3):
            args_calculate_g3_prime=dict(
                    num_rscalars=len(r_scalars), r_scalars=r_scalars, 
                    num_n_n=len(atomic_number_n_n),
                    ind_n_n=ind_neighbor_neighbor,
                    atomic_num_n_n=atomic_number_n_n,
                    pos_n_n=pos_neighbor_neighbor,
                    ind_n=ind_neighbor, 
                    atomic_num_n=atomic_number_n, 
                    pos_n=pos_neighbor,
                    ind_center=ind_center, d=d,
                    atomicn_first=atomicn_company[0], 
                    atomicn_second=atomicn_company[1], 
                    num_rtypes=len(r_types), 
                    r_types=r_types, 
                    rc=cutoff, cutofffn_code=cutofffn_code
                    )
            if p_gamma != None:
                args_calculate_g3_prime['p_gamma'] = p_gamma
            
            ridge[d]  += char_zeta.calculate_g3_prime(**args_calculate_g3_prime)
    ridge = ridge.T
    return ridge

def calculate_g3_angle_prime(
        r_scalars, 
        symbol_neighbor, ind_neighbor, pos_neighbor, 
        symbol_neighbor_neighbor, ind_neighbor_neighbor, pos_neighbor_neighbor,
        des_company, gamma, zeta,
        r_types, 
        cutoff,
        ind_center, 
        cutofffn_code, 
        p_gamma=None):
    atomic_number_n = atomic_numbers[symbol_neighbor]
    atomicn_company = sorted([atomic_numbers[item] for item in des_company])
    atomic_number_n_n = [atomic_numbers[symbol] for symbol in symbol_neighbor_neighbor]
    
    ridge = np.zeros((3, len(r_scalars) * len(r_types))).astype(np.float64)
    if len(atomic_number_n_n) != 0:
        for d in range(3):
            args_calculate_g3_angle_prime=dict(
                    num_rscalars=len(r_scalars), r_scalars=r_scalars, 
                    num_n_n=len(atomic_number_n_n),
                    ind_n_n=ind_neighbor_neighbor,
                    atomic_num_n_n=atomic_number_n_n,
                    pos_n_n=pos_neighbor_neighbor,
                    ind_n=ind_neighbor, 
                    atomic_num_n=atomic_number_n, 
                    pos_n=pos_neighbor,
                    ind_center=ind_center, d=d,
                    atomicn_first=atomicn_company[0], 
                    atomicn_second=atomicn_company[1], 
                    g_gamma=gamma, g_zeta=zeta, 
                    num_rtypes=len(r_types), 
                    r_types=r_types, 
                    rc=cutoff, cutofffn_code=cutofffn_code
                    )
            if p_gamma != None:
                args_calculate_g3_angle_prime['p_gamma'] = p_gamma
            
            ridge[d]  += char_zeta.calculate_g3_angle_prime(**args_calculate_g3_angle_prime)
    ridge = ridge.T
    return ridge

def calculate_g3_angle_partial_prime(
        r_scalars, 
        symbol_neighbor, ind_neighbor, pos_neighbor, 
        symbol_neighbor_neighbor, ind_neighbor_neighbor, pos_neighbor_neighbor,
        des_company, gamma, zeta,
        r_types, 
        cutoff,
        ind_center, 
        cutofffn_code, 
        p_gamma=None):
    atomic_number_n = atomic_numbers[symbol_neighbor]
    atomicn_company = sorted([atomic_numbers[item] for item in des_company])
    atomic_number_n_n = [atomic_numbers[symbol] for symbol in symbol_neighbor_neighbor]
    
    ridge = np.zeros((3, len(r_scalars) * len(r_types))).astype(np.float64)
    if len(atomic_number_n_n) != 0:
        for d in range(3):
            args_calculate_g3_angle_partial_prime=dict(
                    num_rscalars=len(r_scalars), r_scalars=r_scalars, 
                    num_n_n=len(atomic_number_n_n),
                    ind_n_n=ind_neighbor_neighbor,
                    atomic_num_n_n=atomic_number_n_n,
                    pos_n_n=pos_neighbor_neighbor,
                    ind_n=ind_neighbor, 
                    atomic_num_n=atomic_number_n, 
                    pos_n=pos_neighbor,
                    ind_center=ind_center, d=d,
                    atomicn_first=atomicn_company[0], 
                    atomicn_second=atomicn_company[1], 
                    g_gamma=gamma, g_zeta=zeta, 
                    num_rtypes=len(r_types), 
                    r_types=r_types, 
                    rc=cutoff, cutofffn_code=cutofffn_code
                    )
            if p_gamma != None:
                args_calculate_g3_angle_partial_prime['p_gamma'] = p_gamma
            
            ridge[d]  += char_zeta.calculate_g3_angle_partial_prime(**args_calculate_g3_angle_partial_prime)
    ridge = ridge.T
    return ridge

def calculate_g4_prime(
        r_scalars, 
        symbol_neighbor, ind_neighbor, pos_neighbor, 
        symbol_neighbor_neighbor, ind_neighbor_neighbor, pos_neighbor_neighbor,
        des_company, r_types, 
        cutoff, 
        ind_center, 
        cutofffn_code, 
        p_gamma=None):
    atomic_number_n = atomic_numbers[symbol_neighbor]
    atomicn_company = sorted([atomic_numbers[item] for item in des_company])  
    atomic_number_n_n = [atomic_numbers[symbol] for symbol in symbol_neighbor_neighbor]
    
    ridge = np.zeros((3, len(r_scalars) * len(r_types))).astype(np.float64)
    if len(atomic_number_n_n) != 0:
        for d in range(3):
            args_calculate_g4_prime=dict(
                    num_rscalars=len(r_scalars), r_scalars=r_scalars, 
                    num_n_n=len(atomic_number_n_n),
                    ind_n_n=ind_neighbor_neighbor,
                    atomic_num_n_n=atomic_number_n_n,
                    pos_n_n=pos_neighbor_neighbor,
                    ind_n=ind_neighbor, 
                    atomic_num_n=atomic_number_n, 
                    pos_n=pos_neighbor,
                    ind_center=ind_center, d=d, 
                    atomicn_first=atomicn_company[0], 
                    atomicn_second=atomicn_company[1], 
                    atomicn_third=atomicn_company[2], 
                    num_rtypes=len(r_types), 
                    r_types=r_types, 
                    rc=cutoff, cutofffn_code=cutofffn_code
                    )
            if p_gamma != None:
                args_calculate_g4_prime['p_gamma'] = p_gamma
                
            ridge[d] += char_zeta.calculate_g4_prime(**args_calculate_g4_prime)
    ridge = ridge.T
    return ridge

def one_image(image):
    
    des_fp = des_dict['atomic_descriptors']
    cutoff = des_dict['cutoff']
    cutofffn_code = des_dict['cutofffn_code']
    p_gamma = des_dict['p_gamma']
    
    n = NeighborList(
            cutoffs=[cutoff / 2.] * len(image), 
            self_interaction=False, 
            bothways=True, 
            skin=0.
            )
    n.update(image)
    neighborlist = [n.get_neighbors(index) for index in range(len(image))]
    
    '''Compute afp'''
    afp_mat, species_afp = [], []
    for ind in range(len(image)):
        ind_neighbors, offset_neighbors = neighborlist[ind]
        symbol_neighbors = [image[_].symbol for _ in ind_neighbors]
        pos_neighbors = [image.positions[neighbor] + np.dot(offset, image.cell) \
                         for (neighbor, offset) in zip(ind_neighbors, offset_neighbors)]
        Ri, atom_symbol = image[ind].position, image[ind].symbol
        
        species_afp.append(atom_symbol)
        
        piles = []
        for des in des_fp:
            des_type = des[0]
            if des_type == 'g2':
                des_company, r_scalars, r_types = des[1:]
                piles.append(
                        calculate_g2(
                                r_scalars, 
                                atom_symbol, 
                                symbol_neighbors, pos_neighbors,
                                des_company, r_types, 
                                cutoff, 
                                Ri, 
                                cutofffn_code, 
                                p_gamma
                                )
                        )
            if des_type == 'g3':
                 des_company, r_scalars, r_types = des[1:]
                 piles.append(
                         calculate_g3(
                                 r_scalars, 
                                 atom_symbol, 
                                 symbol_neighbors, pos_neighbors,
                                 des_company, r_types, 
                                 cutoff, 
                                 Ri, 
                                 cutofffn_code, 
                                 p_gamma
                                 )
                         )
            if des_type == 'g3_angle':
                gamma, zeta, des_company, r_scalars, r_types = des[1:]
                piles.append(
                        calculate_g3_angle(
                                r_scalars, 
                                atom_symbol, 
                                symbol_neighbors, pos_neighbors,
                                des_company, gamma, zeta,
                                r_types, 
                                cutoff, 
                                Ri, 
                                cutofffn_code, 
                                p_gamma
                                )
                        )
            if des_type == 'g3_angle_partial':
                gamma, zeta, des_company, r_scalars, r_types = des[1:]
                piles.append(
                        calculate_g3_angle_partial(
                                r_scalars, 
                                atom_symbol, 
                                symbol_neighbors, pos_neighbors,
                                des_company, gamma, zeta,
                                r_types, 
                                cutoff, 
                                Ri, 
                                cutofffn_code, 
                                p_gamma
                                )
                        )
            if des_type == 'g4':
                des_company, r_scalars, r_types = des[1:]
                piles.append(
                        calculate_g4(
                                r_scalars, 
                                atom_symbol, 
                                symbol_neighbors, pos_neighbors,
                                des_company, r_types,
                                cutoff,
                                Ri, 
                                cutofffn_code, 
                                p_gamma
                                )
                        )
        afp_mat.append(np.hstack((piles)))
    afp_mat = np.asarray(afp_mat)
    species_afp = np.asarray(species_afp)

    '''Compute derivates'''
    prime_tuple = []
    for ind_c in range(len(image)):
        ind_neighbors = neighborlist[ind_c][0]
        symbol_neighbors = [image[ind].symbol for ind in ind_neighbors]
        prime_tuple.append((ind_c, ind_c, image[ind_c].symbol))
        for ind_n, symbol_n in zip(ind_neighbors, symbol_neighbors):
            prime_tuple.append((ind_c, ind_n, symbol_n))
    prime_tuple = list(set(prime_tuple))
    correct_inds = sorted(range(len(prime_tuple)), key = lambda ind: prime_tuple[ind][1])
    prime_tuple = [prime_tuple[ind] for ind in correct_inds]
    
    der_mat_list, ind_center_list, ind_neighbor_list, species_derafp_list = [], [], [], []
    for item in prime_tuple:
        ind_center, ind_neighbor, symbol_neighbor = item
        '''ingrediants needed for computing the afp of the neighbor'''
        pos_neighbor = image.positions[ind_neighbor]
        ind_neighbor_neighbor, offset_neighbor_neighbor = neighborlist[ind_neighbor]
        symbols_neighbor_neighbor = [image[_].symbol for _ in ind_neighbor_neighbor]
        pos_neighbor_neighbor = [image.positions[index] + np.dot(offset, image.get_cell()) 
                              for index, offset in zip(ind_neighbor_neighbor, 
                                                       offset_neighbor_neighbor)]
        der_mat = calc(
                des_fp, cutoff, 
                ind_center, 
                ind_neighbor, symbol_neighbor, pos_neighbor, 
                ind_neighbor_neighbor, symbols_neighbor_neighbor, pos_neighbor_neighbor, 
                cutofffn_code, p_gamma
                )
        der_mat_list.append(der_mat)
        ind_center_list.append(ind_center)
        ind_neighbor_list.append(ind_neighbor)
        species_derafp_list.append(symbol_neighbor)
        
    ind_center_list = np.array(ind_center_list)
    ind_neighbor_list = np.array(ind_neighbor_list)
    species_derafp_list = np.array(species_derafp_list)
    
    np.save('%s/der_mat/%s' % (db_path, get_hash(image)), der_mat_list)
    np.save('%s/ind_center/%s' % (db_path, get_hash(image)), ind_center_list)
    np.save('%s/ind_neighbor/%s' % (db_path, get_hash(image)), ind_neighbor_list)
    np.save('%s/afp/%s' % (db_path, get_hash(image)), afp_mat)
    np.save('%s/species_afp/%s' % (db_path, get_hash(image)), species_afp)
    np.save('%s/species_derafp/%s' % (db_path, get_hash(image)), species_derafp_list)

if __name__ == '__main__':
    
    import sys
    import zmq
    import multiprocessing as mp
    import logging 
    
    print('<worker_connected>')
    try: 
        server_socket = sys.argv[-1]
        port_number = int(server_socket.split(':')[-1])
        
        context = zmq.Context()
        socket = context.socket(zmq.REQ)
        socket.connect('tcp://%s' % server_socket)
        socket.send_pyobj('ready')
        
        current_path, worker_id, db_path, des_dict = socket.recv_pyobj()    
        logging.basicConfig(filename = '%s/worker_log' % current_path) 
        
        images = Trajectory("%s/part_%i.traj" % (current_path, worker_id))
        
        pool = mp.Pool(mp.cpu_count() - 1)
        pool.map(one_image, images)
        pool.close()
        
        socket.send_pyobj('finished')
        socket.recv_pyobj()  
        
    except Exception as e:
        logging.exception('It is about worker %i' % worker_id)
        logging.exception(e)
        socket.send_pyobj('error')
        socket.recv_pyobj()  