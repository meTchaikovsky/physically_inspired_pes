
import sys
import pickle
from os import remove, makedirs, getcwd

import numpy as np

from socket import gethostname
import zmq
import pexpect
from pexpect import pxssh
from getpass import getuser
import multiprocessing as mp

from ase.io import Trajectory
from ase.io import read 
from ase.neighborlist import NeighborList
from ase.data import atomic_numbers as o_dict
'''Need -1 if no surrounding element is specified.'''
atomic_numbers = o_dict.copy()
atomic_numbers[-1] = -1

from . import interface
from ..utilities import get_hash

python = sys.executable

class compute:
    
    def __init__(
            self, 
            elements, 
            des_fp, fp_dim, 
            cutoff, cutofffn_code, p_gamma, 
            ):
        
        self.elements = elements 
        self.des_fp = des_fp
        self.fp_dim = fp_dim
        self.cutoff = cutoff
        self.cutofffn_code = cutofffn_code
        self.p_gamma = p_gamma
        
    def create_db(self, db_path, images, node_info):
        
        try:
            makedirs(db_path)
        except FileExistsError:
            print('The database exists! Try another path!')
            return 
        
        self.db_path = db_path

        makedirs('%s/afp' % db_path)
        makedirs('%s/der_mat' % db_path)
        makedirs('%s/ind_center' % db_path)
        makedirs('%s/ind_neighbor' % db_path)
        makedirs('%s/species_afp' % db_path)
        makedirs('%s/species_derafp' % db_path)
        
        des_dict = {}
        des_dict['elements'] = self.elements
        des_dict['atomic_descriptors'] = self.des_fp
        des_dict['fp_dim'] = self.fp_dim
        des_dict['cutoff'] = self.cutoff
        des_dict['cutofffn_code'] = self.cutofffn_code
        des_dict['p_gamma'] = self.p_gamma
        with open('%s/des_summary' % db_path, 'wb') as source:
            pickle.dump(des_dict, source)
        
        images = [item for item in Trajectory(images)]
        no_of_images = len(images)
        energies = [item.get_potential_energy(apply_constraint = False) for item in images]
        forces = [item.get_forces(apply_constraint = False) for item in images]
        hashes = [get_hash(item) for item in images]
        ene_dict = {hashes[ind]: energies[ind] for ind in range(no_of_images)}
        for_dict = {hashes[ind]: forces[ind] for ind in range(no_of_images)}
        with open('%s/images_info' % db_path, 'wb') as source:
            pickle.dump([ene_dict, for_dict], source)
            
        # a Trajectory object cannot be pickled
        # so to parallelize computing over multiple nodes, you need to divide the whole Trajectory file to pieces evenly
        images_lengths = [len(item) for item in images]
        inds_parts = [
                [ind for ind, item in enumerate(images_lengths) if item == length] 
                    for length in reversed(list(set((images_lengths))))
                ]
        inds_parts_parts = []
        for item in inds_parts:
            N, n = len(item), len(node_info)
            sizes = [N // n if _ >= N % n else N // n + 1 for _ in range(n)]
            parts = [[item.pop(0) for _ in range(size)] for size in sizes]
            inds_parts_parts.append(parts)
        parts = []
        for ind in range(len(node_info)):
            part = []
            for item in inds_parts_parts:
                part.extend(item[ind])
            parts.append(part)
            
        current_path = getcwd()
        for ind, part in enumerate(parts):
            with Trajectory('%s/part_%i.traj' % (current_path, ind), 'w') as source:
                for item in part:
                    source.write(images[item])
            
        context = zmq.Context()
        socket = context.socket(zmq.REP)
        port = socket.bind_to_random_port("tcp://*")
        server_host = gethostname()
        server_socket = "%s:%s" % (server_host, port)
        worker_command = "%s %s " % (python, interface.__file__) + server_socket
        ssh_command = worker_command + ' &'
            
        children = []
        pid_count = 0
        for compute_node in node_info:
            if compute_node == server_host:
                child = pexpect.spawn(worker_command)
                child.expect('<worker_connected>')
                children.append(child)
                pid_count += 1
            else:
                ssh = pxssh.pxssh()
                ssh.login(compute_node, getuser())
                ssh.sendline(ssh_command)
                ssh.expect('<worker_connected>')
                children.append(ssh)
                pid_count += 1
                
        active = pid_count
        worker_ids = list(range(len(node_info)))
        while True:
            message = socket.recv_pyobj()
            if message == 'ready':
                socket.send_pyobj([current_path, worker_ids.pop(0), db_path, des_dict])
            if message == 'finished':
                socket.send_pyobj('good')
                active -= 1
            if message == 'error':
                socket.send_pyobj('bad')
                raise Exception('Something went wrong with one or more workers, please check worker_log') 
            if active == 0:
                break
        # clean up the worker files
        for _ in range(len(node_info)):
            remove('%s/part_%i.traj' % (current_path, _))
        remove('%s/worker_log' % current_path)
        
    def one_structure(self, image):
        
        if isinstance(image, str):
            image = read(image)
            
        self.image = image
        des_fp, cutoff, cutofffn_code, p_gamma  = self.des_fp, self.cutoff, self.cutofffn_code, self.p_gamma
        n = NeighborList(
                cutoffs=[cutoff / 2.] * len(image), 
                self_interaction=False, 
                bothways=True, 
                skin=0.
                )
        n.update(image)
        neighborlist = [n.get_neighbors(index) for index in range(len(image))]
        self.neighborlist = neighborlist
        
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
                            interface.calculate_g2(
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
                            interface.calculate_g3(
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
                            interface.calculate_g3_angle(
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
                            interface.calculate_g3_angle_partial(
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
                            interface.calculate_g4(
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
        
        pool = mp.Pool(mp.cpu_count() - 1)
        der_mat = np.array((pool.map(self.one_tuple, prime_tuple)))
        pool.close()
        
        ind_center, ind_neighbor, species_derafp = zip(*prime_tuple)
        ind_center = np.asarray(ind_center)
        ind_neighbor = np.asarray(ind_neighbor)
        species_derafp = np.asarray(species_derafp)
        
        return [afp_mat, species_afp, ind_center, ind_neighbor, species_derafp, der_mat]
        
    def one_tuple(self, tuple_in):
        
        image = self.image
        neighborlist = self.neighborlist
        des_fp, cutoff = self.des_fp, self.cutoff
        cutofffn_code, p_gamma = self.cutofffn_code, self.p_gamma
        
        ind_center, ind_neighbor, symbol_neighbor = tuple_in
        '''ingrediants needed for computing the afp of the neighbor'''
        pos_neighbor = image.positions[ind_neighbor]
        ind_neighbor_neighbor, offset_neighbor_neighbor = neighborlist[ind_neighbor]
        symbols_neighbor_neighbor = [image[_].symbol for _ in ind_neighbor_neighbor]
        pos_neighbor_neighbor = [image.positions[index] + np.dot(offset, image.get_cell()) 
                              for index, offset in zip(ind_neighbor_neighbor, 
                                                       offset_neighbor_neighbor)]
        der_mat = interface.calc(
                des_fp, cutoff, 
                ind_center, 
                ind_neighbor, symbol_neighbor, pos_neighbor, 
                ind_neighbor_neighbor, symbols_neighbor_neighbor, pos_neighbor_neighbor, 
                cutofffn_code, p_gamma
                )
        
        return der_mat