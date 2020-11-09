
import sys 
import pickle
import os 
from os import system
from os import listdir
from os import makedirs
from os.path import isdir
from shutil import rmtree

import numpy as np

try:
    from sklearn.model_selection import train_test_split
except ImportError:
    from sklearn.cross_validation import train_test_split
import tensorflow as tf
from tensorflow.contrib.opt import ScipyOptimizerInterface

from . import loader

from ..descriptor.gaussian import compute

from .tools import Convergence_Occurred, Early_stopping
from .tools import weight_bias_initializer
from .tools import check_convergence

# supress tensorflow warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
python = sys.executable

class neural_network:

    def __init__(
            self,
            species_partitions,
            fp_dim, 
            hidden_layers, 
            activation
            ):
        
        self.species_partitions = species_partitions
        self.fp_dim = fp_dim
        self.hidden_layers = hidden_layers
        self.activation = activation

    def engine(self, des_path, saved_model, image):
    
        with open('%s/des_summary' % des_path, 'rb') as source:
            des_dict = pickle.load(source)
        worker = compute(
            elements=des_dict['elements'],
            des_fp=des_dict['atomic_descriptors'], fp_dim=des_dict['fp_dim'],
            cutoff=des_dict['cutoff'],  cutofffn_code=des_dict['cutofffn_code'], p_gamma=des_dict['p_gamma']
            )
        afp, species_afp, ind_center, ind_neighbor, species_derafp, der_mat = worker.one_structure(image)
        
        afp_min, afp_max = np.load('%s/afp_scalar.npy' % saved_model)
        inds_mask = np.zeros(afp.shape[0])
        afp = (afp - afp_min) / (afp_max - afp_min)
        afp_min, afp_max = afp_min.reshape(-1, 1), afp_max.reshape(-1, 1)
        der_mat /= (afp_max - afp_min)

        not_confident = np.logical_and(np.all(afp < 0, axis = 1), np.all(afp > 1, axis = 1))
        wrong_inds = [ind for ind in range(len(not_confident)) if not_confident[ind]]
        if not len(wrong_inds) == 0: print("The computed fingerprints of atom(s) %s fall outside the confidence interval! " % (', '.join(wrong_inds)))
        
        self.weights = [[] for id_s, item_s in enumerate(self.species_partitions)]
        self.bias = [[] for id_s, item_s in enumerate(self.species_partitions)]
        self.construct_graph(xavier_initializer=True)
        with self.graph.as_default():
            sess = tf.Session()
            saver = tf.train.Saver()
            sess.run(tf.global_variables_initializer())
            saver.restore(sess, '%s/model.ckpt' % saved_model)
            
            feed_dict = {}
            for ind_part, ele_names in enumerate(self.species_partitions):
                afp_slice = np.in1d(species_afp, ele_names)
                derafp_slice = np.in1d(species_derafp, ele_names)
    
                feed_dict[self.fingerprints_ph[ind_part]] = afp[afp_slice]
                feed_dict[self.inds_mask_ph[ind_part]] = inds_mask[afp_slice]
    
                feed_dict[self.der_mat_ph[ind_part]] = der_mat[derafp_slice]
                feed_dict[self.ind_center_ph[ind_part]] = ind_center[derafp_slice]
    
                neighbor_slice = ind_neighbor[derafp_slice]
                convert_dict = {item: ind for ind, item in enumerate(np.unique(neighbor_slice))}
                neighbor_slice = np.array([convert_dict[item] for item in neighbor_slice])
                feed_dict[self.ind_neighbor_ph[ind_part]] = neighbor_slice
    
            feed_dict[self.totalNumAtoms_ph] = afp.shape[0]
            feed_dict[self.rate_ph] = 0
            feed_dict[self.batchsize_ph] = 1
            
            energy = sess.run(self.energies, feed_dict=feed_dict)[0]
            forces = sess.run(self.forces, feed_dict=feed_dict)

        return energy, forces

    def model_extractor(self, saved_model):

        self.weights = [[] for id_s, item_s in enumerate(self.species_partitions)]
        self.bias = [[] for id_s, item_s in enumerate(self.species_partitions)]
        self.construct_graph(xavier_initializer=True)
        with self.graph.as_default():
            sess = tf.Session()
            saver = tf.train.Saver()
            # it is essential to have this line above the line for restoring
            sess.run(tf.global_variables_initializer())
            saver.restore(sess, '%s/model.ckpt' % saved_model)
            
            flattened_parameters = []
            for ind in range(len(self.species_partitions)):
                model_param = []
                for _ in range(len(self.weights[ind])):
                    temp_weight = sess.run(self.weights[ind][_])
                    temp_bias = sess.run(self.bias[ind][_])
                    temp_weight = temp_weight.tolist()
                    temp_bias = temp_bias.tolist()
                    temp_weight.append(temp_bias)
                    model_param.extend(np.array(temp_weight).ravel())
                flattened_parameters.append(model_param)

        return flattened_parameters

    def generate_feeddict(self, db_path, selected_hashes, train_params, purpose):
        
        loader_file = loader.__file__
        np.save('./selected.npy', selected_hashes)
        system(
                '%s %s --db_path=%s --model_path=%s --purpose=%s --selected_hashes=%s' % 
                (python, loader_file, db_path, self.model_path, purpose, './selected.npy')
               )
        nAtoms_per_image = np.load('tmp_data/nAtoms_per_image.npy')
        inds_mask = np.load('tmp_data/inds_mask.npy')
        der_mask = np.load('tmp_data/der_mask.npy')
        afp = np.load('tmp_data/afp.npy')
        species_afp = np.load('tmp_data/species_afp.npy')
        species_derafp = np.load('tmp_data/species_derafp.npy')
        ind_neighbor = np.load('tmp_data/ind_neighbor.npy')
        ind_center = np.load('tmp_data/ind_center.npy')
        
        with open('%s/images_info' % db_path, 'rb') as source:
            all_ene, all_forces = pickle.load(source)
        energies = np.array([all_ene[item.split('.')[0]] for item in selected_hashes])
        forces = np.array([all_forces[item.split('.')[0]] for item in selected_hashes])
        
        # generate feeddict
        feed_dict = {}
        feed_dict[self.force_coefficient_ph] = train_params['force_coefficient']
        feed_dict[self.energy_coefficient_ph] = train_params['energy_coefficient']
        feed_dict[self.regularization_strength_ph] = train_params['regularization_strength']
        
        feed_dict[self.energies_in_ph] = energies.reshape(-1, 1)
        feed_dict[self.forces_in_ph] = np.concatenate(forces, axis = 0)
        feed_dict[self.totalNumAtoms_ph] = np.sum(nAtoms_per_image)
        feed_dict[self.natoms_per_image_ph] = nAtoms_per_image.reshape(-1, 1)
        
        # for grabing the force that is the hardest to fit 
        feed_dict[self.forces_mask] = np.concatenate(
                [np.ones(item.shape[0]) * ind for ind, item in enumerate(forces)]
                )
        feed_dict[self.forces_atoms_mask] = np.concatenate([np.arange(item.shape[0]) for item in forces])
        feed_dict[self.selected_hashes] = np.array(selected_hashes)
        
        natoms_offset = np.cumsum(nAtoms_per_image) - nAtoms_per_image
        for ind_part, ele_names in enumerate(self.species_partitions):
            afp_validkeys = np.in1d(species_afp, ele_names)
            der_validkeys = np.in1d(species_derafp, ele_names)    
    
            fingerprints_slice = afp[afp_validkeys]
            # understand why you don't need to reindex inds_mask from how tf.unsorted_segment_sum()
            inds_mask_slice = inds_mask[afp_validkeys]
            ind_neighbor_slice = ind_neighbor[der_validkeys]
            ind_center_slice = ind_center[der_validkeys]
            der_mask_slice = der_mask[der_validkeys]
                
            # reindexing ind_neighbor so it matches with the sliced fingerprints
            for index in range(len(selected_hashes)):
                der_slice = ind_neighbor_slice[der_mask_slice == index]
                unique_items = np.unique(der_slice)
                for info in enumerate(unique_items): der_slice[der_slice == info[1]] = info[0]
                
                ind_neighbor_slice[der_mask_slice == index] = der_slice
            # reindexing  ind_center and ind_neighbor
            for index, offset in enumerate(natoms_offset):
                ind_center_slice[der_mask_slice == index] += offset
            partial_slice = np.array([sum(inds_mask_slice == ind) for ind in range(len(selected_hashes))])
            partial_offset = np.cumsum(partial_slice) - partial_slice
            for index, offset in enumerate(partial_offset):
                ind_neighbor_slice[der_mask_slice == index] += offset
                
            feed_dict[self.fingerprints_ph[ind_part]] = fingerprints_slice
            feed_dict[self.inds_mask_ph[ind_part]] = inds_mask_slice
            feed_dict[self.ind_neighbor_ph[ind_part]] = ind_neighbor_slice
            feed_dict[self.ind_center_ph[ind_part]] = ind_center_slice
        del afp, species_afp, ind_neighbor, ind_center
        
        der_mat = np.load('tmp_data/der_mat.npy')
        for ind_part, ele_names in enumerate(self.species_partitions):
            der_validkeys = np.in1d(species_derafp, ele_names)  
            der_mat_slice = der_mat[der_validkeys]
            np.save('tmp_data/der_mat_%i' % ind_part, der_mat_slice)
            del der_mat_slice
        del der_mat
        for ind_part, ele_names in enumerate(self.species_partitions):
            der_mat_slice = np.load('tmp_data/der_mat_%i.npy' % ind_part)
            feed_dict[self.der_mat_ph[ind_part]] = der_mat_slice
            
        rmtree('tmp_data')
        
        return feed_dict   

    def nn_model(self, ind, hidden_layers, xavier_initializer):

        fp_dim = self.fp_dim
        activation = eval('tf.nn.' + self.activation)
        if not xavier_initializer: weight, bias = weight_bias_initializer(hidden_layers, fp_dim)
        
        nNeurons = hidden_layers[0]
        if xavier_initializer:
            W = tf.get_variable(
                    "w%i_0" % ind, 
                    [fp_dim, nNeurons], initializer=tf.contrib.layers.xavier_initializer())
            b = tf.get_variable(
                    "b%i_0" % ind, 
                    [nNeurons], initializer=tf.contrib.layers.xavier_initializer())
        else:
            W = tf.get_variable("w%i_0" % ind, initializer=weight[1])
            b = tf.get_variable("b%i_0" % ind, initializer=bias[1])
        self.weights[ind].append(W)
        self.bias[ind].append(b)

        layer_activation = tf.nn.dropout(
                activation(tf.matmul(self.fingerprints_ph[ind], W) + b), rate = self.rate_ph)
        l2_regularization = tf.reduce_sum(tf.square(W))

        if len(hidden_layers) > 1:
            for i in range(1, len(hidden_layers)):
                nNeurons = hidden_layers[i]
                nNeuronsOld = hidden_layers[i - 1]

                if xavier_initializer:
                    W = tf.get_variable(
                            "w%i_%i" % (ind, i),
                            [nNeuronsOld, nNeurons], initializer=tf.contrib.layers.xavier_initializer())
                    b = tf.get_variable(
                            "b%i_%i" % (ind, i),
                            [nNeurons], initializer=tf.contrib.layers.xavier_initializer())
                else:
                    W = tf.get_variable("w%i_%i" % (ind, i), initializer=weight[i + 1])
                    b = tf.get_variable("b%i_%i" % (ind, i), initializer=bias[i + 1])
                self.weights[ind].append(W)
                self.bias[ind].append(b)
                
                layer_activation = tf.nn.dropout(
                    activation(tf.matmul(layer_activation, W) + b), rate = self.rate_ph)
                l2_regularization += tf.reduce_sum(tf.square(W)) + tf.reduce_sum(tf.square(b))

        if xavier_initializer:
            W = tf.get_variable(
                    "w%i_out" % ind,
                    [hidden_layers[-1], 1], initializer=tf.contrib.layers.xavier_initializer())
            b = tf.get_variable(
                    "b%i_out" % ind,
                    [1], initializer=tf.contrib.layers.xavier_initializer())
        else:
            W = tf.get_variable("w%i_out" % ind, initializer=weight[len(hidden_layers) + 1])
            b = tf.get_variable("b%i_out" % ind, initializer=bias[len(hidden_layers) + 1])
        self.weights[ind].append(W)
        self.bias[ind].append(b)
        
        y_out = tf.matmul(layer_activation, W) + b
        l2_regularization += tf.reduce_sum(tf.square(W)) + tf.reduce_sum(tf.square(b))
        nn_out = tf.unsorted_segment_sum(y_out, self.inds_mask_ph[ind], self.batchsize_ph)

        # Gives the contribution
        dEjdgj = tf.gradients(y_out, self.fingerprints_ph[ind])[0]
        dEdg_arranged = tf.gather(dEjdgj, self.ind_neighbor_ph[ind])
        dEdg_arranged_expand = tf.expand_dims(dEdg_arranged, 2)
        dEdg_arranged_tile = tf.tile(dEdg_arranged_expand, [1, 1, 3])
        dEdx = tf.reduce_sum(tf.multiply(dEdg_arranged_tile, self.der_mat_ph[ind]), 1)
        dEdx_arranged = - tf.unsorted_segment_sum(dEdx, self.ind_center_ph[ind], self.totalNumAtoms_ph)

        return nn_out, dEdx_arranged, l2_regularization

    def construct_graph(self, xavier_initializer):

        tf.reset_default_graph()

        fp_dim = self.fp_dim
        self.graph = tf.Graph()

        with tf.device('/device:GPU:0'), self.graph.as_default():
            self.energies_in_ph = tf.placeholder('float32', shape=[None, 1], name='energies_in')
            self.forces_in_ph = tf.placeholder('float32', shape=[None, 3], name='forces_in')
            # the rate that activation is made zero
            self.rate_ph = tf.placeholder('float32', name='rate')
            self.batchsize_ph = tf.placeholder('int32', name='batchsize')

            self.force_coefficient_ph = tf.placeholder('float32', name='force_coefficient')
            self.energy_coefficient_ph = tf.placeholder('float32', name='energy_coefficient')
            self.regularization_strength_ph = tf.placeholder('float32', name='regularization_strength')

            self.natoms_per_image_ph = tf.placeholder('float32', shape=[None, 1], name='natoms_per_image')
            self.totalNumAtoms_ph = tf.placeholder('int32', name='total_natoms')
            
            self.forces_mask = tf.placeholder('int32', name='forces_mask')
            self.forces_atoms_mask = tf.placeholder('int32', name='forces_atoms_mask')
            self.selected_hashes = tf.placeholder('string', name='selected_hashes')
            
            self.fingerprints_ph, self.inds_mask_ph = {}, {}
            self.der_mat_ph, self.ind_neighbor_ph, self.ind_center_ph = {}, {}, {}

            for ind_part, _ in enumerate(self.species_partitions):
                self.fingerprints_ph[ind_part] = tf.placeholder('float32', shape=[None, fp_dim], name='afp_%i' % ind_part)
                self.inds_mask_ph[ind_part] = tf.placeholder('int32', shape=[None], name='inds_mask_%i' % ind_part)

                self.der_mat_ph[ind_part] = tf.placeholder('float32', shape=[None, fp_dim, 3], name='der_mat_%i' % ind_part)
                self.ind_neighbor_ph[ind_part] = tf.placeholder('int32', shape=[None], name='ind_neighbor_%i' % ind_part)
                self.ind_center_ph[ind_part] = tf.placeholder('int32', shape=[None], name='ind_center_%i' % ind_part)

            out_list = [self.nn_model(ind, config, xavier_initializer) for ind, config in enumerate(self.hidden_layers)]
            ene_parts, force_parts, l2_parts = zip(*out_list)

            self.energies = tf.add_n(ene_parts)
            self.forces = 6.90 * tf.add_n(force_parts)
            self.l2_regularization = tf.add_n(l2_parts)

            self.energy_loss = tf.reduce_mean(tf.square(tf.subtract(self.energies_in_ph, self.energies)))
            self.force_loss = 1/ 3 * tf.reduce_mean(tf.square(tf.subtract(self.forces_in_ph, self.forces)))
            self.force_loss_huge = 1 / 3 * tf.reduce_mean(tf.square(3 * tf.subtract(self.forces_in_ph, self.forces)))
            
            self.energy_maxresid = tf.reduce_max(tf.abs(tf.subtract(self.energies_in_ph, self.energies)))
            self.force_maxresid = tf.reduce_max(tf.abs(tf.subtract(self.forces_in_ph, self.forces)))
            
            self.energy_mae = tf.reduce_mean(tf.abs(tf.subtract(self.energies_in_ph, self.energies)))
            self.force_mae = tf.reduce_mean(tf.abs(tf.subtract(self.forces_in_ph, self.forces)))
            
            bad_ind = tf.argmax(tf.abs(tf.subtract(self.forces_in_ph, self.forces)))
            self.bf_x_name = tf.gather(self.selected_hashes, self.forces_mask[bad_ind[0]])
            self.bf_x_atom = self.forces_atoms_mask[bad_ind[0]]
            self.bf_y_name = tf.gather(self.selected_hashes, self.forces_mask[bad_ind[1]])
            self.bf_y_atom = self.forces_atoms_mask[bad_ind[1]]
            self.bf_z_name = tf.gather(self.selected_hashes, self.forces_mask[bad_ind[2]])
            self.bf_z_atom = self.forces_atoms_mask[bad_ind[2]]

            self.loss = tf.add_n(
                [self.force_coefficient_ph * self.force_loss_huge,
                 self.energy_coefficient_ph * self.energy_loss,
                 self.regularization_strength_ph * self.l2_regularization]
                )

    def lbfgs(
            self,
            model_path, 
            db_path, portion, test_size, random_state,
            num_iterations, earlystopping_patience,
            train_params,
            convergence,
            xavier_initializer=True,
            retrain=False):
        
        val_recorder = open('val_loss', 'w')
        train_recorder = open('train_loss', 'w')
        names = ['iter', 'loss', 'ene rmse', 'ene mae', 'ene mr', 'for rmse', 'for mae', 'for mr']
        val_recorder.write('{:^10}{:^15}{:^15}{:^15}{:^15}{:^15}{:^15}{:^15}\n'.format(*names))
        train_recorder.write('{:^10}{:^15}{:^15}{:^15}{:^15}{:^15}{:^15}{:^15}\n'.format(*names))
        train_recorder.flush(); val_recorder.flush()
        train_recorder.close(); val_recorder.close()
                
        def step_callbackfun(x): # x is the raveled model parameter
            
            nonlocal val_holder, iteration
            rate_holder = feeddict_train[self.rate_ph]
            feeddict_train[self.rate_ph] = 0
            vals_train = monitor_func_train(x)
            feeddict_train[self.rate_ph] = rate_holder
            
            with open('train_bad', 'a') as train_bad:
                train_bad.write('%i %s %i %s %i %s %i\n' % (iteration,
                                                   vals_train[-6], vals_train[-5], 
                                                   vals_train[-4], vals_train[-3], 
                                                   vals_train[-2], vals_train[-1], ))
                train_bad.flush()            
            converged_train = check_convergence(
                    iteration=iteration, convergence=convergence, func_vals=vals_train[:-6], file_name='train_loss'
                    )
            if converged_train: raise Convergence_Occurred()
            
            with open('val_bad', 'a') as val_bad:
                vals_val = monitor_func_val(x)
                val_bad.write('%i %s %i %s %i %s %i\n' % (iteration,
                                                          vals_val[-6], vals_val[-5], 
                                                          vals_val[-4], vals_val[-3], 
                                                          vals_val[-2], vals_val[-1], ))
                val_bad.flush()            
            _ = check_convergence(
                    iteration=iteration, convergence=convergence, func_vals=vals_val[:-6], file_name='val_loss'
                    )
            if iteration % earlystopping_patience == 0:
                if vals_val[0] - val_holder >= 1: raise Early_stopping()
                val_holder = vals_val[0]
            iteration += 1
        
        if not isdir(model_path):
            makedirs(model_path)
        self.model_path = model_path
        model_dict = {}
        model_dict['fp_dim'] = self.fp_dim
        model_dict['activation'] = self.activation
        model_dict['species_partitions'] = self.species_partitions
        model_dict['hidden_layers'] = self.hidden_layers
        with open('%s/model_config' % model_path, 'wb') as source:
            pickle.dump(model_dict, source)
        
        self.weights = [[] for id_s, item_s in enumerate(self.species_partitions)]
        self.bias = [[] for id_s, item_s in enumerate(self.species_partitions)]
        self.construct_graph(xavier_initializer=xavier_initializer)
        
        all_hashes = listdir('%s/afp' % db_path)[:portion]
        if random_state is None: random_state = np.random.randint(0, 1e4)
        train_hashes, val_hashes = train_test_split(all_hashes, test_size=test_size, random_state=random_state)
        
        feeddict_train = self.generate_feeddict(db_path, train_hashes, train_params, 'train')
        feeddict_val = self.generate_feeddict(db_path, val_hashes, train_params, 'validation')        
        feeddict_train[self.rate_ph] = train_params['rate']
        feeddict_val[self.rate_ph] = 0
        feeddict_train[self.batchsize_ph] = len(train_hashes)
        feeddict_val[self.batchsize_ph] = len(val_hashes)        
        
        with self.graph.as_default():
            sess = tf.Session()
            saver = tf.train.Saver()
            sess.run(tf.global_variables_initializer())

            if retrain: saver.restore(sess, '%s/model.ckpt' % model_path)
            
            optimizer = ScipyOptimizerInterface(
                self.loss,
                method='l-BFGS-b',
                options={
                        'ftol': 1e-5, 'gtol': 1e-8, 'maxls': 1000000, 'eps': 1e-6, 
                        'maxiter': num_iterations
                        }
                )
            monitored_placeholders = [
                    self.loss, 
                    self.energy_loss, self.energy_mae, self.energy_maxresid,
                    self.force_loss, self.force_mae, self.force_maxresid, 
                    self.bf_x_name, self.bf_x_atom,
                    self.bf_y_name, self.bf_y_atom,
                    self.bf_z_name, self.bf_z_atom
                    ]
            monitor_func_train = optimizer._make_eval_func(
                    tensors=monitored_placeholders,
                    session=sess, 
                    feed_dict=feeddict_train, 
                    fetches=[])
            monitor_func_val = optimizer._make_eval_func(
                    tensors=monitored_placeholders,
                    session=sess, 
                    feed_dict=feeddict_val, 
                    fetches=[])
            
            val_holder = 1e10
            iteration = 0
            try:
                optimizer.minimize(session=sess, 
                                   saved_name='%s/model.ckpt' % model_path, 
                                   feed_dict=feeddict_train, step_callback=step_callbackfun)
            except Convergence_Occurred:
                print('Convergence occurred!'); saver.save(sess, '%s/model.ckpt' % model_path); return 
            except Early_stopping:
                print('Early stopped!'); saver.save(sess, '%s/model.ckpt' % model_path); return 
            finally:
                print('Nothing in particular'); saver.save(sess, '%s/model.ckpt' % model_path)
