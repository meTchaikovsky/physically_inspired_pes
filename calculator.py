
import pickle

from ase.calculators.calculator import Calculator, all_changes
from .train.model import neural_network
    

class ky_calc(Calculator):
    
    implemented_properties = ['energy', 'forces']
    nolabel = True

    def __init__(self, model_path, des_path, **kwargs):
        
        Calculator.__init__(self, **kwargs)
        self.model_path = model_path
        self.des_path = des_path

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        
        Calculator.calculate(self, atoms, properties, system_changes)
        saved_model = self.model_path
        des_path = self.des_path
        
        with open('%s/model_config' % saved_model, 'rb') as source:
            model_dict = pickle.load(source)
        fp_dim = model_dict['fp_dim']
        activation = model_dict['activation']
        species_partitions = model_dict['species_partitions']
        hidden_layers = model_dict['hidden_layers']
        
        nn_model = neural_network(
                species_partitions = species_partitions,
                fp_dim = fp_dim, hidden_layers = hidden_layers,
                activation = activation,
                )
        energy, forces = nn_model.engine(des_path, saved_model, atoms)      

        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces