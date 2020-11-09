
import argparse

import numpy as np
from ase.io import read as r 
from kyzhang.calculator import ky_calc

parser = argparse.ArgumentParser()
parser.add_argument('--des_path', help='The path to the database', type=str)
parser.add_argument('--model_paths', help='The path to the models', nargs='+', type=str)

args = parser.parse_args()
des_path = args.des_path
model_paths = args.model_paths

tmp_atoms = r('one.cif')

energies, forces = [], []
for model_path in model_paths:
    trained_calc = ky_calc(model_path=model_path, des_path=des_path)
    tmp_atoms.set_calculator(trained_calc)
    
    energies.append(tmp_atoms.get_potential_energy())
    forces.append(tmp_atoms.get_forces())
    
    del trained_calc
    
np.save('forces', np.mean(forces, axis=0))
print(np.mean(energies))
