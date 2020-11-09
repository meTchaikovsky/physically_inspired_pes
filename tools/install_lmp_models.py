
import argparse

from kyzhang.lammps_md.model_extractor import model_extractor

parser = argparse.ArgumentParser()
parser.add_argument('--des_path', help="The path to the database.", type=str)
parser.add_argument('--model_paths', help="The path to the model.", nargs='+', type=str)

args = parser.parse_args()
des_path = args.des_path
model_paths = args.model_paths

print(des_path)
print(model_paths)

model_extractor(model_paths=model_paths, des_path=des_path)

