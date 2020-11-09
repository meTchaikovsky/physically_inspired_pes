Version physically_inspired_pes 1.0, Nov 9, 2020.

physically_inspired_pes: A Python package for constructing Behler Parrinello neural network potential energy surfaces (BPNN PES) with physically inspired atom-centered symmetry functions (ACSFs). 

References:  

Kangyu Zhang; Lichang Yin; Gang Liu, Computational Materials Science, DOI: 10.1016/j.commatsci.2020.110071


Installation
-------------
This program is written in Python and Fortran 90. The computationally intensive part, i.e., compute the atomic fingerprints and the derivates of the atomic fingerprints w.r.t the coordinates of atoms, is done by Fortran 90. The Fortran source file for computing the atomic fingerprints and the derivates is `descriptor/char_zeta.F90`, it must be compiled by 

    f2py -c char_zeta.F90 -m f_calc
    
to get a python module before calling the module `descriptor`.

A trained PES model can be used to run molecular dynamics (MD) simulations via LAMMPS. We use `kim-api` to interface with `LAMMPS`, therefore, `kim-api-2.0.2 package` must be installed. In order to generate a PES that can be used by `LAMMPS`, run the command below in `lammps_md/model_driver` (see "Install a trained PES model for LAMMPS" for more detail)

    kim-api- collections- management install user driver && kim-api- collections- management install user models

The modules:
- descriptor              For computing atomic fingerprints and the derivates with our physically inspired ACSFs.
- train                   For training PES models with the atomic fingerprints and the derivates. 
- prepare_structures      For creating an `ase` trajectory object from `XDATCAR` and `OUTCAR` of a VASP MD simulation. 
- lammps_md               For installing a trained PES model for LAMMPS.
- tools                   A dir of runnable python script.
- calculator.py           A module for creating an `ase` Calculator object for a trained PES model.
- utilities.py            Miscellaneous tools. 


Train a PES with the physically inspired ACSFs
-------------

See `example_scripts/create_db` and `example_scripts/train_models`.

 Install a trained PES model for LAMMPS
-------------

1. Run `python tools/install_lmp_models.py --des_path $path_to_des --model_paths $path_to_m0 $path_to_m1 ...` to create a text file `model.param` of the physically inspired ACSFs and the trained PES model. 

2. Move `model.param` file to `lammps_md/model_driver/models`

3. In `lammps_md/model_driver`, run the command `kim-api- collections- management install user driver && kim-api- collections- management install user models`. Make sure the `DRIVER_NAME` of `models` is the same as that of `drivers`.


About
-------------
The physically_inspired_pes code was created at Institute of Metal Research (IMR) by Kangyu Zhang. Please feel free to contact Kangyu Zhang (zhangky17s@imr.ac.cn) or open issues in the Github in case of any problems/comments/suggestions in using the code.
