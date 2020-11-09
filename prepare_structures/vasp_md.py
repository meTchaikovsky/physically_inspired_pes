
from ase.io import read
from ase.io.trajectory import TrajectoryWriter


def vasp_md(num_atoms, mode='a', traj_name='sampled.traj'):
    """
    :param num_atoms: The number of atoms in the simulated system
    :param mode: The mode to write Trajectory file, can be w or a.
    :param traj_name: The name of the traj file to write.
    :return: Nothing to return, but a file will be created.
    """
    with open('OUTCAR', 'r') as ene_outcar:
        all_energies = []
        for line in ene_outcar:
            if 'ion-electron   TOTEN  =' in line:
                free_energy = line.split()[4]
                all_energies.append(float(free_energy))
    image_force = []
    with open('OUTCAR', 'r') as force_outcar:
        all_forces = []
        marker, counter = 0, 0
        for lines in force_outcar:
            if marker == 1:
                image_force.append([float(item) for item in lines.split()[3:]])
                counter += 1
            if counter > num_atoms - 1:
                all_forces.append(image_force)
                marker = 0
                counter = 0

            if marker == 0:
                if 'TOTAL-FORCE (eV/Angst)' in lines:
                    marker = 1
                    image_force = []
                    force_outcar.readline()

    md_traj = read('XDATCAR', ':', format='vasp-xdatcar')
    with TrajectoryWriter(traj_name, mode) as source:
        for index, config in enumerate(md_traj):
            atoms_energy = all_energies[index]
            atoms_force = all_forces[index]
            source.write(config, energy=atoms_energy, forces=atoms_force)