import argparse
import math
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory
import pandas as pd
import ase.io.vasp
from rdfpy import rdf
from scipy.ndimage import gaussian_filter
from tqdm import tqdm
from scipy.fftpack import fft, fftfreq
from sklearn.preprocessing import StandardScaler

k_B_ev = 8.617e-5
hbar_ev = 4.136e-15
N_A = 6.0221409e+23

def parse_arguments():
    parser = argparse.ArgumentParser(description="Parser for trajectory data.")
    parser.add_argument('-traj', '--trajectory', type=str, required=True, help="Path to the trajectory file")
    parser.add_argument('-t', '--temp', type=float, help="Temperature in Kelvin (default: 300)")
    parser.add_argument('-pdos', '--phonon_density_of_state', type=bool, help="Generate phonon density of state")
    parser.add_argument('-u', '--units', type=str, help="units options: THz, mev, cm-1")
    return parser.parse_args()

class Traj_Evaluation:
    def __init__(self, Trajectory_path):
        self.df = None
        self.Trajectory_path = Trajectory_path
        self.atoms_data = self._extract_atoms_data()
        self._Loading()

    def _extract_atoms_data(self):
        atoms_data = []
        traj = Trajectory(self.Trajectory_path)
        traj = traj[1:]  # skip first frame

        for atoms in traj:
            atoms_data.append({
                'atom': atoms,
                'temperature': int(atoms.get_temperature()),
                'forces': atoms.get_forces(),
                'energy': atoms.get_potential_energy(),
                'stress': atoms.get_stress(),
                'n_atoms': len(atoms),
                'positions': atoms.get_positions(),
                'scaled_positions': atoms.get_scaled_positions(),
                'velocity': atoms.get_velocities()
            })
        return atoms_data

    def _Loading(self):
        results = []
        for atom_data in self.atoms_data:  # 顺序处理，避免 multiprocessing 占用大量内存
            results.append({
                'Temperature': atom_data['temperature'],
                'Force_per_atom': np.linalg.norm(atom_data['forces']) / atom_data['n_atoms'],
                'Energy_per_atom': atom_data['energy'] / atom_data['n_atoms'],
                'Stress_per_atom': np.linalg.norm(atom_data['stress']) / atom_data['n_atoms'],
                'n_atoms': atom_data['n_atoms'],
                'positions': atom_data['positions'],
                'scaled_positions': atom_data['scaled_positions'],
                'velocity': atom_data['velocity'],
            })
        self.df = pd.DataFrame(results)

    def velocity_autocorrelation_function(self, potim=8):
        num_atoms = self.df['n_atoms'][0]
        num_frames = len(self.df)
        pos = np.array(
            [[self.df['scaled_positions'][i][j] for j in range(num_atoms)] for i in range(num_frames)]
        ).reshape((-1, num_atoms, 3))

        dpos = np.diff(pos, axis=0)
        _lv = self.atoms_data[0]['atom'].cell
        lv = np.array(_lv)
        dpos[dpos > 0.5] -= 1.0
        dpos[dpos < -0.5] += 1.0

        for i in range(len(dpos)):
            dpos[i, :, :] = np.dot(dpos[i], lv) / potim

        velocity = dpos
        VAF = np.zeros((len(velocity) * 2 - 1))
        for i in range(num_atoms):
            for j in range(3):
                VAF += np.correlate(velocity[:, i, j], velocity[:, i, j], 'full')
        VAF /= np.sum(velocity ** 2)
        return VAF

    def get_phonon_dos(self, potim, nblock, unit):
        potim = potim * nblock
        N = len(self.df) - 1
        omega = fftfreq(2 * N - 1, potim * 1e-15) * 1E-12
        if unit.lower() == 'cm-1':
            omega *= 33.356
        elif unit.lower() == 'mev':
            omega *= 4.13567
        VAF2 = self.velocity_autocorrelation_function(potim=potim)
        pdos = np.abs(fft(VAF2 - np.mean(VAF2))) ** 2
        return omega[:N], pdos[:N]

if __name__ == '__main__':
    args = parse_arguments()
    print(f"📂 Loading trajectory: {args.trajectory}")
    traj = Traj_Evaluation(args.trajectory)
    if args.phonon_density_of_state:
        omega, ph_dos = traj.get_phonon_dos(potim=1, nblock=1, unit=args.units)
        ph_dos = gaussian_filter(ph_dos, sigma=2)

        pd.DataFrame([omega, ph_dos]).to_csv('ph_dos.csv')
        scaler = StandardScaler()
        normalized_ph_dos = scaler.fit_transform(ph_dos.reshape(-1, 1))

        plt.plot(omega, normalized_ph_dos, 'b-', alpha=1, label='phonon DOS (MLFF)')
        plt.xlabel('Frequency (THz)')
        plt.ylabel('phonon DOS')
        plt.xlim(0, 50)
        plt.legend()
        plt.tight_layout()
        plt.savefig('phonon_dos.pdf', format="pdf")
        plt.show()
