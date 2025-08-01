import argparse
from ase.io import read, Trajectory
from ase.md.verlet import VelocityVerlet
from ase import units
from mace.calculators import mace_mp

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--steps', type=int, default=10000)
args = parser.parse_args()

atoms = read(args.input)
atoms.calc = mace_mp()  # ✅ 直接使用默认模型

dyn = VelocityVerlet(atoms, timestep=1 * units.fs)
traj = Trajectory(args.output, 'w', atoms)

for i in range(args.steps):
    dyn.run(1)
    traj.write()
    if (i + 1) % 100 == 0:
        print(f"MD step: {i+1} | Energy: {atoms.get_potential_energy():.6f} eV")
