from ase.io import read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from mace.calculators import mace_mp

# 1. Load structure
atoms = read("outputs/structure.xyz")

# 2. Load pretrained model
calc = mace_mp(model="medium-mpa-0", device="cpu")
atoms.calc = calc

# 3. Set temperature
MaxwellBoltzmannDistribution(atoms, temperature_K=300)

# 4. MD simulation
dyn = VelocityVerlet(atoms, timestep=1.0 * units.fs)

# 5. Define callback function for printing
def print_status(a=atoms):
    epot = a.get_potential_energy()
    print(f"MD step: {dyn.nsteps:3d} | Energy: {epot:.6f} eV")

# 6. Run 100 steps, call print every step
dyn.attach(print_status, interval=1)
dyn.run(100)

# 7. Save final structure
write("outputs/final_md_structure.xyz", atoms)
