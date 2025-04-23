from ase.io import read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from mace.calculators import mace_mp  

atoms = read("outputs/structure.xyz")

calc = mace_mp(model="medium-mpa-0", device="cpu")
atoms.set_calculator(calc)

MaxwellBoltzmannDistribution(atoms, temperature_K=300)

dyn = VelocityVerlet(atoms, timestep=1.0 * units.fs)

for step in range(100):
    dyn.run(1)
    energy = atoms.get_potential_energy()
    print(f"Step {step + 1}: Energy = {energy:.6f} eV")

write("outputs/final_md_structure.xyz", atoms)
