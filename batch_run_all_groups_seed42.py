import os
import subprocess

# Define directories
base_dir = os.getcwd()
structures_dir = os.path.join(base_dir, "structures")
md_dir = os.path.join(base_dir, "md_outputs")
pdos_dir = os.path.join(base_dir, "pdos_outputs")

os.makedirs(md_dir, exist_ok=True)
os.makedirs(pdos_dir, exist_ok=True)

# Collect all seed42 xyz files
xyz_files = sorted([f for f in os.listdir(structures_dir) if f.endswith("seed42.xyz")])
if not xyz_files:
    print("No seed42 structure files found in 'structures' directory.")
    exit()

for xyz in xyz_files:
    name = xyz.replace(".xyz", "")
    xyz_path = os.path.join(structures_dir, xyz)
    traj_path = os.path.join(md_dir, f"{name}.traj")
    pdos_csv = os.path.join(pdos_dir, f"ph_dos_{name}.csv")
    pdos_pdf = os.path.join(pdos_dir, f"ph_dos_{name}.pdf")

    print(f"\n>>> Running MD for {xyz}")
    subprocess.run([
        "python", "run_md_custom.py",
        "--input", xyz_path,
        "--output", traj_path,
        "--steps", "10000"
    ], check=True)

    print(f">>> Calculating PDOS for {name}")
    subprocess.run([
        "python", "velocity_atuocorrelation_function3.py",
        "-traj", traj_path,
        "-pdos", "True",
        "-u", "THz"
    ], check=True)

    # Rename PDOS outputs if they exist
    if os.path.exists("ph_dos.csv"):
        os.rename("ph_dos.csv", pdos_csv)
    if os.path.exists("phonon_dos.pdf"):
        os.rename("phonon_dos.pdf", pdos_pdf)

    print(f"*** Completed processing {name} ***")

print("\nAll seed42 MD and PDOS calculations finished successfully.")
