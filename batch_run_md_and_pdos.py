import subprocess

# 只处理 seed 43（你可以换成 44）
seed = 43

print(f"\n🚀 Running MD simulation for seed {seed} ...")
subprocess.run([
    "python",
    "run_md_custom.py",
    "--input", f"sqs_seed_{seed}.xyz",
    "--output", f"md_outputs/seed_{seed}.traj"
], check=True)

print(f"\n📊 Calculating PDOS for seed {seed} ...")
subprocess.run([
    "python",
    "velocity_atuocorrelation_function2.py",
    "-traj", f"md_outputs/seed_{seed}.traj",
    "-pdos", "True",
    "-u", "THz"
], check=True)
