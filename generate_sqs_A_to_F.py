# generate_sqs_A_to_F_compatible.py
from ase.build import bulk, make_supercell
from ase.io import write
from icet import ClusterSpace
from icet.tools.structure_generation import generate_sqs
import os
import numpy as np
import time

# ==== 定义组 A–F ====
group_definitions = {
    "A": ["Co", "Cr", "Fe", "Mn", "Ni"],
    "B": ["Mg", "Co", "Ni", "Cu", "Zn"],
    "C": ["La", "Y", "Ce", "Pr", "Nd"],
    "D": ["Zr", "Nb", "Ti", "V", "Hf"],
    "E": ["Ba", "Sr", "Ca", "La", "Ce"],
    "F": ["Ce", "Zr", "Hf", "Th", "U"]
}

# ==== 定义晶体模板（用已有能运行的结构近似） ====
def build_template(group):
    if group in ["A", "B", "C"]:  # 使用 rocksalt 2x2x2
        atoms = bulk("MgO", crystalstructure="rocksalt", a=4.2)
        atoms *= (2, 2, 2)
    elif group == "D":  # BCC 2x2x2
        atoms = bulk("V", crystalstructure="bcc", a=3.3)
        atoms = make_supercell(atoms, np.eye(3) * 2)
    elif group == "E":  # perovskite 近似为 rocksalt
        atoms = bulk("MgO", crystalstructure="rocksalt", a=4.2)
        atoms *= (2, 2, 2)
    elif group == "F":  # fluorite 近似为 rocksalt
        atoms = bulk("MgO", crystalstructure="rocksalt", a=4.2)
        atoms *= (2, 2, 2)
    return atoms

# ==== 输出目录 ====
output_dir = os.path.join("..", "structures")
os.makedirs(output_dir, exist_ok=True)

# ==== 生成所有组和 seed ====
for group, elements in group_definitions.items():
    atoms = build_template(group)
    cutoffs = [7.0, 5.0, 4.0]
    chemical_symbols = [elements for _ in range(len(atoms))]
    cluster_space = ClusterSpace(atoms, chemical_symbols=chemical_symbols, cutoffs=cutoffs)
    target_concentrations = {el: 1 / len(elements) for el in elements}

    print(f"\n🔹 Generating SQS for Group {group} ({elements})")

    for seed in [42, 43, 44]:
        try:
            print(f"   🔄 Seed {seed} ...")
            start = time.time()
            structure = generate_sqs(
                cluster_space=cluster_space,
                target_concentrations=target_concentrations,
                random_seed=seed,
                n_steps=10000,
                max_size=64  # ✅ 保留旧版参数
            )
            fname = os.path.join(output_dir, f"{group}_seed{seed}.xyz")
            write(fname, structure)
            print(f"   ✅ Saved {fname} | {len(structure)} atoms | Time {time.time()-start:.1f}s")
        except Exception as e:
            print(f"   ❌ Failed seed {seed} for Group {group}: {e}")

print("\n🎯 All A–F groups generated successfully!")
