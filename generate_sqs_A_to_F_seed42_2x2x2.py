# generate_sqs_A_to_F_seed42_2x2x2_fixed.py
from ase.build import bulk, make_supercell
from ase.io import write
from icet import ClusterSpace
from icet.tools.structure_generation import generate_sqs
import os, time, numpy as np, threading


groups = {
    "A": ["Co","Cr","Fe","Mn","Ni"],
    "B": ["Mg","Co","Ni","Cu","Zn"],
    "C": ["La","Y","Ce","Pr","Nd"],
    "D": ["Zr","Nb","Ti","V","Hf"],
    "E": ["Ba","Sr","Ca","La","Ce"],
    "F": ["Ce","Zr","Hf","Th","U"]
}


def build_template(group):
    if group in ["A","B","C","E","F"]:
        atoms = bulk("MgO","rocksalt",a=4.2)*(2,2,2)   # ✅ rocksalt 2×2×2
    elif group == "D":
        atoms = make_supercell(bulk("V","bcc",a=3.3),np.eye(3)*2)  # ✅ bcc 2×2×2
    return atoms


output_dir = "structures"
os.makedirs(output_dir, exist_ok=True)


for g, els in groups.items():
    atoms = build_template(g)

    
    if g in ["A","B","C","E","F"]:
        chemical_symbols = []
        for atom in atoms:
            if atom.symbol != "O":        
                chemical_symbols.append(els)
            else:
                chemical_symbols.append(["O"])  
    else:  
        chemical_symbols = [els] * len(atoms)

    try:
        cluster_space = ClusterSpace(atoms, chemical_symbols=chemical_symbols, cutoffs=[7,5,4])
    except Exception as e:
        print(f"\n❌ Failed to create ClusterSpace for Group {g}: {e}")
        continue

    target = {el: 1/len(els) for el in els}
    print(f"\n🔹 Generating SQS for Group {g} ({els}) using 2x2x2, seed42")

    seed = 42
    fname = os.path.join(output_dir, f"{g}_seed{seed}.xyz")
    if os.path.exists(fname):
        print(f"   ⏭️  {fname} already exists, skipping.")
        continue

    try:
        print(f"   🔄 Running seed {seed} ... [n_steps=10000]")

        # === 进度条 ===
        running = True
        def show_progress():
            t0 = time.time()
            while running:
                elapsed = int(time.time() - t0)
                print(f"\r   ⏳ Generating {g} seed {seed} ... {elapsed}s elapsed", end="")
                time.sleep(5)
        t = threading.Thread(target=show_progress)
        t.start()

        # === 调用 icet ===
        start = time.time()
        structure = generate_sqs(
            cluster_space=cluster_space,
            target_concentrations=target,
            random_seed=seed,
            n_steps=10000,
            max_size=64
        )

        running = False
        t.join()
        print(f"\r   ✅ Finished {g} seed {seed} in {time.time()-start:.1f}s{' '*20}")

        write(fname, structure)
        print(f"   📦 Saved {fname} ({len(structure)} atoms)")

    except Exception as e:
        running = False
        print(f"\r   ❌ Failed seed {seed} for Group {g}: {e}")
