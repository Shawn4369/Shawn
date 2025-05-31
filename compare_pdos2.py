import pandas as pd
import matplotlib.pyplot as plt

# 定义文件路径和标签
files = {
    'seed_42': 'ph_dos_42.csv',
    'seed_43': 'ph_dos_43.csv',
    'seed_44': 'ph_dos_44.csv'
}

# 设置颜色
colors = {
    'seed_42': 'red',
    'seed_43': 'green',
    'seed_44': 'blue'
}

plt.figure(figsize=(10, 7))

for label, file in files.items():
    try:
        df = pd.read_csv(file, header=None)
        freq = df.iloc[1, :].values.astype(float)
        dos = df.iloc[2, :].values.astype(float)

        # 只保留频率在 0 到 10 THz 的部分
        mask = (freq >= 0) & (freq <= 10)
        plt.plot(freq[mask], dos[mask], label=label, color=colors[label])
        print(f"✅ Loaded and plotted (0–10 THz): {file}")
    except Exception as e:
        print(f"❌ Error reading {file}: {e}")

plt.title('Comparison of Phonon DOS for Different Seeds (0–10 THz)')
plt.xlabel('Frequency (THz)')
plt.ylabel('Phonon DOS')
plt.xlim(0, 10)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('compare_pdos_output.png')
print("📊 Plot saved as: compare_pdos_output.png")
