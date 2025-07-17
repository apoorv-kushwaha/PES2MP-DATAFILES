import os
import re
import numpy as np
import pandas as pd

# Folders
carbon_folders = [f"C{n}" for n in range(1, 7)]

# Output collection
output = []

for folder in carbon_folders:
    molscat_path = os.path.join(folder, "MOLSCAT")
    if not os.path.isdir(molscat_path):
        print(f"Skipping {folder}: MOLSCAT directory missing")
        continue

    # === Read redm from 3_molrates_dex.py ===
    molrates_path = os.path.join(molscat_path, "3_molrates_dex.py")
    if not os.path.exists(molrates_path):
        print(f"Skipping {folder}: 3_molrates_dex.py not found")
        continue

    with open(molrates_path) as f:
        content = f.read()
    redm_match = re.search(r"redm\s*=\s*([0-9.]+)", content)
    if not redm_match:
        print(f"Skipping {folder}: redm not found in script")
        continue
    redm = float(redm_match.group(1))

    print(f"\n==== Processing {folder} | redm = {redm} ====")

    # === Run auto_rates.py with updated redm ===
    # Create a modified script from the base template
    with open("auto_rates_template.py") as f:
        lines = f.readlines()

    modified_lines = []
    for line in lines:
        if line.strip().startswith("redm"):
            modified_lines.append(f"redm = {redm}\n")
        else:
            modified_lines.append(line)

    temp_script_path = os.path.join(molscat_path, "run_temp.py")
    with open(temp_script_path, "w") as f:
        f.writelines(modified_lines)

    # Execute the script
    os.chdir(molscat_path)
    os.system("python3 run_temp.py > runlog.txt")
    os.remove("run_temp.py")
    os.chdir("../../")  # Back to root

    # === Read and compare k_dex.dat and k_exc.dat ===
    k_dex_path = os.path.join(molscat_path, "k_dex.dat")
    k_exc_path = os.path.join(molscat_path, "k_exc.dat")
    if not (os.path.exists(k_dex_path) and os.path.exists(k_exc_path)):
        print(f"Skipping {folder}: Rate files missing")
        continue

    # Read both data files
    dex = pd.read_csv(k_dex_path, delim_whitespace=True, comment="#")
    exc = pd.read_csv(k_exc_path, delim_whitespace=True, comment="#")

    # Matching columns and temperatures
    common_cols = [c for c in dex.columns if c in exc.columns and c != 'T']
    for col in common_cols:
        ratio = dex[col] / exc[col].replace(0, np.nan)  # Avoid /0
        for t, r in zip(dex['T'], ratio):
            if np.isnan(r) or np.isinf(r):
                continue
            label = col  # like "0->1"
            output.append([t, folder[1:], r, label])

# === Save final output ===
df = pd.DataFrame(output, columns=["T", "n", "Ratio", "label"])
df.to_csv("rate_ratios_all.csv", index=False)
print("\n✅ Output saved to rate_ratios_all.csv")

