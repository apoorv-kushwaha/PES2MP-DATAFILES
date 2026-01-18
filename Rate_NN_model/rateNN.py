import numpy as np
import os
import joblib
from tensorflow.keras.models import load_model
from pes2mp_driver import create_TrainableActivation, create_CustomDecayLayer

# ============ USER INPUT ============
file_list = [
    ('C2Nm.dat',   0.3943, 3.621205, -1),
    ('C3.dat',     0.4305, 3.602107,  0),
    ('C5Hp.dat',   0.080453792, 3.756167, +1),
    ('C6H-.dat',   0.045927, 3.794568, -1),
    ('HC3N.dat',   0.15174, 3.711386,  0),
    ('NCCNH+.dat', 0.148, 3.721617, +1),
    ('AlCN.dat',   0.168, 3.7224, 0),
    ('AlNC.dat',   0.2,   3.7224, 0),
    ('HCOp.dat',   1.487501, 3.5171996, 1),
    ('CO2.dat',    0.390219, 3.668845, 0),
    ]

symmetric_mols = {'C3', 'C5', 'CO2', 'pH2_C3', 'oH2_C3'}

current_dir = os.getcwd()
os.makedirs('results', exist_ok=True)
os.makedirs('predicted_data', exist_ok=True)

# ===== Load model =====
feature_scaler = joblib.load('feature_scaler.pkl')
target_scaler = joblib.load('target_scaler.pkl')

TrainableActivation = create_TrainableActivation()
CustomDecayLayer = create_CustomDecayLayer()

model = load_model(
    'NN_model_en.keras',
    custom_objects={
        'TrainableActivation': TrainableActivation,
        'CustomDecayLayer': CustomDecayLayer
    }
)

batch_size = 5120

#  Extract valid ranhge from BASECOL (Cleaned files)

def is_valid_transition(i, j, symmetric):
    """
    Transition selection rules:
    - Non-symmetric: i > j, i < 11, j < 11
    - Symmetric:     i > j, i < 6,  j < 6
    """
    limit = 6 if symmetric else 11
    return i > j and i < limit and j < limit


def extract_basecol_data(filepath, symmetric):
    with open(filepath) as f:
        lines = f.readlines()

    # Locate temperature header
    for idx, line in enumerate(lines):
        try:
            floats = [float(x) for x in line.split()]
            if len(floats) > 2:
                temp_idx = idx
                break
        except ValueError:
            continue
    else:
        raise RuntimeError("Temperature header not found")

    T_full = np.array(lines[temp_idx].split(), dtype=float)
    valid = T_full <= 100
    T_vals = T_full[valid]

    data = {}
    for line in lines[temp_idx + 1:]:
        p = line.split()
        if len(p) < 6:
            continue

        try:
            i, j = int(p[0]), int(p[1])
        except ValueError:
            continue

        rates = np.array(p[4:], float)
        if rates.size != T_full.size:
            continue

        if is_valid_transition(i, j, symmetric):
            data[(i, j)] = rates[valid]

    return T_vals, data

#  Predict NN rates
def predict_nn(T_vals, transitions, B, mu, charge, symmetric):
    preds = {}
    for i, j in transitions:
        if symmetric:
            X = [[T, B, mu, (i - 1) * 2, (j - 1) * 2, charge] for T in T_vals]
        else:
            X = [[T, B, mu, i - 1, j - 1, charge] for T in T_vals]

        X = feature_scaler.transform(np.asarray(X))
        y_scaled = model.predict(X, batch_size=batch_size, verbose=0)
        preds[(i, j)] = np.exp(
            target_scaler.inverse_transform(y_scaled)
        ).ravel()

    return preds

#  calculate MAEF
def compute_factor(ref, pred):
    lines, all_f = [], []

    for (i, j), ref_r in ref.items():
        pr = pred[(i, j)]
        mask = ref_r >= 1e-100

        if mask.any():
            f = pr[mask] / ref_r[mask]
            f = np.maximum(f, 1 / f)
        else:
            f = np.ones_like(pr)

        all_f.extend(f)
        lines.append(
            f"{i}\t{j}\t" + "\t".join(f"{x:.3f}" for x in f) + "\n"
        )

    return lines, np.mean(all_f) if all_f else 1.0


# = driver file ==

data_folder = 'basecol_data'  # your folder with .dat files


compiled = os.path.join(current_dir, 'compiled_maef.dat')
with open(compiled, 'w') as f:
    f.write("File\tB\tmu\tCharge\tMeanAbsoluteFactorError\n")

for fname, B, mu, charge in file_list:
    base = os.path.splitext(fname)[0]
    symmetric = base in symmetric_mols

    print(f"Processing {fname}")

    T_vals, ref = extract_basecol_data(os.path.join('basecol_data', fname), symmetric)
    transitions = list(ref.keys())

    # Cleaned reference
    with open(f"predicted_data/{base}_cleaned.dat", "w") as f:
        f.write("I\tJ\t" + "\t".join(f"{T:.3f}" for T in T_vals) + "\n")
        for (i, j), r in ref.items():
            f.write(f"{i}\t{j}\t" + "\t".join(f"{x:.2E}" for x in r) + "\n")

    # NN prediction
    preds = predict_nn(T_vals, transitions, B, mu, charge, symmetric)

    with open(f"predicted_data/{base}_NN.dat", "w") as f:
        f.write("I\tJ\t" + "\t".join(f"{T:.3f}" for T in T_vals) + "\n")
        for (i, j), r in preds.items():
            f.write(f"{i}\t{j}\t" + "\t".join(f"{x:.2E}" for x in r) + "\n")

    # Error
    err_lines, mape = compute_factor(ref, preds)

    with open(f"results/{base}_error.dat", "w") as f:
        f.write("I\tJ\t" + "\t".join(f"{T:.3f}" for T in T_vals) + "\n")
        f.writelines(err_lines)

    with open(compiled, "a") as f:
        f.write(f"{fname}\t{B:.6f}\t{mu:.6f}\t{charge}\t{mape:.2f}\n")

print("\n All files processed successfully.")
