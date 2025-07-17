import numpy as np
import os
import joblib
from tensorflow.keras.models import load_model
from pes2mp_driver import create_TrainableActivation, create_CustomDecayLayer

# ============ USER TO FILL ============
file_list = [
    ('C2Hm_2.dat',    1.389, 3.450359, -1, 1),        # C2H– (anion)
    ('C2Hm.dat',      1.389, 3.450359, -1, 1),        # C2H– (anion)
    ('C2Nm.dat',      0.3943, 3.621205, -1, 1),       # C2N– (anion)
    ('C3.dat',        0.4305, 3.602107,  0, 1),       # C3 (neutral)  symm
    ('C5Hp.dat',      0.080453792, 3.756167, +1, 1),  # C5H+ (cation, radical)
    ('C6H-.dat',      0.045927, 3.794568, -1, 1),     # C6H– (anion)
    ('CFp.dat',       1.72, 3.544877, +1, 1),         # CF+ (closed-shell cation)
    ('CHp.dat',       13.9534, 3.060779, +1, 1),      # CH+ (radical cation)
    ('HC3N.dat',      0.15174, 3.711386,  0, 1),      # HC3N (neutral)
    ('HCN.dat',       1.478216, 3.486026,  0, 1),     # HCN (neutral)
    ('HNC.dat',       1.512106, 3.486026,  0, 1),     # HNC (neutral)
    ('NCCNH+.dat',    0.148, 3.721617, +1, 1),        # NCCNH+ (closed-shell cation)
    ('AlCN.dat',     0.168, 3.7224, 0, 1),        # NCCNH+ (closed-shell cation)
    ('AlNC.dat',     0.2, 3.7224, 0, 1),        # NCCNH+ (closed-shell cation)
    ('CNm.dat',     1.87239, 3.4687, -1, 1),        # NCCNH+ (closed-shell cation)
    ('CO.dat',     1.925, 3.5014, 0, 1),        # NCCNH+ (closed-shell cation)
    ('HCOp.dat',     1.487501, 3.5171996, 1, 1),        # NCCNH+ (closed-shell cation)
    ('CO2.dat',     0.390219, 3.668845, 0, 1),        # NCCNH+ (closed-shell cation) symm
]

symmetric_mols = {'C5', 'C3', 'CO2'}  # add more if needed


# Create directories
os.makedirs('results', exist_ok=True)
os.makedirs('extra', exist_ok=True)

# Load model & scalers
feature_scaler = joblib.load('feature_scaler.pkl')
target_scaler = joblib.load('target_scaler.pkl')
TrainableActivation = create_TrainableActivation()
CustomDecayLayer = create_CustomDecayLayer()
model = load_model('NN_model_en.keras', custom_objects={
    'TrainableActivation': TrainableActivation,
    'CustomDecayLayer': CustomDecayLayer
})
batch_size = 5120

# ===== Functions =====

#def extract_basecol_data(filepath):
#    """Extract filtered (i,j), T values, and rates from BASECOL format."""
#    with open(filepath, 'r') as f:
#        lines = f.readlines()
#
#    # Locate temperature header line
#    temp_line_idx = None
#    for idx, line in enumerate(lines):
#        try:
#            floats = [float(x) for x in line.strip().split()]
#            if len(floats) > 5:
#                temp_line_idx = idx
#                break
#        except ValueError:
#            continue
#    if temp_line_idx is None:
#        raise ValueError("Temperature header not found.")
#
#    T_values = np.array([float(x) for x in lines[temp_line_idx].strip().split()])
#    valid_temp_indices = T_values <= 101
#    T_values = T_values[valid_temp_indices]
#
#    data = {}
#    for line in lines[temp_line_idx + 1:]:
#        parts = line.strip().split()
#        if len(parts) < 6:
#            continue
#        try:
#            i = int(parts[0])
#            j = int(parts[1])
#        except ValueError:
#            continue
#        if i > j and i < 11 and j < 11:
#            rates = np.array([float(x) for x in parts[4:]])[valid_temp_indices]
#            data[(i, j)] = rates
#    return T_values, data


def extract_basecol_data(filepath):
    """Extract filtered (i,j), T values, and rates from BASECOL format."""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Locate the temperature header line
    temp_line_idx = None
    for idx, line in enumerate(lines):
        try:
            floats = [float(x) for x in line.strip().split()]
            if len(floats) > 2:  # allow even 3 temperatures
                temp_line_idx = idx
                break
        except ValueError:
            continue
    if temp_line_idx is None:
        raise ValueError("Temperature header not found.")

    T_full = np.array([float(x) for x in lines[temp_line_idx].strip().split()])
    valid_temp_indices = T_full <= 101
    T_values = T_full[valid_temp_indices]

    data = {}
    for line in lines[temp_line_idx + 1:]:
        parts = line.strip().split()
        if len(parts) < 6:
            continue
        try:
            i = int(parts[0])
            j = int(parts[1])
        except ValueError:
            continue
        rates_all = np.array([float(x) for x in parts[4:]])
        if rates_all.shape[0] != T_full.shape[0]:
            continue  # skip malformed lines
        if i > j and i < 11 and j < 11:
            rates = rates_all[valid_temp_indices]
            data[(i, j)] = rates

    return T_values, data



#def predict_nn(T_values, transitions, B, mu, charge, multiplicity):
#    """Generate NN predictions for all transitions and T."""
#    all_preds = {}
#    for (i, j) in transitions:
#        # Convert N → J by subtracting 1
#        inputs = [[T, B, mu, i - 1, j - 1, charge, multiplicity] for T in T_values]
#        X = feature_scaler.transform(np.array(inputs))
#        pred_scaled = model.predict(X, batch_size=batch_size, verbose=0)
#        pred = np.exp(target_scaler.inverse_transform(pred_scaled)).flatten()
#        all_preds[(i, j)] = pred  # keys remain as (N, N') = (i, j)
#    return all_preds




def predict_nn(T_values, transitions, B, mu, charge, multiplicity, is_symmetric):
    """Generate NN predictions for all transitions and T."""
    all_preds = {}
    for (i, j) in transitions:
        if is_symmetric:
            inputs = [[T, B, mu, (i - 1)*2, (j - 1)*2, charge, multiplicity] for T in T_values]
            X = feature_scaler.transform(np.array(inputs))
            pred_scaled = model.predict(X, batch_size=batch_size, verbose=0)
            pred = 2 * np.exp(target_scaler.inverse_transform(pred_scaled)).flatten()
        else:
            inputs = [[T, B, mu, i - 1, j - 1, charge, multiplicity] for T in T_values]
            X = feature_scaler.transform(np.array(inputs))
            pred_scaled = model.predict(X, batch_size=batch_size, verbose=0)
            pred = np.exp(target_scaler.inverse_transform(pred_scaled)).flatten()
        all_preds[(i, j)] = pred
    return all_preds


def compute_error(reference, predicted, T_values):
    """Compute % error for each (i,j,T) and return per-line + overall mean."""
    error_lines = []
    all_rel_errors = []
    for (i, j), ref_vals in reference.items():
        if (i, j) not in predicted:
            continue
        pred_vals = predicted[(i, j)]
        with np.errstate(divide='ignore', invalid='ignore'):
            rel_err = np.abs((pred_vals - ref_vals) / ref_vals) * 100
            rel_err = np.nan_to_num(rel_err)
        all_rel_errors.extend(rel_err)
        line = f"{i}\t{j}\t" + "\t".join(f"{e:.2f}" for e in rel_err) + "\n"
        error_lines.append(line)
    MAPE = np.mean(all_rel_errors) if all_rel_errors else 0.0
    return error_lines, MAPE

# ===== Main Loop =====

compiled_error_path = os.path.join('extra', 'compiled_error.dat')
with open(compiled_error_path, 'w') as ce:
    ce.write("File\tB\tmu\tCharge\tMultiplicity\tMeanAbsolutePercentError\n")

for fname, B, mu, charge, multiplicity in file_list:
    base_name = os.path.splitext(os.path.basename(fname))[0]
    print(f"\nProcessing: {fname}")

    # Step 1: Reference data (filtered)
    T_vals, ref_data = extract_basecol_data(fname)
    transitions = list(ref_data.keys())

    cleaned_file = os.path.join('results', f"{base_name}_cleaned.dat")
    with open(cleaned_file, 'w') as f:
        f.write("I\tJ\t" + "\t".join(f"{T:.3f}" for T in T_vals) + "\n")
        for (i, j), rates in ref_data.items():
            f.write(f"{i}\t{j}\t" + "\t".join(f"{r:.2E}" for r in rates) + "\n")
    print(f"Saved cleaned file: {cleaned_file}")

    # Step 2: NN predictions
    is_symmetric = base_name in symmetric_mols
    pred_data = predict_nn(T_vals, transitions, B, mu, charge, multiplicity, is_symmetric)

#    pred_data = predict_nn(T_vals, transitions, B, mu, charge, multiplicity)

    nn_file = os.path.join('results', f"{base_name}_NN.dat")
    with open(nn_file, 'w') as f:
        f.write("I\tJ\t" + "\t".join(f"{T:.3f}" for T in T_vals) + "\n")
        for (i, j), rates in pred_data.items():
            f.write(f"{i}\t{j}\t" + "\t".join(f"{r:.2E}" for r in rates) + "\n")
    print(f"Saved predicted NN rates: {nn_file}")

    # Step 3: Error computation
    error_lines, MAPE = compute_error(ref_data, pred_data, T_vals)

    error_file = os.path.join('extra', f"{base_name}_error.dat")
    with open(error_file, 'w') as f:
        f.write("I\tJ\t" + "\t".join(f"{T:.3f}" for T in T_vals) + "\n")
        f.writelines(error_lines)
    print(f"Saved detailed error file: {error_file}")

    # Step 4: Log overall error
    with open(compiled_error_path, 'a') as ce:
        ce.write(f"{fname}\t{B:.6f}\t{mu:.6f}\t{charge}\t{multiplicity}\t{MAPE:.2f}\n")
