import numpy as np

def calc_rotational_constant(bond_lengths, atomic_mass=11.9989, unit='angstrom'):
    """
    Calculate rotational constant B (cm^-1) for a linear carbon chain.
    
    Parameters:
    - bond_lengths: list of bond lengths (in Angstrom or Bohr, specify unit)
    - atomic_mass: atomic mass of carbon in amu (default 12.011)
    - unit: 'angstrom' or 'bohr' (default 'angstrom')
    
    Returns:
    - B_rot: rotational constant in cm^-1
    """
    amu_to_kg = 1.66053906660e-27
    h = 6.62607015e-34
    c = 2.99792458e10
    bohr_to_m = 0.52917721092e-10
    angstrom_to_m = 1e-10

    if unit.lower() == 'angstrom':
        bond_lengths_m = np.array(bond_lengths) * angstrom_to_m
    elif unit.lower() == 'bohr':
        bond_lengths_m = np.array(bond_lengths) * bohr_to_m
    else:
        raise ValueError("Unit must be 'angstrom' or 'bohr'")

    n_atoms = len(bond_lengths) + 1
    masses = np.full(n_atoms, atomic_mass) * amu_to_kg

    positions = np.zeros(n_atoms)
    for i in range(1, n_atoms):
        positions[i] = positions[i-1] + bond_lengths_m[i-1]

    total_mass = masses.sum()
    com = np.sum(masses * positions) / total_mass
    I = np.sum(masses * (positions - com)**2)

    B = h / (8 * np.pi**2 * c * I)

    return B

# Input bond lengths (in Angstrom)
C2_bond_len = [1.2425]  # e.g. [1.2431]
C3_bond_len = [1.2981,1.2981]  # e.g. [1.2966, 1.2966]
C4_bond_len = [1.3135,1.2938,1.3135]  # e.g. [1.2966, 1.2858, 1.2966]
C5_bond_len = [1.2966,1.2858,1.2858,1.2966]  # e.g. [1.2966, 1.2858, 1.2858, 1.2966]
C6_bond_len = [1.3049,1.2905,1.278,1.2905,1.3049]  # e.g. [1.2966, 1.2858, 1.2858, 1.2858, 1.2966]

# Example fill-ins, uncomment to use:
# C2_bond_len = [1.2431]
# C3_bond_len = [1.2966, 1.2966]
# C4_bond_len = [1.2966, 1.2858, 1.2966]
# C5_bond_len = [1.2966, 1.2858, 1.2858, 1.2966]
# C6_bond_len = [1.2966, 1.2858, 1.2858, 1.2858, 1.2966]

for name, bl in zip(['C2', 'C3', 'C4', 'C5', 'C6'],
                    [C2_bond_len, C3_bond_len, C4_bond_len, C5_bond_len, C6_bond_len]):
    if bl:
        B = calc_rotational_constant(bl)
        print(f"Rotational constant B for {name} (cm^-1): {B:.4f}")
    else:
        print(f"No bond lengths provided for {name}, skipping...")

