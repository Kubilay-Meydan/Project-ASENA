def calc_pI(sequence):
    # Define pK values for individual amino acids
    pK_values = {
        'K': 10.8,
        'R': 12.5,
        'H': 6.5,
        'D': 3.9,
        'E': 4.1,
        'C': 8.5,
        'Y': 10.1,
        'N': 8.1,
        'Q': 5.5,
        'W': 5.4,
        'S': 5.7,
        'T': 6.3,
        'M': 5.7,
    }

    pI = 0
    for pH in range(0, 14):
        charge = 0
        for aa in sequence:
            if aa in pK_values:
                pK = pK_values[aa]
                if pH < pK:
                    charge += 10**(pK - pH)
                elif pH > pK:
                    charge += 10**(pH - pK) - 1
        if charge == 0:
            pI = pH
            break

    return pI



def calc_hydrophobicity(sequence):
    # Define hydrophobicity values for individual amino acids
    hydrophobicity_values = {
        'A': 1.8,
        'R': -4.5,
        'N': -3.5,
        'D': -3.5,
        'C': 2.5,
        'Q': -3.5,
        'E': -3.5,
        'G': -0.4,
        'H': -3.2,
        'I': 4.5,
        'L': 3.8,
        'K': -3.9,
        'M': 1.9,
        'F': 2.8,
        'P': -1.6,
        'S': -0.8,
        'T': -0.7,
        'W': -0.9,
        'Y': -1.3,
        'V': 4.2
    }

    hydrophobicity = 0
    for aa in sequence:
        if aa in hydrophobicity_values:
            hydrophobicity += hydrophobicity_values[aa]

    return hydrophobicity