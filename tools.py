sequence_AA = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"

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

print("L'hydrophobicité est de", calc_hydrophobicity(sequence_AA))

frequency = "att"
string = "attgttgattgccgtattgatt"

def find_and_write_all_frequency_recurrences(string, frequency):
    """
    This function searches for all recurrences of a given frequency of nucleotids in the input string of nucleotides and writes the results to a file.
    """
    start = 0
    result = []
    while True:
        index = string.find(frequency, start)
        if index == -1:
            break
        result.append(index)
        start = index + 1
    
    with open("frequency_recurrence_result.txt", "w") as file:
        if result:
            file.write("The frequency '{}' was found at the following positions in the input string:\n".format(frequency))
            for i, pos in enumerate(result):
                file.write("Position {}: {}\n".format(i+1, pos))
        else:
            file.write("The frequency '{}' was not found in the input string.".format(frequency))

#find_and_write_all_frequency_recurrences(string, frequency)

def calc_pI(sequence):
    # Define pK ranges for individual amino acids
    pK_ranges = {
        'K': (10.4, 10.8),
        'R': (12.0, 12.5),
        'H': (6.0, 7.4),
        'D': (3.8, 4.4),
        'E': (4.0, 4.8),
        'C': (8.0, 8.5),
        'Y': (9.6, 10.1),
        'N': (7.5, 8.7),
        'Q': (5.0, 6.0),
        'W': (5.2, 5.4),
        'S': (5.6, 5.9),
        'T': (5.6, 6.1),
        'M': (5.5, 5.8),
        'A': (6.0, 7.6),
        'F': (4.0, 4.5),
        'G': (9.0, 10.0),
        'I': (6.0, 6.5),
        'L': (6.0, 6.5),
        'P': (7.2, 8.8),
        'V': (6.0, 6.5),
    }

    # Define function to calculate net charge at given pH
    def calc_net_charge(pH, pK_ranges, sequence):
        net_charge = 0
        for aa in sequence:
            if aa in pK_ranges:
                pK_min, pK_max = pK_ranges[aa]
                if pH <= pK_min:
                    net_charge += 1
                elif pH >= pK_max:
                    net_charge -= 1
                else:
                    net_charge += (10**(pK_min - pH) - 10**(pK_max - pH)) / (1 + 10**(pK_min - pH) + 10**(pH - pK_max))
        return net_charge

    # Use bisection method to find pI with 0.01 pH unit precision
    a = 0
    b = 14
    while (b - a) > 0.01:
        c = (a + b) / 2
        if calc_net_charge(c, pK_ranges, sequence) > 0:
            a = c
        else:
            b = c

    return (a + b) / 2

print("Le point isoélectrique est d'environ", calc_pI(sequence_AA), ", ce résultat peut varier en fonction des conditions environnementales")
