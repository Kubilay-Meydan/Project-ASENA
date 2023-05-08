from Bio import Entrez, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO, Phylo
from Bio.Align.Applications import MuscleCommandline
import subprocess
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from ete3 import Tree


def search_a_pattern(seq, pattern):
    positions = []
    count = 0
    for i in range(len(seq)):
        if seq[i:i+len(pattern)] == pattern:
            positions.append(i+1)
            count += 1
    return pattern, seq, positions, count

def search_sequence(sequence, pattern):
        # Search for the pattern in the sequence
        positions = []
        count = 0
        for i in range(len(sequence)):
            if sequence[i:i+len(pattern)] == pattern:
                positions.append(i+1)
                count += 1

        # Return the number of occurrences and positions of the pattern in the sequence
        return count, positions

def is_valid_sequence(seq):
    amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    return all(aa in amino_acids for aa in seq)

def is_valid_enter_DNA(seq):
    nucleotids = set("ATGCatgc")
    return all(nc in nucleotids for nc in seq)

def is_valid_enter_RNA(seq):
    nucleotids = set("AUGCaugc")
    return all(nc in nucleotids for nc in seq)

def DNA_to_RNA(entry):
    ans = ''
    for char in entry.lower():
        if char == 'a':
            ans+=('a')
        if char == 'c':
            ans+=('c')
        if char == 'g':
            ans+=('g') 
        if char == 't':
            ans+=('u')
        if char == 'A':
            ans+=('A')
        if char == 'C':
            ans+=('C')
        if char == 'G':
            ans+=('G') 
        if char == 'T':
            ans+=('U')
    return ans

def get_genbank_info(gene_id):
    with open('results/'+str(gene_id), 'w') as f:
        handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        cds_sequence = []
        numberexon = 1
        numbercds = 1
        numberintron = 1
        f.write(record.name)
        f.write('\n')
        f.write('\n')
        f.write(record.description)
        f.write('\n')
        f.write('\n')
        for feature in record.features:
            if feature.type == "exon":
                f.write('Exon ')
                f.write(str(numberexon))
                f.write('\n')
                f.write(str(feature.location))                
                f.write('\n')
                f.write(str(feature.extract(record).seq))
                f.write('\n')
                f.write('\n')
                numberexon += 1
            if feature.type == "CDS":
                f.write('CDS ')
                f.write(str(numbercds))
                f.write('\n')
                f.write(str(feature.location))
                f.write('\n')
                f.write(str(feature.extract(record).seq))
                f.write('\n')
                f.write('\n')
                numbercds += 1
            if feature.type == "intron":
                f.write('Intron ')
                f.write(str(numberintron))
                f.write('\n')
                f.write(str(feature.location))
                f.write('\n')
                f.write(str(feature.extract(record).seq))
                f.write('\n')
                f.write('\n')
                numberintron += 1

def get_sequence(acc, seq_type='protein'):
    Entrez.email = 'your_email@example.com' # required by NCBI
    handle = Entrez.efetch(db='protein', id=acc, rettype='gb', retmode='text')
    record = SeqIO.read(handle, 'genbank')
    if seq_type == 'protein':
        return record.seq
    elif seq_type == 'DNA':
        return record.seq.transcribe()
    else:
        return None

def all_sequences(accessions):
    sequences = []
    for e in accessions:
        sequences.append(str(get_sequence(e)))
    return sequences

def write_fasta(sequences, names):
    with open("results/output.fasta", "w") as fasta:
        for seq, name in zip(sequences, names):
            fasta.write(">" + name + "\n")
            fasta.write(seq + "\n")

def Align_muscle(sequences, output):
    """
    Aligns sequences using MUSCLE and outputs an aligned fasta file.

    Parameters:
    sequences (str): The path to the file containing sequences in fasta format.
    output (str): The path to the file to write the aligned sequences in fasta format.

    Returns:
    str: The path to the aligned fasta file.
    """

    # Align the sequences using MUSCLE
    muscle_exe = "muscle.exe"
    cline = MuscleCommandline(muscle_exe, input=sequences, out=output, clw=True)
    cline()

    # Convert the output to aligned fasta format
    alignment = AlignIO.read(output, "clustal")
    AlignIO.write(alignment, output + ".fasta", "fasta")

    return output

def run_bmge_on_alignment(muscle_file_path, output_file_path):
    # Run BMGE on the input alignment file
    bmge_command = ['java', '-jar', 'BMGE.jar', '-i', muscle_file_path, '-t', 'aa', '-o', output_file_path]
    bmge_output = subprocess.check_output(bmge_command, universal_newlines=True)

    print('BMGE filtering complete. Filtered alignment saved to: ' + 'results/'+output_file_path)

def make_phylo_tree_newick(fasta_file):
    # Load the FASTA file and perform multiple sequence alignment
    alignment = AlignIO.read(fasta_file, "fasta")

    # Calculate pairwise distance matrix
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)

    # Construct a neighbor-joining tree from the distance matrix
    constructor = DistanceTreeConstructor()
    
    nj_tree = constructor.nj(dm)

    # Convert the tree to a Newick format string
    newick_string = nj_tree.format("newick")

    # Save the Newick format tree to a file
    with open("results/my_tree_file.txt", "w") as f:
        f.write(newick_string)

    # Load the Newick format tree file into ETE toolkit
    tree = Tree("results/my_tree_file.txt", format = 1)
    phylo_tree = Phylo.read("results/my_tree_file.txt", "newick")
    for clade in phylo_tree.find_clades():
        if clade.name.startswith("Inner"):
            clade.name = ""
    tree.ladderize()
    Phylo.draw(phylo_tree, branch_labels=lambda c: c.branch_length)
    return tree


def edit_nested_plain(S1, P1, S2):
    # S1 : séquence primaire
    # P1 : structure imbriquée (sous forme de liste de paires d'indices)
    # S2 : seconde séquence

    n = len(S1)  # longueur de la séquence primaire
    m = len(S2)  # longueur de la seconde séquence

    # Initialisation du tableau de programmation dynamique
    DP = np.full((n+2, n+2, m+2, m+2), np.inf)

    # Paramètres de scoring
    wm = 1  # pénalité de substitution
    wd = 1  # pénalité de suppression
    wa = 1.75  # pénalité de modification d'arc
    wb = 1.5  # pénalité de rupture d'arc
    wr = 2  # pénalité de suppression d'arc

    def mismatch(i, j):
        return 1 if S1[i - 1] != S2[j - 1] else 0

    # Initialisation des cas de base (1)
    for i in range(1, n+1):
        for j in range(1, m+1):
            DP[i, i - 1, j, j - 1] = abs(i - j) * wd

    # Boucle principale
    for k in range(1, n+1):
        for i in range(1, n - k + 2):
            i_prime = i + k - 1
            for j in range(1, m+1):
                for j_prime in range(j, m+1):

                    # Cas 3 : l'arc actuel fait partie de la structure imbriquée
                    if (i, i_prime) in P1:
                        DP[i, i_prime, j, j_prime] = min(
                            DP[i+1, i_prime-1, j+1, j_prime-1] + wb + (mismatch(i, j) + mismatch(i_prime, j_prime)) * wm if j < j_prime else np.inf,
                            DP[i+1, i_prime-1, j, j_prime-1] + wa + mismatch(i_prime, j_prime) * wm,
                            DP[i+1, i_prime-1, j+1, j_prime] + wa + mismatch(i, j) * wm,
                            DP[i+1, i_prime-1, j, j_prime] + wr,
                            DP[i, i_prime, j, j_prime-1] + wd,
                            DP[i, i_prime, j+1, j_prime] + wd,
                        )
                    # Cas 4 : l'arc actuel ne fait pas partie de la structure imbriquée
                    elif i_prime not in [u[1] for u in P1]:
                        DP[i, i_prime, j, j_prime] = min(
                            DP[i, i_prime - 1, j, j_prime - 1] + mismatch(i_prime, j_prime) * wm,
                            DP[i, i_prime - 1, j, j_prime] + wd,
                            DP[i, i_prime, j, j_prime - 1] + wd,
                        )
                    # Cas 4 : l'arc actuel est adjacent à un arc de la structure imbriquée
                    else:
                        min_val = np.inf
                        for j_double_prime in range(j, j_prime + 1):
                            cur_val = DP[i, [u[0] for u in P1 if u[1] == i_prime][0] - 1, j, j_double_prime - 1] + \
                                      DP[[u[0] for u in P1 if u[1] == i_prime][0], i_prime, j_double_prime, j_prime]
                            min_val = min(min_val, cur_val)
                            DP[i, i_prime, j, j_prime] = min_val
                            
    return DP

def backtrace(seq1, seq2, DP):
    def match_score(a, b):
        if a == b:
            return 1
        else:
            return -1

    m = len(seq1)
    n = len(seq2)
    score = [[0 for x in range(n + 1)] for y in range(m + 1)]

    dp = DP[0][0]

    for i in range(m + 1):
        score[i][0] = -i

    for j in range(n + 1):
        score[0][j] = -j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[i - 1], seq2[j - 1])
            delete = score[i - 1][j] - 1
            insert = score[i][j - 1] - 1
            score[i][j] = max(match, delete, insert)

    align1 = ""
    align2 = ""

    i = m
    j = n

    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i - 1][j - 1] + match_score(seq1[i - 1], seq2[j - 1]):
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif i > 0 and score[i][j] == score[i - 1][j] - 1:
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    return align1, align2

def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file