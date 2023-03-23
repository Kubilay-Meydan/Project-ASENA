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


def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file