import requests
import Bio
import matplotlib
from Bio import Entrez, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO
from Bio import Align, AlignIO, Phylo
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
''''
entry = ''
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
    return ans
DNA_to_RNA('AATTTTTGGGGGCCCCC')

Entrez.email = "theicebreakerofficial@gmail.com"  # replace with your email address

def get_genbank_info(gene_id):
    with open(str(gene_id)+'.txt', 'w') as f:
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

gene_id = 'L15440.1'
#get_genbank_info(gene_id)

# List of protein accession numbers
accessions = ['P66928','P14949']

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

sequences = []
for e in accessions:
    sequences.append(str(get_sequence(e)))

'''
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
def Align_muscle(sequences,output):
# Align the sequences using muscle
    muscle_exe = "muscle.exe"
    cline = MuscleCommandline(muscle_exe,input=sequences, out=output, clw=True)
    cline()
    alignment = AlignIO.read(output, "clustal")
    return alignment
Align_muscle('sequences.fasta','aligned')


from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

def make_phylogenetic_tree_bof(alignment_file):
    # Load the multiple sequence alignment from file
    alignment = AlignIO.read(alignment_file, "clustal")
    # Calculate the pairwise distances between sequences
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    # Construct the phylogenetic tree using the UPGMA method
    constructor = DistanceTreeConstructor(calculator, 'upgma')
    tree = constructor.build_tree(alignment)
    # Draw and show the tree
    Phylo.draw(tree)

make_phylogenetic_tree_bof('aligned')


'''






Bootsrap???




from Bio import AlignIO
from dendropy import Tree
from dendropy.math import euclidean_distance

def make_phylogenetic_tree(alignment_file, bootstrap_replicates=100):
    # Load the multiple sequence alignment from file
    alignment = AlignIO.read(alignment_file, "clustal")
    # Compute the pairwise distances between sequences in the alignment
    distances = []
    for i, seq1 in enumerate(alignment):
        for j in range(i + 1, len(alignment)):
            seq2 = alignment[j]
            distances.append(euclidean_distance(seq1.seq, seq2.seq))
    # Construct the phylogenetic tree using the UPGMA method
    tree = Tree.from_matrix(distances, "upgma")
    # Apply bootstrapping to the tree
    bootstrap_trees = tree.bootstrap_trees(replicates=bootstrap_replicates)
    # Create a consensus tree from the bootstrapped trees
    consensus_tree = tree.consensus(min_freq=0.5)
    # Draw and show the tree
    consensus_tree.print_plot()



make_phylogenetic_tree('aligned')


'''