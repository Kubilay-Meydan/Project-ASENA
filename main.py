from Bio import Entrez, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO, Phylo
from Bio.Align.Applications import MuscleCommandline

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

Entrez.email = "your@email.com"  # replace with your email address

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
    with open("output.fasta", "w") as fasta:
        for seq, name in zip(sequences, names):
            fasta.write(">" + name + "\n")
            fasta.write(seq + "\n")



def Align_muscle(sequences,output):
# Align the sequences using muscle
    muscle_exe = "muscle.exe"
    cline = MuscleCommandline(muscle_exe,input=sequences, out=output, clw=True)
    cline()
    alignment = AlignIO.read(output, "clustal")
    return alignment



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

'''


# List of protein accession numbers
accessions = ['P66928','P14949','P10599','P34723','P0A4L3','P0AA25','P08629','P10639','P42115']
# List of protein names 
names = ['THIO_HELPY', 'THIO_BACSU', 'THIO_HUMAN', 'THIO_PENCH', 'THIO_LISMO','THIO_ECOLI', 'THIO_CHICK', 'THIO_MOUSE', 'THIO_NEUCR']


#
#       Gene id for human insulin, gets Exon and Intron details:
#

gene_id = 'L15440.1'
get_genbank_info(gene_id)

#
#       Makes Phylo tree from accession:
#

#   Writes protein sequences to a fasta file using the uniprot accessions above and the names list given by user:
sequences = all_sequences(accessions)
write_fasta(sequences,names)

#   Uses MUSCLE alignement to make a Multiple Sequence Alignement File
Align_muscle('output.fasta','aligned')

#   Makes a simple Phylo tree, no bootstrap yet
make_phylogenetic_tree_bof('aligned')



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