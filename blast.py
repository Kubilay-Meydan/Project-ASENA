from Bio.Blast import NCBIWWW

def blast(sequence, database='nt', program='blastn', num_alignments=10):
    """
    Perform a BLAST search using the given sequence and parameters.
    
    sequence: the query sequence to search.
    database: the name of the database to search (default: 'nt').
    program: the type of BLAST program to use (default: 'blastn').
    num_alignments: the number of alignments to return (default: 10).
    """
    result_handle = NCBIWWW.qblast(program, database, sequence, hitlist_size=num_alignments)
    return result_handle.read()
result = blast("ATGGCCATG")
print(result)


