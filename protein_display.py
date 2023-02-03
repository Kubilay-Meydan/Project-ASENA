
import requests
import re
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import IPythonConsole

def display_protein(uniprot_id):
    # Use UniProt's RESTful API to retrieve the PDB ID associated with the UniProt identifier
    response = requests.get("https://www.uniprot.org/uniprot/" + uniprot_id + ".txt")
    txt = response.text
    pdb_id = re.search("PDB; (.+?);", txt).group(1)

    # Download the PDB file for the protein structure
    response = requests.get("https://files.rcsb.org/download/" + pdb_id + ".pdb")
    pdb_file = response.text

    # Use RDKit to display the protein structure
    mol = Chem.MolFromPDBBlock(pdb_file, removeHs=False)
    rdDepictor.Compute2DCoords(mol)
    IPythonConsole.molToPNG(mol)

# Example usage
display_protein("P12345")