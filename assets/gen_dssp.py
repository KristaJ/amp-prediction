from Bio.PDB import PDBParser, PDBList
from Bio.PDB.DSSP import DSSP
import warnings
import re
import requests
from Bio import SeqIO
import os.path



def getFastaFromPDB(ID=None, filename=None, chain=None):
    if ID:
        PDBFile = './PDB_FILES/pdb{}.ent'.format(ID.lower()[:4])
        if not os.path.exists(PDBFile):
            pdbl = PDBList()
            pdbl.retrieve_pdb_file(ID, file_format = 'pdb', pdir ='./PDB_FILES' )
            if not os.path.exists(PDBFile):
                return (1, "PDB entry not found")
    elif filename:
        PDBFile = filename
    fasta_dict = {}
    for record in SeqIO.parse(PDBFile, 'pdb-atom'):
        pdbid = record.id[:4]
        chain_key = record.id[-1]
        fasta_dict[chain_key] = record.seq._data
    try:
        return (0, fasta_dict[chain])
    except KeyError:    
        if (chain == "") & (len(fasta_dict.keys()) > 1):
            return_string = "This entry contains multiple chains.  Please select one to continue\n"
            for key in fasta_dict.keys():
                return_string += f'\nChain {key}:\t{fasta_dict[key]}'
            return_string += '\n Please select a chain to continue'
            return (1,return_string)
        elif (chain == ""):
            chain = list(fasta_dict.keys())[0]
        else:
            return (1,"Selected chain is not present in this PDB file")


def getDSSPfromPDB(ID=None, filename=None, chain=None):
    if chain == "":
        chain = 'A'
    Helix = ['H', 'G', 'I']
    Strand = ['E', 'B']
    Loop = ['S', 'T', 'C', '-']
    
    p = PDBParser()
    
    if ID:
        PDBFile = './PDB_FILES/pdb{}.ent'.format(ID.lower()[:4])
        if not os.path.exists(PDBFile):
            return "PDB entry not found"
    if filename:
        PDBFile=filename
    structure = p.get_structure(ID, PDBFile)
    model = structure[0]
    dssp = DSSP(model, PDBFile, dssp='mkdssp')
    
    keys = dssp.keys()
    struct = []
    sequence = []
    for key in keys:
        if key[0] == chain:
            struct.append(dssp[key][2])
    struct = ['H' if (x in Helix) else x for x in struct]
    struct = ['B' if (x in Strand) else x for x in struct]
    struct = ['C' if (x in Loop) else x for x in struct]
    struct = ''.join(struct)
    return struct
