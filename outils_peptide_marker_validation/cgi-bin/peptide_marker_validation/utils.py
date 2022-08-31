
from Bio import Phylo
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Align
from Bio.Align import substitution_matrices
from molmass import *
from peptides import *

import os
import csv
import re

AA = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T'] #List of AA to check if the peptide sequences given are composed of AA


def fasta_format_checker(file_name):
    """Checks if the given file is a valid FASTA format

    Args:
        file_name (str): the name of the file with the sequences in fasta

    Returns:
        bool: whether or not the file is formatted correctly
    """
    #we check before hand if the file is formatted correctly
    valid_response = False
    try:
        with open(file_name, "r") as filin: #we read the file
            l = filin.readline() #we read the 1st line
            if l[0] == ">": #we check if the 1st character is ">" as required in a fasta format
                valid_response = True
                return valid_response
            else:
                print("Your file is not formatted correctly, please try another")
    except FileNotFoundError:
        print("Please verify the file's name")


def read_multi_fasta(file_name):
    """Read a multifasta file and return a dictionary with the sequence and the name of the name of the sequence.

    Args:
        file_name (str): the name of the file with the sequences used.

    Returns:
        dict: dictionnary mapping the name and its sequences
    """
    if fasta_format_checker(file_name):
        seq_dict = {}
        with open(file_name, "r") as filin: #we read the file
            l = filin.readline() #we begin to read the 1st line
            while l != "":
                ligne_entete = l.rstrip() 
                name = ligne_entete[1:] #we define the 1st line as the name of the 1st sequence without the ">"
                if "|" in name:
                    name = name.split("|")[0].strip()
                seq = ""
                s = filin.readline() #we begin to read the 2nd line
                while s != "" and s[0] != ">": 
                    seq += s.rstrip() #we add the line to the sequence
                    s = filin.readline()
                if name not in seq_dict:
                    seq_dict[name] = [seq]
                else:
                    seq_dict[name].append(seq)
                l = s

    return seq_dict

def rna_ro_prot(seq_dict):
    """Takes a dictionnary mapping the name and the RNA sequences and translate it to proteins

    Args:
        seq_dict (dict): dictionnary mapping the name and the RNA sequences from the multifasta sequence input

    Returns:
        dict: dictionnary mapping the name and the sequences translated to proteic sequences
    """
    rna_dict = {}
    for seq_name, liste_seq in seq_dict.items():
        for sequence in liste_seq:
            seq_rna = Seq(sequence)
            seq_prot_not_final = str(seq_rna.translate()) #utilisation de la fonction translate de bio python
            myList = seq_prot_not_final.split("*") #on enlève des "*" des codons stop
            seq_prot = "".join(myList)
            
            if seq_name not in rna_dict:
                rna_dict[seq_name] = [seq_prot]
            else:
                rna_dict[seq_name].append(seq_prot)
    return rna_dict

def dna_to_prot(seq_dict):
    """Takes a dictionnary mapping the name and the DNA sequences and translate it to proteins in 3 frames

    Args:
        seq_dict (dict): dictionnary mapping the name and the DNA sequences
    
    Returns:
        dict: dictionnary mapping the name and the sequences translated to proteic sequences 
        example : name_of_sequence = [[seq1 translated frame 1, seq1 translated frame 2, seq2 translated frame 3], [seq2 translated frame 1, seq2 translated frame 2, seq2 translated frame 3]]
    """
    dna_dict = {}
    for seq_name, liste_seq in seq_dict.items():
        for sequence in liste_seq:
            seq_dna = str(Seq(sequence))
            seq_prot_list = []
            for framestart in range(3): # 3 cadres de lecture
                seq_prot_not_final = translate(seq_dna[framestart:]) #on utilise la fonction translate sur les 3 cadres de lecture
                myList = seq_prot_not_final.split("*") #on enlève les "*" des codons stop
                seq_prot = "".join(myList)
                seq_prot_list.append(seq_prot)
            if seq_name not in dna_dict:
                dna_dict[seq_name] = [seq_prot_list]
            else:
                dna_dict[seq_name].append(seq_prot_list)
    return dna_dict


def prot_mature_alignment(seq, seq_model):
    """Function to map the proteic sequence with the model chosen in order to have the mature sequence of the protein

    Args:
        seq (str): the proteic sequence not mature
        seq_mature (str): the sequence of the model chosen

    Returns:
        str: the mature sequence of the protein
    """
    aligner = Align.PairwiseAligner() #aligner de bio python
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62") #on charge la matrice BLOSUM62
    #paramètres de l'aligner :
    aligner.open_gap_score = -10   
    aligner.extend_gap_score = -0.5
    aligner.target_end_gap_score = -10
    aligner.query_end_gap_score = -0.5
    alignments = aligner.align(seq_model, seq) #retourne une liste des alignements 
    alignment = alignments[0] #on prend le premier (qui est le seul)
    max_len = 0
    for subalignment in list(alignment.aligned[1]):
        len_alignment = subalignment[1]-subalignment[0]
        if len_alignment > max_len: #on regarde le match le plus grand (avec la plus grande séquence d'alignement entre les 2 séquences)
            max_len = len_alignment
            max_len_alignment = subalignment
    mature_prot = seq[max_len_alignment[0]:max_len_alignment[1]]
    return mature_prot


def prot_immature_to_mature(dict_seq, seq_model):
    """Function that takes the dictionnary containing the not mature sequences to return another dictionnary with the mature sequences

    Args:
        dict_seq (dict): the dictionnary containing the not mature sequences, mapping the name and the not mature sequences
        seq_model (str): the model sequence chosen

    Returns:
        dict: dictionnary mapping the name and the not mature sequences
    """
    dict_seq_mature = {}
    for species, list_seq in dict_seq.items():
        for seq in list_seq:
            prot_mature = prot_mature_alignment(seq, seq_model)
            if species not in dict_seq_mature:
                dict_seq_mature[species] = [prot_mature]
            else:
                dict_seq_mature[species].append(prot_mature)
    return dict_seq_mature


def pep_parser_csv(csv_file, threshold = None):
    """CSV parser for the format expected in input

    Args:
        csv_file (str): The path to the csv file.
        threshold (str): Number of decimal digits for the peptide mass.
                         Default to None. (Only used when -m is used).

    Returns:
        list: list of Peptide obj corresponding to the peptide in the csv file.

    """
    _, ext = os.path.splitext(csv_file)
    if ext != ".csv":
        raise NameError("the input file must be a csv file")
    list_pep = []
    f= open (csv_file)
    myReader = csv.reader(f)
    next(myReader)
    for row in myReader:
        if row != []:
            tmp = True
            seq = row[0]
            for elt in seq:
                if elt not in AA:
                    tmp = False
            if tmp :
                if threshold != None:
                    f = Formula('peptide({})'.format(seq)) #utilisation du script pour calculer la masse (script molmass)
                    mass = [format(f.isotope.mass, threshold)]
                else:
                    mass = row[1]
                node = [] #node == nom des séquence où le peptide doit apparaître => pour l'instant vide
                label = row[5]
                if list_pep != []:
                    list_pep_seq = [pep.get_seq() for pep in list_pep]
                    if seq not in list_pep_seq:
                        list_pep.append(Peptides(seq, mass, node, label)) #création des obj peptides
                else:
                    list_pep.append(Peptides(seq, mass, node, label))
    
    return list_pep

def create_fasta(NameOutput, dict_seq):
    """ Create a fasta file with the sequences in dict_seq (with the mature sequences)

    Args:
        NameOutput (str): Name of the output file
        dict_seq (dict): the dictionnary containing the mature sequences, mapping the name and the mature sequences
    """

    list_seq_record = []
    fasta_path = NameOutput + "mature_prot.fasta"
    for species, list_seq in dict_seq.items():
        for seq in list_seq:
            record = SeqRecord( #on créer un obj SeqRecord pour créer un fichier fasta en passant par le module SeqIO de bio python
                Seq(seq),
                id= species,
                description= "|"
                )
            list_seq_record.append(record)
    SeqIO.write(list_seq_record, fasta_path, "fasta")


def pep_parser_rpg(csv_file, threshold=None):
    """Function used to convert RPG csv file to the csv file explained in README.

    Args:
        csv_file (str): The path of the RPG file to convert.
        threshold (str): Number of decimal digits for the peptide mass.
                         Default to None. (Only used when -m is used).
    
    Returns:
        list: List containing sublists where each sublist is the row to insert to the converted CSV file.
    """ 
    _, ext = os.path.splitext(csv_file)
    if ext != ".csv":
        raise NameError("the input file must be a csv file")
    list_pep_to_insert = []
    f= open (csv_file)
    myReader = csv.reader(f)
    next(myReader)
    for row in myReader:
        tmp = True
        seq = row[7]
        for elt in seq:
            if elt not in AA:
                tmp = False
        if tmp :
            if threshold != None:
                f = Formula('peptide({})'.format(seq))
                mass = [format(f.isotope.mass, threshold)]
            else: 
                mass = row[5]
            specie = row[0]
            if "|" in specie:
                specie = specie.split("|")[0].strip()
            position = row[3]
            ptm = None
            label = None
            comment = None
            list_pep_to_insert.append([seq, mass, specie, position, ptm, label, comment])
    return list_pep_to_insert


def pretty_print_peptides_csv_from_rpg(csv_file, NameOutput):
    """Function to write the converted RPG file to the CSV format expected.
    Args:
        csv_file (str): The path to the RPG file to convert.
    """
    list_pep_to_insert = pep_parser_rpg(csv_file)

    csv_file_converted =  NameOutput + "rpg_converted.csv"
    with open(csv_file_converted, 'w', newline='') as results:
        writer = csv.writer(results)
        writer.writerow(["Sequence", "Mass", "Species", "Position", "PTM", "Label", "Comment"])
        for list_to_insert in list_pep_to_insert:
            rowToInsert = list_to_insert
            writer.writerow(rowToInsert)