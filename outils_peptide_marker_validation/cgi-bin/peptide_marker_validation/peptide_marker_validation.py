import csv
from Bio import Align
from decimal import *

from utils import *



def compare_pep_rpg_pep_input(list_pep_rpg, list_pep_input, dict_seq, threshold = None):
    """Function to compare the peptides given by rpg and the peptides given in input

    Args:
        list_pep_rpg (list): list containing the peptides given by rpg
        list_pep_input (list): list containing the peptides given in input
        dict_seq (dict): dictionnary mapping the name and the mature sequences 
        threshold (str, optional): the threshold chosen only of we want to compare the mass. Defaults to None.

    Returns:
        dict: dictionnary mapping the input peptides and the peptides given by rpg that are equal (mass or sequence)
    """

    dict_pep_input = {}

    list_pep_rpg_complete = search_pep_in_seq(list_pep_rpg, dict_seq, threshold) #fonction pour savoir où les peptides rpg apparaissent 
    for pep_input in list_pep_input:
        list_pep_rpg_equal_pep_input = []
        for pep_rpg in list_pep_rpg_complete:
            if threshold != None: #si on regarde la masse
                for mass_pep_input in pep_input.get_mass():
                    for mass_pep_rpg in pep_rpg.get_mass():
                        if mass_pep_input == mass_pep_rpg:
                            if pep_rpg not in list_pep_rpg_equal_pep_input:
                                list_pep_rpg_equal_pep_input.append(pep_rpg)
            elif pep_input.get_seq() == pep_rpg.get_seq(): #sinon on regarde la séquence
                list_pep_rpg_equal_pep_input.append(pep_rpg)
        if pep_input not in dict_pep_input:
            dict_pep_input[pep_input] = list_pep_rpg_equal_pep_input
        else:
            print("this peptide already exists : {}".format(pep_input.get_seq()))
    

    for pep_input, list_pep_rpg_equal in dict_pep_input.items():
        list_pep_nodes = [] #les noms des espèces 
        list_pep_position = []
        list_target_seq = []
        list_target_mass = []
        for pep_rpg in list_pep_rpg_equal:
            for nodes in pep_rpg.get_nodes():
                if nodes not in list_pep_nodes:
                    list_pep_nodes.append(nodes)
            for positions in pep_rpg.get_positions():
                if positions not in list_pep_position:
                    list_pep_position.append(positions)
            for target_seq in pep_rpg.get_target_seq():
                if target_seq not in list_target_seq:
                    list_target_seq.append(target_seq)
            for target_mass in pep_rpg.get_target_mass():
                if target_mass not in list_target_mass:
                    list_target_mass.append(target_mass)
        pep_input.set_nodes(list_pep_nodes) 
        pep_input.set_positions(list_pep_position)
        pep_input.set_target_seq(list_target_seq)
        pep_input.set_target_mass(list_target_mass)
    
    return dict_pep_input




def search_pep_in_seq(list_pep, dict_seq, threshold= None):
    """Function to search peptides in different proteic sequence.

    Args:
        list_pep (list): List of peptides.
        dict_seq (dict): dict of proteic sequences. Keys = name of proteic sequence, values = proteic sequence.

    Returns:
        list: list with Peptide obj with their attribute "nodes", "positions", "target_seq" and "target_mass" filled.

    """
    dict_pep = {}
    dict_indices = {}

    for pep in list_pep:
        pep_seq = pep.get_seq()
        liste_species = []
        for specie, list_seq in dict_seq.items():
            for seq in list_seq:
                match = re.search(pep_seq, seq)
                if match != None:
                    liste_species.append(specie)
                    for match in re.finditer(pep_seq, seq):
                        if threshold == None: #on calcule la masse du target par défault pour l'afficher dans les résultat même si on ne compare pas les masses
                            f = Formula('peptide({})'.format(seq[match.start():match.end()]))
                            mass_match = [format(f.isotope.mass, ".3f")] #pas utile
                        else:
                            f = Formula('peptide({})'.format(seq[match.start():match.end()]))
                            mass_match = [format(f.isotope.mass, threshold)]
                        if pep not in dict_indices:
                            dict_indices[pep] = [((match.start(), match.end()), specie, (seq[match.start():match.end()], specie), (mass_match, specie))] 
                        else:
                            dict_indices[pep].append(((match.start(), match.end()), specie, (seq[match.start():match.end()], specie), (mass_match, specie)))

        if liste_species == []:
            print("This pep is not present in any of the sequences : {}".format(pep_seq))
        if pep not in dict_pep:
            dict_pep[pep] = liste_species
        else:
            print("This peptide already exists :{}".format(pep_seq))
    
    for pep, species in dict_pep.items():
        pep.set_nodes(species) # on met en mémoire les noms des séquences dans l'attribut nodes 
    
    for pep, list_match in dict_indices.items():
        liste_position = []
        liste_target_seq = []
        liste_target_mass = []
        for tuple_match in list_match:
            liste_position.append((tuple_match[0], tuple_match[1])) #pour mettre en mémoire dans l'attribut position les positions de début et de fin et le nom de la séquence
            liste_target_seq.append(tuple_match[2]) #pour mettre en mémoire dans l'attribut target_seq la séquence target et le nom de la séquence target
            liste_target_mass.append(tuple_match[3]) #pour mettre en mémoire dans l'attribut target_mass la masse des séquences target et le nom de la séquence target
        pep.set_positions(liste_position)
        pep.set_target_seq(liste_target_seq)
        pep.set_target_mass(liste_target_mass)

            
    return list(dict_pep.keys())


#Function used when -m is used
#not final
def search_pep_in_seq_with_mass(list_pep, dict_seq, threshold):
    """Function to see if the mass of the peptides appears in the proteic sequence.

    Args:
        list_pep (list): List of peptides.
        dict_seq (dict): dict of proteic sequence. Keys = name of proteic sequence, values = proteic sequence.
        threshold (str): Number of decimal digits for the peptide mass.

    Returns:
        tuple: First : a dict with keys = Peptide and values = species they appear in, 
               second : a dict with keys = Peptide and values = species they appear in ONLY MASS WIZE
    """
    dict_pep = {}
    dict_pep_mass_only = {}
    
    for pep in list_pep: #for each peptide
        pep_seq = pep.get_seq()
        pep_mass = pep.get_mass()
        list_species_mass_seq = [] #the peptides that match sequence and mass wize
        list_species_mass_only = [] #the peptides that only match mass wize
        for mass in pep_mass: #for each mass (it could include ptm)
            substring_unmatched = [] #already seen substrings that didn't match the peptide mass (we can pass them)
            for specie, list_seq in dict_seq.items(): 
                for seq in list_seq: #for each sequence
                    for i in range(len(seq)):
                        for j in range(i+1, len(seq)+1):
                            substring = seq[i:j]
                            tmp = True
                            for elt in substring:
                                if elt not in AA: #We check if the substring is made only of AA
                                    tmp = False #If not we pass it
                            if tmp :
                                if substring in substring_unmatched:
                                    continue
                                else:
                                    f = Formula('peptide({})'.format(substring))
                                    mass_substring = format(f.isotope.mass, threshold)
                                    if Decimal(mass_substring).compare(Decimal(mass)) == Decimal("0"): #if mass_subtring = pep_mass
                                        if pep_seq == substring:
                                            list_species_mass_seq.append(specie)
                                        else:
                                            list_species_mass_only.append(specie)
                                    elif Decimal(mass_substring).compare(Decimal(mass)) == Decimal("1"): #if mass_substring > pep_mass
                                        substring_unmatched.append(substring)
                                        break
                                    else:
                                        substring_unmatched.append(substring)

        if list_species_mass_seq == [] and list_species_mass_only == []:
            print("This peptide is not found in any of the sequences : {}".format(pep_seq))
        
        if pep not in dict_pep:
            list_all_species = list_species_mass_seq + list_species_mass_only
            dict_pep[pep] = list(dict.fromkeys(list_all_species)) #remove duplicate species
        
        if pep not in dict_pep_mass_only:
            dict_pep_mass_only[pep] = list_species_mass_only
    
    for pep, species in dict_pep.items():
        pep.set_nodes(species)
    
    return dict_pep, dict_pep_mass_only


def pretty_print_peptides_csv(list_pep_with_species, output_file):
    """Function to generate the csv file with the peptides.

    Args:
        list_pep_with_species (list): list of peptides obj.
        output_file (str): the name of the output file (with path). 
    """
    with open(output_file, 'w', newline='') as results:
        writer = csv.writer(results)
        writer.writerow(["Sequence", "Mass", "Nodes", "Position", "Label", "Target_seq", "Target_mass"])
        for peptide in list_pep_with_species:
            node_liste = [specie for specie in peptide.get_nodes()]
            if len(peptide.get_mass()) == 1:
                rowToInsert = [peptide.get_seq(), peptide.get_mass()[0], node_liste, peptide.get_positions(), peptide.get_label(), peptide.get_target_seq(), peptide.get_target_mass()]
            else:
                rowToInsert = [peptide.get_seq(), peptide.get_mass(), node_liste, peptide.get_positions(), peptide.get_label(), peptide.get_target_seq(), peptide.get_target_mass()]
            writer.writerow(rowToInsert)

#if -m is used
def pretty_print_peptides_csv_with_mass(dict_pep_with_species, dict_pep_with_species_mass_only, output_file):
    """Function to generate the csv file with the peptides and their species.

    Args:
        dict_pep_with_species (dict): keys = Peptide, values = list of species they appear in.
        dict_pep_with_species_mass_only (dict): keys = Peptide, values = list of species they appear in MASS WIZE ONLY.
        output_file (str): the name of the output file (with path). 
    """
    with open(output_file, 'w', newline='') as results:
        writer = csv.writer(results)
        writer.writerow(["Sequence", "Mass", "Nodes", "Comment"])
        for peptide in list(dict_pep_with_species.keys()):
            node_liste = [specie for specie in peptide.get_nodes()]
            if dict_pep_with_species_mass_only[peptide] != []:
                if len(peptide.get_mass()) == 1:
                    rowToInsert = [peptide.get_seq(), peptide.get_mass()[0], node_liste, "warning: mass only for {}".format(dict_pep_with_species_mass_only[peptide])]
                else:
                    rowToInsert = [peptide.get_seq(), peptide.get_mass(), node_liste, "warning: mass only for {}".format(dict_pep_with_species_mass_only[peptide])]
            else:
                if len(peptide.get_mass()) == 1:
                    rowToInsert = [peptide.get_seq(), peptide.get_mass()[0], node_liste, ""]
                else:
                    rowToInsert = [peptide.get_seq(), peptide.get_mass(), node_liste, ""]
            writer.writerow(rowToInsert)

def pretty_print_pep_rpg_equal_pep_input(dict_pep_input_equal_pep_rpg, name_file):
    """fonction to generate the csv result file in the case when RPG is used

    Args:
        dict_pep_input_equal_pep_rpg (dict): dict that maps the peptide of the input and the list of peptides from rpg that are equal (mass or sequence)
        name_file (str): the name of the csv file
    """

    with open(name_file, "w", newline="") as results:
        writer = csv.writer(results)
        writer.writerow(["Sequence", "Mass", "Nodes", "Position", "Label", "Target_seq", "Target_mass"])

        for pep_input, list_pep_rpg_equal in dict_pep_input_equal_pep_rpg.items():
            if len(pep_input.get_mass()) == 1: #dans le cas où pas de PTM
                rowToInsert = [pep_input.get_seq(), pep_input.get_mass()[0], pep_input.get_nodes(), pep_input.get_positions(), pep_input.get_label(), pep_input.get_target_seq(), pep_input.get_target_mass()] #get_nodes = les espèces
            else: #utilité de "pep_input.get_mass()" soit une liste pour gérer les PTM
                rowToInsert = [pep_input.get_seq(), pep_input.get_mass(), pep_input.get_nodes(), pep_input.get_positions(), pep_input.get_label(), pep_input.get_target_seq(), pep_input.get_target_mass()]
            writer.writerow(rowToInsert)