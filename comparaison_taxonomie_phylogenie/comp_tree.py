from Bio import Phylo
import os
import pandas as pd
import csv
import argparse
acides_amines = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']


def remove_not_common_clades(tree1, tree2):
    """Remove clades from tree1 not present in tree2 .

    Args:
        tree1 (Tree): tree with the species that needs to be removed
        tree2 (Tree): the tree we're doing the comparison with
    """
    for clade in tree1.get_terminals():  # On regarde toutes les espèces (feuilles) de l'arbre 1
        if not tree2.find_any(clade.name): #si l'espèce n'est pas retrouver dans l'arbre 2
            tree1.collapse(clade) # on supprime la branche correspondant à cette espèce dans l'arbre 1


def get_dico_clade_non_terminal(tree_taxo):
    """Returns a dictionary mapping all the non terminal nodes to their child leaves (species) in the taxonomy.

    Args:
        tree_taxo (Tree): taxonomic tree

    Returns:
        dict: dictionary mapping all the non terminal nodes to their child leaves (species) in the taxonomy
    """
    dico_clade_taxo = {} # le dictionnaire qui va mapper les clades (les noeud de l'arbre taxonomique) avec leur espèces
    intern_nodes = tree_taxo.get_nonterminals() #on ne va prendre que les noeuds de l'arbre
    intern_nodes += tree_taxo.find_clades("Eutheria") # on ajoute le noeud correspondan à la racine
    for clade in intern_nodes:
        node_childs = clade.clades #on regarde les enfants du sous arbre correspondant à la clade dans l'arbre taxonomique
        for child in node_childs:
            if not child.is_terminal(): #si le fils n'est pas une feuille
                node_childs += get_leaf_nodes(child) #on ajoute les feuilles de ce sous arbre fils
        dico_clade_taxo[clade.name] = list(dict.fromkeys([x for x in node_childs if x.is_terminal()])) #on mappe les feuilles à la clade
    return dico_clade_taxo

def get_clade_mrca_phylo(tree_phylo, dico_clade_taxo):
    """Get the mrca (most recent common ancestor) subtree in the phylogenetic tree for all species of a clade (a non terminal node). For all clade in the taxonomy.

    Args:
        tree_phylo (Tree): the phylogenetic tree
        dico_clade_taxo (dict): the dictionary mapping all the non terminal nodes to their child leaves (species) in the taxonomy

    Returns:
        dict: dictionary mapping all the non terminal node (clade) from the taxonomy with the mrca subtree (the mrca of all the species of the clade in the taxonomy) in the phylogenetic tree.
    """
    dico_mrca = {} # Le dictionnaire mappant la clade dans la taxonomie avec les espèces de l'arbre qui représente le dernier ancêtre commun
    for clade in dico_clade_taxo: #pour toutes les clades présents dans la taxonomie
        liste_mrca = []
        mrca = tree_phylo.common_ancestor([x.name for x in dico_clade_taxo[clade]]) # most recent common ancestors => plus petit sous arbres commun avec toutes les espèces de la clade
        if mrca.clades==[]: #si le dernier ancêtre commun est une feuille
            liste_mrca.append(mrca) #on ajoute la feuille dans le dictionnaire
            dico_mrca[clade] = liste_mrca
        else:
            for node in mrca.clades: #pour les fils du sous arbre (dernier ancêtre commun)
                if not node.is_terminal():
                    mrca.clades += get_leaf_nodes(node) #on ajoute les feuilles des sous arbres (des noeuds fils)
            dico_mrca[clade] = list(dict.fromkeys([x for x in mrca.clades if x.is_terminal()])) #on mappe les feuilles à la clades 
    return dico_mrca

def get_leaf_nodes(node):
    """Get all leaf nodes of a node (clade).

    Args:
        node (Clade): the clade (subtree) we want to have all the leaves from

    Returns:
        list: list containing all the leaves of the clade (subtree)
    """
    leaves = []
    collect_leaf_nodes(node,leaves)
    return leaves #la liste avec toutes les feuilles du sous arbre

def collect_leaf_nodes(clade, leaves):
    """Collect all the leaves from a clade recursively

    Args:
        clade (Clade): the clade we want to collect the leaves from
        leaves (list): the list containing all the leaves collected from the clade
    """
    if clade is not None:
        if len(clade.clades) == 0:
            leaves.append(clade) #on ajoute l'object clade dans la liste si c'est une feuille
        for n in clade.clades:
            collect_leaf_nodes(n, leaves) #on continue dans les sous arbre fils recursivement


def count_misplaced_clade(dico_clade_taxo, dico_mrca, species_grouped, tree_phylo):
    """Counts the misplaced clade in the phylogenetic tree and gives a percentage of clades without problems. Compare the species of a clade in the taxonomy with the mrca subtree of a phylogenetic tree. 
    Only for problematic clades, prints the ascii representation of the mrca, prints also the correct subtrees found in the mrca, gives the distance of the subtrees to the root of the mrca tree. 
    Also gives the name of the species missplaced in the phylogenetic tree. Gives the number of species in the taxonomy for the clade and the number of species in the mrca too.

    Args:
        dico_clade_taxo (dict): the dictionary mapping all the non terminal nodes to their child leaves (species) in the taxonomy
        dico_mrca (dict): the dictionary mapping all the non terminal nodes (clades) from the taxonomy with the mrca subtree (the mrca of all the species of the clade in the taxonomy) in the phylogenetic tree

    Returns:
        tuple: a tuple containing the number of misplaced clade and the list of missplaced clade
    """
    #dictionnaire pour avoir chaque clade mappé avec une liste de sous arbres du dernier ancêtre commun cohérent avec la taxonomie
    dico_smallest_tree = get_smallest_good_tree(species_grouped, dico_test_taxo, tree_phylo) 

    count_false=0
    clade_misplaced = []
    for clade in dico_clade_taxo: #Pour toutes les clades de la taxonomie
        liste_taxo = sorted([x.name for x in dico_clade_taxo[clade]]) #liste de toutes les espèces de cette clade dans la taxonomie
        liste_phylo = sorted([x.name for x in dico_mrca[clade]]) #liste de toutes les espèces du dernier ancêtre commun des espèces de la clade dans l'arbre phylogénétique
        if liste_taxo != liste_phylo: #s'il y a une différence entre ces listes alors il y a des espèces mal placées
            count_false += 1 #on ajoute 1 au conteur des clades problématiques
            clade_misplaced.append(clade) #on ajoute la clade dans la liste des clades avec un problème
            title = "|"+"For clade {} :".format(clade) + "|"
            for i in range(len(title)):
                tiret = "-"*i
            print(tiret)
            print("|"+"For clade {} :".format(clade) + "|")
            print(tiret)

            print("Nb of species in the taxonomy : {}".format(len(dico_clade_taxo[clade])))
            print("Nb of species in the LCA found in the phylogenetic tree : {}".format(len(dico_mrca[clade])))
            species_missplaced = set(liste_taxo) ^ set(liste_phylo) #on fait un XOR sur les 2 listes d'espèces pour avoir la différence
            print("Nb of misplaced species : {}".format(len(list(species_missplaced))))
            print("Misplaced species:")
            print(list(species_missplaced))
            
            mrca = tree_phylo.common_ancestor([x.name for x in dico_clade_taxo[clade]]) #on retrouve le dernier ancêtre commun dans l'arbre phylogénétique
            print("LCA in ascii :")
            Phylo.draw_ascii(mrca) #print en ascii du dernier ancêtre commun

            
            nb_subtree = 0 #compteur de sous arbre
            liste_distance = []
            for subtree in dico_smallest_tree[clade]: #pour chaque sous arbre du dernier ancêtre commun cohérents avec la taxonomie 
                liste_distance.append(int(mrca.depths(unit_branch_lengths=True)[subtree])) #la méthode pour calculer la distance du sous arbre avec la racine du dernier ancêtre commun
                nb_subtree += 1
                try: #il se peut que le print ascii ne marche pas
                    print("Subtree nb {}:".format(nb_subtree))
                    Phylo.draw_ascii(subtree) #méthode pour print en ascii
                except: #si c'est la cas on peut quand même visualiser l'arbre
                    if subtree.clades == []: #si le sous arbre est une feuille
                        print(subtree) #on affiche alors juste la feuille
                    else: #on utilise ici matplotlib si la 1ere méthode échoue
                        print("matplotlib figure")
                        Phylo.draw(subtree, title={"label" : "Clade {}, subtree {}".format(clade, nb_subtree)})
    
            print("Nb of subtrees taxonomically correct: {}".format(nb_subtree))

            for i in range(nb_subtree):
                print("Distance of subtree {} from LCA's root : {}".format(i+1, liste_distance[i]))

            print("=========================================")
    if count_false == 0:
        print("All LCA trees match with the clade in the taxonomy")
        print("=========================================")

    pourcentage = ((len(dico_clade_taxo)-count_false)/len(dico_clade_taxo))*100 #calcule du pourcentage de clade sans aucun problème
    print("Percentage of clade in the taxonomy that was correctly found in the phylogenetic tree : {}%".format(pourcentage))
    print("=========================================")

    return count_false, clade_misplaced


def terminal_child_node_duplicate_remover(mrca):
    """Remove duplicate species sometimes found in the list generated by ".clades" attribute

    Args:
        mrca (Clade): The clade corresponding to the mrca

    Returns:
        list: the childs of the mrca but without duplicate leaves (species)
    """
    liste_child = []
    for node in mrca.clades: #pour tous les fils trouvés dans l'attribut ".clades"
        if not node.is_terminal(): #si c'est un noeud
            liste_child += get_leaf_nodes(node) #on ajoute toutes les feuilles de ce sous arbre
    for node in [x for x in mrca.clades if x.is_terminal()]: #pour tous les fils qui sont des feuilles
        if node in liste_child: #s'il est déjà présent dans la liste des fils
            mrca.clades.remove(node)  #on l'enlève
    return mrca.clades


def get_species_grouped(dico_clade_taxo, tree_phylo):
    """Get species grouped by subtree. For a clade of the taxonomy, map the mrca found in the phylogenetic tree. But this time, the species are grouped in sublists
    corresponding to their subtree.

    Args:
        dico_clade_taxo (dict): the dictionary mapping the clade's name to its species in the taxonomy

    Returns:
        dict: the dictionnary with keys being the clade's name and the values are sublists corresponding to the species found in the mrca subtree.
        Each sublist correspond to the leaves of a subtree
    """
    dico_smallest_sub_tree = {} #le dictionnaire mappant les clades avec leur espèces mais cette fois ci groupées par sous arbres
    for clade in dico_clade_taxo:
        liste_mrca = []
        # On retrouve le dernier ancêtre commun des espèces de la clade dans l'arbre phylogénétique
        mrca = tree_phylo.common_ancestor([x.name for x in dico_clade_taxo[clade]]) 

        mrca_clade = terminal_child_node_duplicate_remover(mrca) #on enlèce les espèces doublon qu'on peut trouver dans les fils du dernier ancêtre commun

        if mrca_clade==[]: #si le dernier ancêtre commun de la clade est une feuille
            liste_mrca.append([mrca])   #on l'ajoute aux dictionnaire
            dico_smallest_sub_tree[clade] = liste_mrca
        else:
            for node in mrca_clade:
                if not node.is_terminal():
                    # on ajoute dans la liste les sous listes correspondant aux espèces des sous arbres du dernier ancêtre commun
                    liste_mrca.append(list(dict.fromkeys(get_leaf_nodes(node)))) 
                else:
                    liste_mrca.append([node])

            #si toutes les sous listes ont une longueur de 1, alors ce sont les fils d'un noeud pré terminal, on peut alors les fusionner les sous listes        
            if all(len(sublist) == 1 for sublist in liste_mrca): 
                total = []
                for sublist in liste_mrca:
                    total += sublist
                liste_mrca = total

            dico_smallest_sub_tree[clade] = liste_mrca #on ajoute la liste au dictionnaire

    return dico_smallest_sub_tree

def get_smallest_good_tree(dico_species_grouped, dico_clade_taxo, tree_phylo):
    """Returns the smallest subtree of a MRCA that is taxonomicaly correct (every leaves of the subtree is a specie corresponding to a clade, a subtree without species that do not belong to the clade)

    Args:
        dico_species_grouped (dict): a dictionnary mapping the clade with sublists corresponding to the leaves of the 2 subtrees (the species grouped by subtrees)
        dico_clade_taxo (dict): the dictionnary mapping the clade with the species of this clade in the taxonomy
        tree_phylo (Tree): the tree object corresponding to the phylogenetic tree

    Returns:
        dict: a dictionnary mapping the clade with a list containing Tree objects corresponding to the smallest correct subtrees (subtrees containing only species of the clade found in the taxonomy)
    """
    dico_sublist = {} # le dictionnaire mappant les clades avec les plus petits sous arbres du dernier ancêtre commun qui sont correct par rapport à la taxonomie
    dico_clade_taxo_copy = dico_clade_taxo.copy() #on fait ici une copie du dictionnaire de base...
    for clade in dico_clade_taxo_copy:
        dico_clade_taxo_copy[clade] = [x.name for x in dico_clade_taxo_copy[clade]] #...pour pouvoir remplacer les objects Clade par les noms d'espèces qu'ils représentent
        
    for clade in dico_species_grouped:
        dico_sublist[clade] = []
        if type(dico_species_grouped[clade][0]) == list: #on vérifie que le 1er object dans la liste est bien une liste aussi (si on regarde bien des feuilles d'un sous arbres du dernier ancêtre commun)
            for subtree in dico_species_grouped[clade]:
                if subtree != []:
                    the_good_subtree = []
                    mrca_subtree = tree_phylo.common_ancestor([x for x in subtree]) #on retrouve le dernier ancêtre commun pour ces espèces de la clade dans l'arbre phylogénétique
                    get_smallest_child_tree(mrca_subtree, dico_clade_taxo_copy[clade], the_good_subtree) # la fonction qui permet de trouver les sous arbres cohérents avec la taxonomie récursivement
                    dico_sublist[clade] += the_good_subtree #on ajoute les bon sous arbres dans une liste qui est elle même ajoutée au dictionnaire
        else:
            #si le 1er object n'est pas une liste (le dernier ancêtre commun est ici un noeud pré terminal et ses fils sont des feuilles et pas des sous arbres)
            dico_sublist[clade].append(tree_phylo.common_ancestor([i for i in dico_species_grouped[clade] if i.name in dico_clade_taxo_copy[clade]])) # on ajoute seulement les espèces présentes dans la taxonomie
    return dico_sublist

def get_smallest_child_tree(mrca_subtree, dico_taxo_clade, the_good_subtree):
    """This is the function that allows to see if the MRCA subtree contains species that do not belong to the clade. If all the species belong, it takes the Tree object corresponding to the subtree et appends
    it to a list. Then recursively it checks the childs subtrees until all the species of the subtree are correct. This list is afterward passed to get_smallest_good_tree() to map it with a clade in the taxonomy.

    Args:
        mrca_subtree (Clade): the subtree corresponding to the MRCA in the phylogenetic tree
        dico_taxo_clade (dict): the dictionnary mapping all the clades with their species in the taxonomy
        the_good_subtree (list): list containing the subtrees of the MRCA that are taxonomy correct that is recursively filled (in the beginning the list is empty)

    Returns:
        list: list containing the subtrees of the MRCA that are taxonomy correct
    """
    def is_one_leaf_not_ok(dico_taxo_clade, childs):
        """Check if at least one of the leaf corresponds to a specie that doesn't belong to the clade.

        Args:
            dico_taxo_clade (dict): the dictionnary mapping all the clade with their species found in the taxonomy (this time the species are in the form of strings representing their name)
            childs (list): a list containing all the leaves of the MRCA. All the species found in the MRCA in the phylogenetic tree

        Returns:
            bool: whether if there is a specie that do not belong to the clade or not
        """
        if any(leaf.name not in dico_taxo_clade for leaf in childs): 
            return True
        else:
            return False

    childs = get_leaf_nodes(mrca_subtree) #on trouve les feuilles de ce sous arbres
    if mrca_subtree != None:
        if childs != []:
            if is_one_leaf_not_ok(dico_taxo_clade, childs): #si au moins une espèce du sous arbre du dernier ancêtre commun n'est pas retrouvée dans la taxonomie
                liste_child_mrca = terminal_child_node_duplicate_remover(mrca_subtree) #on regarde les sous arbres fils
                if liste_child_mrca != []:
                    for n in liste_child_mrca:
                        get_smallest_child_tree(n, dico_taxo_clade, the_good_subtree) #on revérifie la même condition récursivement sous les sous arbres
                
            else: #si toutes les espèces sont cohérentes avec la taxonomie on peut l'ajouter dans une liste
                the_good_subtree.append(mrca_subtree) 
                return the_good_subtree

def fasta_format_checker(file_name):
    """Checks if the given file is a valid FASTA format

    Args:
        file_name (str): the name of the file with the proteic sequences used to build the phylogenetic tree in fasta format

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
    """Read a multifasta file and return a dictionary with the sequence and the name of the species.

    Args:
        file_name (str): the name of the file with the proteic sequences used to build the phylogenetic tree in fasta format

    Returns:
        dict: dictionnary mapping the name of the specie and its proteic sequence
    """
    if fasta_format_checker(file_name):
        seq_dict = {}
        with open(file_name, "r") as filin: #we read the file
            l = filin.readline() #we begin to read the 1st line
            while l != "":
                ligne_entete = l.rstrip() 
                name = ligne_entete[1:] #we define the 1st line as the name of the 1st sequence without the ">"
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

def newick_to_xml_with_seq(seq_dict, tree_taxo):
    """This function converts a Tree in Newick to a phyloXML tree and links the proteic sequence to the clade in the tree.

    Args:
        seq_dict (dict): dictionnary mapping the name of the specie and its proteic sequence
        tree_taxo (Tree): the Tree object build from the taxonomic tree in Newick format

    Raises:
        Exception: if a specie in the tree doesn't have a proteic sequence linked with it

    Returns:
        tuple: (Tree object in phyloxml, name of the phyloxml file)
    """
    tree_taxo_xml = tree_taxo.as_phyloxml() #converti l'object Tree fait en Newick en phyloxml
    for clade in tree_taxo_xml.get_terminals(): #on prend toutes les espèces (les feuilles de l'arbre)
        key = clade.name
        if key in seq_dict: #si une espèce dans l'arbre se retrouve dans le dictionnaire (et donc a une séquence protéique)
            mol_seq = Phylo.PhyloXML.MolSeq(seq_dict[key][0]) #on crée un object MolSeq à partir de la séquence protéique
            sequence = Phylo.PhyloXML.Sequence(mol_seq=mol_seq) #on crée un object Sequence qui a pour attribut l'objet MolSeq
            clade.sequences.append(sequence) #on ajoute l'object Sequence dans l'objet Clade de l'arbre phyloxml
        else:
            raise Exception("The species {} don't have a sequence".format(clade.name)) #prévient si une espèce de l'arbre n'a pas de séquence protéique associée
    #On sauvergarde l'arbre phyloxml dans un fichier
    name = "taxonomy.xml" #nom donné par défaut
    Phylo.write(tree_taxo_xml, name, "phyloxml")

    return tree_taxo_xml, name 


def pep_list_parser(file_name_pep):
    """Parse a list of peptides from a csv file

    Args:
        file_name_pep (str): the name of the csv file containing the peptides

    Raises:
        NameError: if the extension is not csv

    Returns:
        list: each sublist corresponds to a line with peptides from the csv file
    """
    _, ext = os.path.splitext(file_name_pep)
    if ext != ".csv": #on regarde si le fichier est bien un csv
        raise NameError("the input file must be a csv file")
    list_pep = []
    f= open(file_name_pep)
    myReader = csv.reader(f)
    for row in myReader:
        if row != []:
            tmp = True
            seq_pep = row
            for pep in seq_pep:
                for elt in pep:
                    if elt not in acides_amines: #on vérifie que les séquences des peptides sont bien des séquences protéiques
                        tmp = False
            if tmp :
                list_pep.append(seq_pep) #on ajoute la sous liste de peptides correspondant à une ligne dans le fichier dans la grande liste
    return list_pep


def unique_pep_finder(phyloxml_taxo, unique_pep_list):
    """Finds peptide sequences that are present in the proteic sequence of the species in the tree .

    Args:
        phyloxml_taxo (Phylogeny): the taxonomic tree in the phyloxml format
        unique_pep_list (list): the list containing sublists of peptides from a line in the csv file. 
                                The peptides that we want to check if there are caracteristic of a clade or not.

    Returns:
        dict: dictionnary mapping the clade's name and the list of the peptide found in their sequence (None if the peptide was not found)
    """
    dict_clade_pep_found = {}
    tree_taxo_phyloxml = Phylo.read(phyloxml_taxo, 'phyloxml') #on lit le fichier phyloxml 
    for liste in unique_pep_list: #pour chaque sous liste (= ligne du fichier csv)
        for clade in tree_taxo_phyloxml.get_terminals(): #pour toutes les espèces
            current_sequence = str(clade.sequences[0].mol_seq) #on reprend la séquence protéique associée à l'espèce qu'on regarde
            pep_present = [None if pep not in current_sequence else pep for pep in liste] #les peptides non présent dans la séquence sont remplacés par None
            if clade.name not in dict_clade_pep_found:
                dict_clade_pep_found[clade.name] = [pep_present] #on ajoute au dictionnaire le nom de l'espèce et la liste avec des sous listes (= ligne fichier csv) avec les peptides présents 
            else:
                dict_clade_pep_found[clade.name].append(pep_present)

    return dict_clade_pep_found


def unique_pep_in_tree(dict_clade_pep_found, tree_taxo_phyloxml, liste_pep):
    """Prints the species that have all peptides from a list (the list being a line in the csv file) in their sequence. Also finds the LCA of those species and says whether or 
    not one clade in the LCA doesn't share all the peptides from the list.
    If no species in the tree have all peptides from the list, it prints the species that contains at least one peptide from the list and find their LCA. Also says if one specie
    in the LCA doesn't have any peptide from the list in their sequence.

    Args:
        dict_clade_pep_found (dict): dictionnary mapping the clade's name and the list of the peptide found in their sequence (None if the peptide was not found)
        tree_taxo_phyloxml (Phylogeny): the taxonomic tree in the phyloxml format
        liste_pep (list): each sublist corresponds to a line with peptides from the csv file
    """
    df_pep = pd.DataFrame.from_dict(dict_clade_pep_found, orient="index") #on convertit en dataframe avec pandas
    for i in range(df_pep.shape[1]): #pour chaque colonne du dataframe (une colonne est une ligne du fichier csv)
        pep_complete = liste_pep[i] #la liste complète de peptide dans cette ligne du csv
        liste_clade_name_complete = []
        dict_clade_not_complete = {}
        for j in range(df_pep.shape[0]): #pour chaque ligne du csv (pour chaque espèce)
            liste_pep_present = []
            for pep in df_pep.iloc[j,i]:
                if pep != None: #si le peptide n'est pas None c'est qu'il est présent dans la séquence de l'espèce
                    liste_pep_present.append(pep)
            if len(liste_pep_present) == len(pep_complete): #si on a le même nb d'éléments, c'est que tout les peptides de la liste sont présents dans la séquence
                liste_clade_name_complete.append(df_pep.index[j]) #on ajoute le nom de l'espèce dans la liste des espèces avec tout les peptides 
            elif any(pep in pep_complete for pep in liste_pep_present): #sinon on vérifie si on retrouve au moins 1 peptide dans la liste
                liste_pep_not_found = [pep for pep in pep_complete if pep not in liste_pep_present] #liste des peptides non trouvés 
                if df_pep.index[j] not in dict_clade_not_complete:
                    dict_clade_not_complete[df_pep.index[j]] = liste_pep_not_found  #on ajoute au dictionnaire le nom de la clade et la liste des peptides non trouvés
                else:
                    dict_clade_not_complete[df_pep.index[j]].append(liste_pep_not_found)
        print("For the peptides {} :".format(pep_complete))
        if liste_clade_name_complete != []: #s'il y a des espèces avec tous les peptides de la liste
            print("Clades that have all the peptides :")
            print(liste_clade_name_complete)
            mrca = tree_taxo_phyloxml.common_ancestor(liste_clade_name_complete) #on regarde leur dernier ancêtre commun
            if len(liste_clade_name_complete) == 1:
                print("=========================================================")
            else:
                liste_LCA_leaf_not_complete = []
                liste_LCA_leaf_incomplete = []
                for leaf in get_leaf_nodes(mrca): #pour toutes les espèces se trouvant dans l'arbre du dernier ancêtre commun
                    if leaf.name not in liste_clade_name_complete: #si l'espèce ne fait pas partie des espèces avec tous les peptides dans leur séquence
                        if leaf.name in dict_clade_not_complete: #si l'espèce se trouve dans le dictionnaire des espèces avec au moins 1 peptide de la liste dans leur séquence
                            liste_LCA_leaf_not_complete.append(leaf.name) #on ajoute le nom de l'espèce dans une liste
                        else:
                            liste_LCA_leaf_incomplete.append(leaf.name) #sinon c'est que l'espèce n'a aucun peptide de la liste dans sa séquence et on ajoute le nom de l'espèce dans une autre liste
                if liste_LCA_leaf_not_complete != []: #si on a trouvé des espèces dans l'arbre du LCA avec au moins 1 peptide de la liste présent dans leur séquence
                    print("The LCA contains species that don't have all the peptides")
                    print("The LCA's clades that don't have all the peptides and the peptides missing are :")
                    for leaf_name in liste_LCA_leaf_not_complete:
                        print("{} : {}".format(leaf_name, dict_clade_not_complete[leaf_name]))
                if liste_LCA_leaf_incomplete != []: #si on a trouvé des espèces dans l'arbre du LCA qui n'avait aucun des peptides de la liste dans leur séquence
                    print("LCA contains clades that don't have any peptides from the list ")
                    print("The clades that don't have any peptide from the list are :")
                    print(liste_LCA_leaf_incomplete)
        elif dict_clade_not_complete != {}: #s'il n'y a aucune espèce qui a tous les peptides de la liste dans sa séquence
            print("No clades have all the peptide in the list")
            print("Clades that have at least one peptide of the list + the peptides not found :")
            print(dict_clade_not_complete) #on affiche celle qui ont au moins 1 peptide de la liste
            mrca = tree_taxo_phyloxml.common_ancestor(dict_clade_not_complete.keys()) #on trouve leur LCA

            liste_LCA_leaf_incomplete = []
            for leaf in get_leaf_nodes(mrca):
                if leaf.name not in dict_clade_not_complete: # on regarde dans le l'arbre du LCA s'il y a des espèces qui n'ont aucun peptides de la liste dans leur séquence
                    liste_LCA_leaf_incomplete.append(leaf.name)
            if liste_LCA_leaf_incomplete != []:
                print("LCA contains clades that don't have any peptides from the list ")
                print("The clades that don't have any peptide from the list are :")
                print(liste_LCA_leaf_incomplete) #on affiche les espèces dans l'arbre du LCA qui n'ont aucun des peptides de la liste dans leur séquence
        else: #si aucune séquence n'a aucun peptide de la liste
            print("The peptides are not found in any clade's sequence")
            mrca = None
        if mrca != None:
            if mrca.clades != []: 
                Phylo.draw_ascii(mrca)
                print("All the clades are parts of the {}".format(mrca.name))
                print("=========================================================")


##########################################################################################################################################################################
def peptide_signature_unique_pep_parser(peptide_signature_output, clade):
    """specific method to read output file from peptide_signature and have a csv file for test

    Args:
        peptide_signature_output (str): the output of peptides_signature
        clade (str): the specific taxononic level we used in peptides_signature

    Raises:
        NameError: if the input file is not a csv
    """
    root, ext = os.path.splitext(peptide_signature_output)
    if ext != ".csv":
        raise NameError("the input file must be a csv file")

    tmp =  True
    current_prot_name = ""
    dict_pep = {}
    df = pd.read_csv(peptide_signature_output, skipfooter=1, engine='python')
    for i in range(len(df)):
        seq = df.iloc[i,-1]
        for elt in seq:
            if elt not in acides_amines:
                tmp = False
        if tmp:
            if clade == "specie":
                col_index = 2
            elif clade == "genus":
                col_index = 1
            else:
                col_index = 0
            prot_name = df.iloc[i, col_index]
            if prot_name != current_prot_name:
                current_prot_name = prot_name
            if prot_name not in dict_pep:
                dict_pep[prot_name] = [seq]
            else:
                dict_pep[prot_name].append(seq)

    for key in dict_pep:
        dict_pep[key] = list(dict.fromkeys(dict_pep[key]))
    
    name = root.split("/")[2] + "_reformated.csv"
    #'peptide_unique_test.csv'
    with open(name, 'w') as csv_file:  
        writer = csv.writer(csv_file)
        for value in dict_pep.values():
            writer.writerow(value)
    return name
##############################################################################################################################################################################                


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--phylogeneticTree", required=True, help="The name of the phylogenetic tree in Newick", type=str)
    parser.add_argument("-t", "--taxonomicTree", required=True, help="The name of the taxonomic tree in Newick", type=str)
    parser.add_argument("-s", "--sequenceFasta", required=True, help="The name of the fasta file containing the proteic sequences used to build the phylogeny", type=str)
    parser.add_argument("-e", "--peptideSeq", required=True, help="The name of the cvs file containing the peptides", type=str)
    args = parser.parse_args()

    # L'arbre généré par un outils de phylogénie comme phyML par exemple avec une séquence servant de marqueur phylogénétique ici le collagène
    tree_phylo = Phylo.read(args.phylogeneticTree, "newick")
    tree_phylo.rooted = True
    # L'arbre généré via Lifemap en prenant les espèces trouvé à partir de Blast
    tree_taxo = Phylo.read(args.taxonomicTree, "newick")
    tree_taxo.rooted = True

    remove_not_common_clades(tree_taxo, tree_phylo) #on enlève les espèces présentes dans l'arbre taxonomique qui ne sont pas présentes dans la phylogénie

    seq_dict = read_multi_fasta(args.sequenceFasta)

    tree_taxo_xml, tree_taxo_xml_name = newick_to_xml_with_seq(seq_dict, tree_taxo) 
    
    dico_test_taxo = get_dico_clade_non_terminal(tree_taxo) #On génère un dictionnaire qui match les clades avec leur espèces dans la taxonomie à partir de l'arbre taxonomique

    #On génère un dictionnaire mappant les clades avec les espèces trouvées dans l'arbre phylogénétique (correspondant au dernier ancêtre commun dans l'arbre phylogénétique)
    dico_mrca_phylo = get_clade_mrca_phylo(tree_phylo, dico_test_taxo) 
    
    #ici on génère un nouveau dictionnaire en reprenant le dictionnaire de get_clade_mrca_phylo() mais cette fois ci en groupant les espèces appartenant à un même sous arbre fils. 
    # Le sous arbre correspondant au dernier ancêtre commun dans l'arbre phylogénétique a 2 sous arbres fils (l'arbre phylogénétique est binaire) et 
    # on regroupe les espèces en fonction de leur appartenance à un sous arbre.
    species_grouped = get_species_grouped(dico_test_taxo, tree_phylo) 
   
    # La fonction qui permet de générer le print des clades problèmatiques en comparant le sous arbre de la clade dans l'arbre taxonomique et le sous arbre du dernier ancêtre commun dans l'arbre phylogénétique. 
    # print de l'arbre correspondant au dernier ancêtre commun + les plus petits sous arbres cohérants avec la taxonomie ainsi que leur distance avec la racine du dernier ancêtre commun
    # Print également les espèces mal placées s'il y en a pour chaque clades.
    # Donne également le nombre d'espèces de la clade dans la taxonomie et nb d'espèce dans l'arbre correspondant au dernier ancêtre commun de ces espèces dans l'arbre phylogénétique
    count_misplaced_clade(dico_test_taxo, dico_mrca_phylo, species_grouped, tree_phylo)

    #csv_reformated = peptide_signature_unique_pep_parser("peptides_signature/results_taxonomie_complete/unique_pep_family.csv", "family")

    #liste_pep_unique = pep_list_parser(csv_reformated)
    liste_pep_unique = pep_list_parser(args.peptideSeq)
   
    dict_clade_pep_found = unique_pep_finder(tree_taxo_xml_name, liste_pep_unique)
    unique_pep_in_tree(dict_clade_pep_found, tree_taxo_xml, liste_pep_unique)


    