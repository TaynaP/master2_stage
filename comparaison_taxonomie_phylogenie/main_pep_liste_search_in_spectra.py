import argparse
import os 

from pep_parser import *
import sys














if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--peptidesCSVFile", required=True, help="Input CSV file created by rpg (required)", type=str)
    parser.add_argument("-o", "--resultsPath", required=True, help="The directory (path) wich will contain all the result files (required)", type=str)
    parser.add_argument("-t", "--taxonomicTree", required=True, help="The name of the taxonomic tree in Newick (required)", type=str)
    parser.add_argument("-s", "--sequenceFasta", required=True, help="The name of the fasta file containing the proteic sequences (required)", type=str)
    parser.add_argument("-p", "--peakMS", required=True, help="The name of the csv file with the peaks of MS (required)", type=str)
    parser.add_argument("-r", "--ptm", required=False, help="The name of the file given by PTM_mass_simulator", type=str)
    parser.add_argument("-m", "--threshold", help="the threshold decided to find equivalent peptide mass in the database (default '.3f')", type=str, default='.3f')
    args = parser.parse_args()


    if not os.path.exists(args.resultsPath):
        os.makedirs(args.resultsPath)
    else:
        if os.listdir(args.resultsPath):
            raise Exception("This directory already exists. Please give another directory.")
            

    tree_taxo = Phylo.read(args.taxonomicTree, "newick") #taxonomie au format Newick
    tree_taxo.rooted = True

    seq_dict = read_multi_fasta(args.sequenceFasta) #parser des séquences protéiques au format multifasta
  
    tree_taxo_xml, tree_taxo_xml_path = newick_to_xml_with_seq(seq_dict, tree_taxo, args.resultsPath) #converti la taxonomie du format newick au format xml
    
    # Create a list of Peptide objects
    pep = pep_parser_rpg(args.peptidesCSVFile, args.threshold)
    
    list_pep = which_species_pep_present(pep, tree_taxo_xml_path) #on trouve dans quelles espèces les peptides apparaissent
    
    list_pep_modified = get_smallest_good_tree(tree_taxo_xml_path, list_pep) #on remonte dans la taxonomie pour avoir les LCA des espèces dans lesquelles les peptides apparaissent
    
    if args.ptm: #si on utilise PTM_mass_simulator pour avoir les PTM 
        list_pep_ptm = PTM_csv_parser(args.ptm, list_pep_modified, args.threshold) #parser pour le fichier donné par PTM_mass_simulator (on s'attend à avoir les mêmes peptides)
        pretty_print_peptides_csv(list_pep_ptm, args.resultsPath) #génération du fichier résultat csv
        pretty_print_unique_peptides_fasta(list_pep_ptm, args.resultsPath) #génération du fichier résultat fasta 
    else: #si on utilise pas de PTM
        pretty_print_peptides_csv(list_pep_modified, args.resultsPath)
        pretty_print_unique_peptides_fasta(list_pep_modified, args.resultsPath)
    

    list_peak = main_peak_parser(args.peakMS, args.threshold) #on parse les fichiers correspondant au spectre de masse
    
    if args.ptm:
        dict_peak_found, peak_not_found = list_pep_compare_ptm(list_peak, list_pep_ptm) #on compare la masse des pics et des peptides avec PTM
    else:
        dict_peak_found, peak_not_found = list_pep_compare_ptm(list_peak, list_pep_modified) #on compare la masse des pics et des peptides sans PTM
    
    tree_taxo_xml_peak, tree_taxo_xml_peak_name = peak_to_clade(dict_peak_found, tree_taxo_xml_path, args.resultsPath) # on modifie l'arbre taxonomique pour associer les pics avec les noeud de l'arbre

    terminal_nodes = terminal_node_with_peak(tree_taxo_xml_peak) #on trouve les noeuds terminaux

    pretty_print_term_nodes_csv(terminal_nodes, args.resultsPath, tree_taxo_xml_peak) #génère le fichier résultat pour les noeuds terminaux dans l'arbre taxonomique

    with open(os.path.join(args.resultsPath, 'peak_found.txt'), 'w') as f:
        sys.stdout = f
        print(dict_peak_found) #génère un fichier pour pour les pics associés à un peptide
    
    with open(os.path.join(args.resultsPath, 'peak_orphan.txt'), 'w') as file:
        sys.stdout = file
        print(peak_not_found) #génère un fichier pour les pics orphelins
 