import argparse
import subprocess

from peptide_marker_validation import *

from utils import *












if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Search peptides in different proteic sequences")
    parser.add_argument("-type", "--TypeOfSeq", required=True, help="Type of sequences used, 3 choices : DNA, RNA or proteins", type=str)
    parser.add_argument("-fpep","--NamePep", required=True, help="Input CSV file with species known or not (not required if using -e)", type=str)
    parser.add_argument("-fprot", "--NameSeq", required=True, help="The name of the fasta file containing the sequences (required)", type=str)
    parser.add_argument("-output", "--NameOutput", required=True, help="The name of the output file (required)", type=str)
    parser.add_argument("-model", "--NameModel", required=True, help="The name of the model sequence (required)", type=str)
    parser.add_argument(
        "-e", "--Enzyme", required=False, 
        help="Use this option if you want to only see peptides from a digestion with a given enzyme. It's the enzyme number in RPG (optional).", 
        type=int
    )
    parser.add_argument(
        "-m", "--MassThreshold", required=False, 
        help="Use this if you want to check the mass equality between peptides in addition to the sequence equality. It's the number of decimal digits to check for mass equality (optional).",
        type=int
    )
    parser.add_argument("-PTM", "--MassPTM", required=False, help="Use this if you want to take into account PTM mass. Use only with -m")
    args = parser.parse_args()

    if args.MassPTM != None and args.MassThreshold == None:
        raise Exception("Use -PTM only with -m")

    if args.NamePep == None and args.Enzyme == None:
        raise Exception("Use -fpep if you're not using -e")
    
    if args.TypeOfSeq not in ["DNA", "RNA", "proteins"]:
        raise Exception("please type one of 3 : DNA, RNA or proteins")

    if args.TypeOfSeq == "DNA" and args.MassThreshold != None:
        raise Exception("We can't do mass search for DNA sequence")
    
    if not os.path.exists(args.NameOutput):
        os.makedirs(args.NameOutput)
    else:
        if os.listdir(args.NameOutput):
            raise Exception("This directory already exists. Please give another directory.")
    
    dict_model = read_multi_fasta(args.NameModel) #le mod??le utilis?? pour avoir la prot??ine mature
    
    if args.TypeOfSeq == "proteins": #si on a des s??quences prot??iques
        dict_prot_immature = read_multi_fasta(args.NameSeq) #les s??quences fasta immatures
        dict_prot = prot_immature_to_mature(dict_prot_immature, list(dict_model.items())[0][1][0]) #fonction qui permet d'avoir les prot??ines matures

    elif args.TypeOfSeq == "RNA": #si on a des s??quences d'ARN
        dict_rna = read_multi_fasta(args.NameSeq) #parse les s??quences d'ARN
        dict_prot_immature = rna_ro_prot(dict_rna) #les s??quences d'ARN traduites en prot??ines immatures
        dict_prot = prot_immature_to_mature(dict_prot_immature, list(dict_model.items())[0][1][0]) #fonction qui permet d'avoir les prot??ines matures

    else: #si on utilise des s??quences d'ADN
        dict_dna = read_multi_fasta(args.NameSeq) #parse les s??quences d'ADN
        dict_prot = dna_to_prot(dict_dna) #les s??quences d'ADN traduites en prot??ines dans les 3 cadres de lecture

    #When -e and -m are omitted
    #Le cas o?? on donne des peptides en entr??e et on ne regarde pas la masse
    if args.Enzyme == None and args.MassThreshold == None:
        list_pep = pep_parser_csv(args.NamePep) #parser pour les peptides (format csv)
        if args.TypeOfSeq == "proteins" or args.TypeOfSeq == "RNA": #si le type est ARN ou prot??ines
            list_pep_with_species = search_pep_in_seq(list_pep, dict_prot) #fonction pour rechercher les peptides dans les s??quences (sans ??galit?? de masse)
            name_file = args.NameOutput + "peptides.csv"
            pretty_print_peptides_csv(list_pep_with_species, name_file) #fonction pour g??n??rer le fichier csv r??sultat
        else:
            pass
            #Le cas de l'ADN ?? finir pour chercher les peptides dans les s??quences


    #When only -m is omitted
    #Le cas o?? on g??n??re les peptides par rpg et on ne regarde pas la masse
    elif args.MassThreshold == None:
        #Use of RPG
        csv_file_rpg = args.NameOutput + "rpg.csv"
        create_fasta(args.NameOutput, dict_prot) #on cr??e un fichier fasta avec les s??quences prot??iques matures (dict_prot contient les s??quences prot??iques matures)
        mature_prot_fasta = args.NameOutput + "mature_prot.fasta"
        subprocess.run(
            ["rpg", "-i", mature_prot_fasta, "-o", csv_file_rpg, "-f", "csv", "-e", str(args.Enzyme)] #on lance rpg
        )
        pretty_print_peptides_csv_from_rpg(csv_file_rpg, args.NameOutput) #converts RPG file to CSV format accepted
        csv_file_rpg_converted = args.NameOutput + "rpg_converted.csv"
        list_pep_rpg = pep_parser_csv(csv_file_rpg_converted) #on parse les peptides g??n??r??s par rpg
        list_pep_input = pep_parser_csv(args.NamePep) #parser pour les peptides donn??s (format csv)
        dict_pep_input_equal_pep_rpg = compare_pep_rpg_pep_input(list_pep_rpg, list_pep_input, dict_prot) #fonction pour comparer les peptides inputs et les peptides rpg
        name_file = args.NameOutput + "peptides.csv"
        pretty_print_pep_rpg_equal_pep_input(dict_pep_input_equal_pep_rpg, name_file) #fonction pour g??n??rer le fichier csv r??sultat
    
    #When only -e is omitted
    #When -m is used the masses are replaced
    elif args.Enzyme == None:
        threshold = ".{}f".format(args.MassThreshold)
        list_pep_with_mass = pep_parser_csv(args.NamePep, threshold)
        dict_pep_with_species_mass, dict_pep_with_species_mass_only = search_pep_in_seq_with_mass(list_pep_with_mass, dict_prot, threshold) #too long
        #pas fini
    
    #Le cas o?? utilise -e et -m 
    #on g??n??re les peptides par rpg et on regarde la masse
    else:
        threshold = ".{}f".format(args.MassThreshold) #la pr??cision soit 1 ou 3 chiffres apr??s la virgule
        csv_file_rpg = args.NameOutput + "rpg.csv"
        create_fasta(args.NameOutput, dict_prot) #on cr??e un fichier fasta avec les s??quences prot??iques matures (dict_prot contient les s??quences prot??iques matures)
        mature_prot_fasta = args.NameOutput + "mature_prot.fasta"
        #Use of RPG
        subprocess.run(
            ["rpg", "-i", mature_prot_fasta, "-o", csv_file_rpg, "-f", "csv", "-e", str(args.Enzyme)] 
        )
        pretty_print_peptides_csv_from_rpg(csv_file_rpg, args.NameOutput) #converts RPG file to CSV format accepted in README
        csv_file_rpg_converted = args.NameOutput + "rpg_converted.csv"
        list_pep_rpg = pep_parser_csv(csv_file_rpg_converted, threshold=threshold) #on parse les peptides g??n??r??s par rpg mais cette fois en rempla??ant la masse par la masse calcul??e par molmass
        list_pep_input = pep_parser_csv(args.NamePep, threshold=threshold) #parser pour les peptides donn??s (format csv) en rempla??ant la masse
        dict_pep_input_equal_pep_rpg = compare_pep_rpg_pep_input(list_pep_rpg, list_pep_input, dict_prot, threshold=threshold) #on compare les peptides input et les peptides rpg avec la masse
        name_file = args.NameOutput + "peptides.csv"
        pretty_print_pep_rpg_equal_pep_input(dict_pep_input_equal_pep_rpg, name_file) #fonction pour g??n??rer le fichier csv r??sultat