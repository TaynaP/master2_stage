#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import cgitb
import cgi
import json
import re
import sys
import csv
import common
import molmass
import peptides

#Vérification du mail
def verif_mail(mail):
    """ check the validity of email address """

    res = mail.strip()
    if res.find(' ') > 0:
        return False
    a = res.find('@')
    if a <= 0:
        return False
    if res.find('@', a+1) > 0:
        return False
    return True

#Vérifie si les séquences sont au bon format
def valid_sequence(seq):
    """ check the validity of the sequence seq """

    if seq is None:
        return False

    pattern =  """
        ^                               # beginning of string
        (                               # one sequence 
        \ *>.+[\r\n]+                   # line of the name
        ([-\.\ ATCGU0-9]+[\r\n]+)+      # line of the primary structure
        )+$                             # end of one sequence                  
        """
    if re.search(pattern, seq, re.VERBOSE) is None:
        return False
    else:
        return True

def valid_proteins(prot):
    """ check the validity of the proteic sequence prot """

    if prot is None:
        return False

    pattern =  """^(\ *>.+[\r\n]+([-\.\ GPAVLIMCFYWHKRQNEDST0-9]+[\r\n]+)+)+$"""
    if re.search(pattern, prot, re.VERBOSE) is None:
        return False
    else:
        return True


#vérification si les peptides sont au bon format
def valid_peptide(pep):
    """check the validity of the peptide"""

    if pep is None:
        return False
    
    pattern = "^(([a-zA-Z]+)=(\[([GPAVLIMCFYWHKRQNEDST,\ ]+)\])[\ ,\r\n]+)+$" #regex pattern pour le format des peptides 

    if re.search(pattern, pep, re.VERBOSE) is None:
        return False
    else:
        return True

#si les peptides sont formatés correctement, parse pour avoir le nom et les peptides
def format_peptide(pep):
    dict_pep = {}
    try:
        list_pep_input = re.split(r"[\[,\]]", pep)
        for i in range(len(list_pep_input)):
            if re.search(r"=", list_pep_input[i]):
                pep_name = list_pep_input[i].rstrip(list_pep_input[i][-1])
                list_pep = []
                j = i+1
                while all((list_pep_input[j] != None, list_pep_input[j] != "")):
                    list_pep.append(list_pep_input[j].strip())
                    j += 1 
                dict_pep[pep_name] = list_pep
        return dict_pep
    except:
        return None
    

#Vérification si le format donné est bien fasta
def format_fasta(seq):
    try:
        list_seq = seq.splitlines()
        temp_seq = ''
    
        # remove spaces for non-title line
        for x in list_seq:
            if re.search('^ *>', x) is None:
                x = x.upper()
                x = x.strip()
            else:
                x = x.replace(' ', '_')
            temp_seq = temp_seq + x + '\n'

        return temp_seq
    except:
        return None


#Création du fichier csv avec les peptides
def create_csv_pep(run_id, pep, check_mass, resolution):
    """ 
    create a csv file with pep
    """

    csv_pep_path = common.csv_pep_path_result(run_id)
    csv_pep = open(csv_pep_path, "w")
    pep_dict = pep

    list_pep = []
    for label in pep_dict.keys():
        for peptide in pep_dict[label]:
            if check_mass != None:
                f = molmass.Formula('peptide({})'.format(peptide))
                if resolution == "0.1":
                    resolution = ".1f"
                elif resolution == "0.001":
                    resolution = ".3f"
                mass = [format(f.isotope.mass, resolution)]
            else:
                mass = []

            node = []
            if list_pep != []:
                list_pep_seq = [pep.get_seq() for pep in list_pep]
                if peptide not in list_pep_seq:
                    list_pep.append(peptides.Peptides(peptide, mass, node, label)) #pour ne pas avoir des peptides en double
            else:
                list_pep.append(peptides.Peptides(peptide, mass, node, label))


    writer = csv.writer(csv_pep)
    writer.writerow(["Sequence", "Mass", "Specie", "Position", "PTM", "label", "Comment"])
    for pep in list_pep:
        writer.writerow([pep.get_seq(), pep.get_mass(), pep.get_nodes(), None, None, pep.get_label(), None])


#Création du fichier fasta avec les séquences soumises
def create_fasta(run_id, seq, type):
    """ 
    create a fasta file with the sequences seq 
    """

    if type != None:
        fasta_path = common.fasta_path_result(run_id, type)
    ffasta = open(fasta_path, 'w')
    seq_list = seq.split ('>')
    nb_seq = len(seq_list)
    name_seq = []
    for i in range(1, len(seq_list)):
        split_seq = seq_list[i].split('\n')
        name_seq.append(split_seq[0])
        ffasta.write('>' + split_seq[0] + '\n')
        for x in range(1, len(split_seq)):
            ch = re.sub('[ 0-9-]', '', split_seq[x])
            ffasta.write(ch)
    ffasta.write('\n')        
    ffasta.close()


#Fonction permettant de vérifier les données soumises dans le formulaire, retour dans un dictionnaire (req)
#Fonction qui sera modifiée pour le traitement de vos données

def extract_request(form):
    """ check the validity of the different fields of form """

    error_found = False
    error_messages = []
    req = {}
    
    # check peptides
    # if peptides are given in the textarea "peptides"
    # rajouter le fait de pouvoir donner un fichier en entrée
    if (form.has_key('peptides') and not (len(form['peptides'].value.strip()) is 0) ):
        pep = form['peptides'].value
    else:
        pep = None
        error_found = True
        error_messages.append('Missing sequences.')
    if pep is not None:
        # valid the peptide format
        pep_test = valid_peptide(pep)
        pep = format_peptide(pep)

        if not pep_test:
            error_found = True
            error_messages.append("Sequences are not in expected format")
        else:
            req["pep"] = pep

    ## pour le nom des séquences
    if form.has_key('seq_name') and len(form['seq_name'].value) is not 0:
        req["seq_name"] = form['seq_name'].value
    else:
        req["seq_name"] = None

    #check proteins
    # if both "proteins" and "file_prot" are filled => error
    if (form.has_key('proteins') and form.has_key('name_file')
        and not((len(form['proteins'].value.strip()) is 0) or 
                (len(form['file_prot'].value) is 0))
        ):
        error_found = True
        error_messages.append('Paste sequences <strong>OR</strong>' +
                              'upload a file.')
       
    else:
        # if a file is uploaded
        if form.has_key('file_prot') and not (len(form['file_prot'].value) is 0):
            seq_prot = form['file_prot'].value
            
        # if sequences are given in the textarea "proteins"
        elif (form.has_key('proteins') and
              not (len(form['proteins'].value.strip()) is 0) ):
            seq_prot = form['proteins'].value
        else:
            seq_prot = None
            
        if seq_prot is not None:

            # valid the sequence format
            seq_prot = format_fasta(seq_prot)
            seq_prot_test = valid_proteins(seq_prot)
 
            if not seq_prot_test:
                error_found = True
                error_messages.append("Sequences are not in Fasta format")
            else:
                req["seq_prot"] = seq_prot
    
    #check rnas
    # if both "rna" and "file_rna" are filled => error
    if (form.has_key('rna') and form.has_key('name_file')
        and not((len(form['rna'].value.strip()) is 0) or 
                (len(form['file_rna'].value) is 0))
        ):
        error_found = True
        error_messages.append('Paste sequences <strong>OR</strong>' +
                              'upload a file.')
       
    else:
        # if a file is uploaded
        if form.has_key('file_rna') and not (len(form['file_rna'].value) is 0):
            seq_rna = form['file_rna'].value
            
        # if sequences are given in the textarea "rna"
        elif (form.has_key('rna') and
              not (len(form['rna'].value.strip()) is 0) ):
            seq_rna = form['rna'].value
        else:
            seq_rna = None
            
        if seq_rna is not None:

            # valid the sequence format
            seq_rna = format_fasta(seq_rna)
            seq_rna_test = valid_sequence(seq_rna)
 
            if not seq_rna_test:
                error_found = True
                error_messages.append("Sequences are not in Fasta format")
            else:
                req["seq_rna"] = seq_rna

    #check dnas
    # if both "dna" and "file_dna" are filled => error
    if (form.has_key('dna') and form.has_key('name_file')
        and not((len(form['dna'].value.strip()) is 0) or 
                (len(form['file_dna'].value) is 0))
        ):
        error_found = True
        error_messages.append('Paste sequences <strong>OR</strong>' +
                              'upload a file.')
       
    else:
        # if a file is uploaded
        if form.has_key('file_dna') and not (len(form['file_dna'].value) is 0):
            seq_dna = form['file_dna'].value
            
        # if sequences are given in the textarea "dna"
        elif (form.has_key('dna') and
              not (len(form['dna'].value.strip()) is 0) ):
            seq_dna = form['dna'].value
        else:
            seq_dna = None
            
        if seq_dna is not None:

            # valid the sequence format
            seq_dna = format_fasta(seq_dna)
            seq_dna_test = valid_sequence(seq_dna)
 
            if not seq_dna_test:
                error_found = True
                error_messages.append("Sequences are not in Fasta format")
            else:
                req["seq_dna"] = seq_dna
    
    if seq_prot == None and seq_rna == None and seq_dna == None:
        error_found = True
        error_messages.append('Missing sequences.')
 
    #check radio button model
    if form.has_key('model'):
        req["model"] = form['model'].value
    else:
        req["model"] = None
    
    #if a taxonomy file is uploaded
    if form.has_key('file_taxonomy') and not (len(form['file_taxonomy'].value) is 0):
        taxonomy = form['file_taxonomy'].value
        req["taxonomy"] = taxonomy
    
    #check radio button digestion
    if form.has_key('digestion'):
        req["digestion"] = form['digestion'].value
    else:
        req["digestion"] = None

    #check choose of enzyme
    if form.has_key('enzyme'):
        req["enzyme"] = form['enzyme'].value
    else:
        req["enzyme"] = None
    
    #check mass_comparison
    if form.has_key('check_mass'):
        req["check_mass"] = form.has_key("check_mass")

        if form.has_key('resolution'):
            if form['resolution'].value == "0.1":
                req["resolution"] = "1"
            elif form['resolution'].value == "0.001":
                req["resolution"] = "3"
        else:
            error_messages.append("Choose a resolution")

    else:
        req["check_mass"] = None

   # check email value
    if not form.has_key('email'):
        req["email"] = None
    else:
        email = form['email'].value.strip()
        if email is not '':
            if not verif_mail(email):
                error_found = True
                error_messages.append('Invalid <strong>email ' +
                                      'address</strong>.')
            else:
                req["email"] = email
        else:
            req["email"] = None
                            
    req["error_messages"] = error_messages
    
    return req    

#Lance le programme en fonction des données contenant dans req
#Fonction à modifier pour l'adapter à votre programme
def launch_software(run_id, req, model):
    #-e None and -m None (proteins)
    if req["seq_prot"] != None and req["pep"] != None and model != None and req["digestion"] == None and req["enzyme"] == None and req["check_mass"] == None:
        os.system("main_add_species_to_peptide.py -type %s -fpep %s -fprot %s -output %s -model %s"%("proteins", common.csv_pep_path_result(run_id), common.fasta_path_result(run_id, "proteins"), common.result_dir(run_id), model))
    #-e not None and -m None  (proteins)
    elif req["seq_prot"] != None and req["digestion"] != None and req["enzyme"] != None and req["check_mass"] == None:
        os.system("main_add_species_to_peptide.py -type %s -fprot %s -output %s -model %s -e %"%("proteins", common.fasta_path_result(run_id, "proteins"), common.result_dir(run_id), model, req["enzyme"]))
    #-e not None and -m not None (proteins)
    elif req["seq_prot"] != None and req["digestion"] != None and req["enzyme"] != None and req["check_mass"] != None:
        os.system("main_add_species_to_peptide.py -type %s -fprot %s -output %s -model %s -e % -m %"%("proteins", common.fasta_path_result(run_id, "proteins"), common.result_dir(run_id), model, req["enzyme"], req["resolution"]))

    # -e None and -m None (rna)
    elif req["seq_rna"] != None and req["pep"] != None and model != None and req["digestion"] == None and req["enzyme"] == None and req["check_mass"] == None:
        os.system("main_add_species_to_peptide.py -type %s -fpep %s -fprot %s -output %s -model %s"%("proteins", common.csv_pep_path_result(run_id), common.fasta_path_result(run_id, "RNA"), common.result_dir(run_id), model))
    #-e not None and -m None (rna)
    elif req["seq_rna"] != None and req["digestion"] != None and req["enzyme"] != None and req["check_mass"] == None:
        os.system("main_add_species_to_peptide.py -type %s -fprot %s -output %s -model %s -e %"%("proteins", common.fasta_path_result(run_id, "RNA"), common.result_dir(run_id), model, req["enzyme"]))
    #-e not None and -m not None (rna)
    elif req["seq_rna"] != None and req["digestion"] != None and req["enzyme"] != None and req["check_mass"] != None:
        os.system("main_add_species_to_peptide.py -type %s -fprot %s -output %s -model %s -e % -m %"%("proteins", common.fasta_path_result(run_id, "RNA"), common.result_dir(run_id), model, req["enzyme"], req["resolution"]))

#Vérifie si il n'y a pas d'erreur au retour des données soumises par le formulaire
def process_request(form):
    run_id = form['run_id'].value
    error=0
    error_page=""
    req = extract_request(form)
    
    if len(req["error_messages"])>0:
        error_page=req["error_messages"]
        error=1
    else :
        if req["pep"] != None:
            create_csv_pep(run_id, req["pep"], req["check_mass"], req["resolution"])

        if req["seq_prot"] != None:
            create_fasta(run_id, req['seq_prot'], "proteins")

        if req["seq_rna"] != None:
            create_fasta(run_id, req['seq_rna'], "RNA")

        if req["seq_dna"] != None:
            create_fasta(run_id, req['seq_dna'], "DNA")

        #séquence modèle pour la protéine mature
        if req["model"] == "col1a2": 
            model = "files/model_human_col1a2.fasta"  
        elif req["model"] == "col1a1":
            model = "files/model_human_col1a1.fasta"
        else:
            model = None
        #rajouter la possibilité d'entrer une séquence modèle 

        #lance le programme
        launch_software(run_id, req, model)

    return (error, error_page)



def main():
    fs = cgi.FieldStorage()
    
    error, error_page=process_request(fs)

    sys.stdout.write("Content-Type: application/json")

    sys.stdout.write("\n")
    sys.stdout.write("\n")

    result={}
    if error==0:
        result['success'] = True
        result['run_id'] = fs['run_id'].value
        result['path_result'] = common.HTML_PATH+"/result/"+fs['run_id'].value+"/results.php"
    else:
        result['success'] = False

    sys.stdout.write(json.dumps(result, indent=1))
    sys.stdout.write("\n")

main()
